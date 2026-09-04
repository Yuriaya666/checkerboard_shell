# Checkerboard Shell

## 1. Overview

这是一个面向后续开发的最小 Checkerboard Shell forward simulation framework。它以位置为唯一壳体状态（position-only），在规则结构化网格上计算薄壳的 B/W Checkerboard 几何、能量、残量和非线性静态变形。

当前工程的输入是彼此独立的：

```text
geometry + material + boundary conditions + external load
                               ↓
                       CheckerboardShell.solve()
                               ↓
                    deformation / displacement
```

当前已经实现：

- 规则矩形 Checkerboard 网格和 quad connectivity；
- B 面/W 点的 midpoint 几何、第一基本形式和第二基本形式；
- 有限边界的 moment-exact W 闭合；
- 论文冻结的 `X_M` 膜补全和参考相对 `X_W_rel` 弯曲兼容项；
- 各向同性 StVK/Koiter 能量分块；
- 自动微分 residual 和 consistent tangent；
- Hessian 特征分解、下降方向和线搜索 Newton 静态求解；
- 节点死载荷、固定自由度、变形网格绘图和结果保存。

这里的目标是提供可复用的 forward engine，而不是论文复现仓库。当前没有实现 inverse design、动力学、接触、电磁载荷、专用压力载荷类、载荷延续驱动或非结构化网格求解器。

## 2. Installation

### 环境要求

- Python 3.10 或更高版本；
- NumPy；
- PyTorch；
- Matplotlib（只在调用绘图时使用）。

依赖在 `pyproject.toml`（工程安装与依赖配置文件）中声明。推荐在项目根目录执行 editable install（可编辑安装）：

```bash
python -m pip install -e .
```

然后可以运行 demo（示例脚本）：

```bash
python examples/plate_demo.py
```

也可以直接从源码 checkout 运行 `examples/plate_demo.py`（矩形板示例脚本）；脚本会自行把 `src/`（源码目录）加入 Python 路径。

运行测试：

```bash
python -m pytest -q
```

## 3. Project Structure

```text
checkerboard_shell/
├── src/checkerboard_shell/
│   ├── geometry/
│   ├── mechanics/
│   ├── material/
│   ├── solver/
│   └── visualization/
├── examples/
├── tests/
├── docs/
├── README.md
└── pyproject.toml
```

### `src/checkerboard_shell/`（可安装的 Python 包）

顶层 `__init__.py`（公共 API 导出文件）导出最常用的 `CheckerboardMesh`、`CheckerboardShell`、`IsotropicMaterial`、`NodalLoad`、`BoundaryConditions` 和 `NewtonSettings`。因此普通使用者不必从内部文件逐层导入。

### `geometry/`（网格与离散几何）

- `geometry/mesh.py`（结构化网格表示）：定义 `CheckerboardMesh`。`CheckerboardMesh.rectangle()`（规则矩形网格构造函数）生成 `[ny + 1, nx + 1, 3]` 的 `reference_vertices`（参考节点坐标）、`hx`/`hy`（网格步长）、row-major `quads`（quad 拓扑）以及 `bottom`、`top`、`left`、`right` 四条边界节点编号。
- `geometry/core.py`（B/W 基础几何）：`torch_geometry()`（B/W 几何计算函数）生成 B 面对角切向、B 法向、第一基本形式、W 点平均切向、未对称/对称第二基本形式、skew diagnostic 和 `chi`。`edge_midpoints()`（物理边中点构造函数）提供水平/垂直边中点。
- `geometry/boundary.py`（有限域边界闭合）：`boundary_geometry()`（含边界的 W 几何计算函数）在内部使用中心样本，在边界和角点使用 moment stencils；`targets_from_reference()`（参考状态缓存函数）和 `koiter_targets()`（能量目标与重构算子缓存函数）建立参考 metric、面积权重、`gamma_w` 和边界算子；`fast_boundary_geometry()`（向量化边界几何函数）供能量评估使用；`transported_boundary_x()`（参考 frame 传输函数）计算参考相对 `X_W_rel`。
- `geometry/intrinsic.py`（内禀兼容量）：`intrinsic_connection_field()`（B-frame connection 场）和 `intrinsic_m_coordinate()`（归一化 `X_M` 坐标）提供膜补全项。

当前几何算子是结构化 Checkerboard 算子，不是通用 unstructured mesh API。自定义几何必须保持规则的 `[ny + 1, nx + 1, 3]` 张量布局，并且边界闭合要求 `nx >= 3`、`ny >= 3`。

### `mechanics/`（能量和外载荷）

- `mechanics/energy.py`（论文冻结能量）：`make_targets()`（建立参考目标）返回 `EnergyTargets`；`paper_components()`（能量分块评估函数）组合 `membrane_b`、`membrane_xm`、`bending` 和 `bending_compatibility`，同时返回 `X_M`、`X_W_rel` 和几何诊断；`potential_from_positions()`（底层全位置势能函数）可用节点力张量计算内能减去死载荷势，普通用户通常通过 `CheckerboardShell.potential()` 调用。
- `mechanics/loads.py`（节点外载荷）：`NodalLoad` 保存与 mesh 同形状的 `[ny + 1, nx + 1, 3]` 节点力。`NodalLoad.zero()`（零载荷构造函数）建立无载荷；`NodalLoad.total_force()`（节点合力均分函数）把一个 resultant force 均匀分配到给定节点；`scaled()`（载荷缩放函数）用于载荷因子。
- `mechanics/model.py`（forward 模型组装）：`CheckerboardShell` 将 mesh、material、边界条件、load 和 Newton settings 连接起来，提供 `positions()`、`energy_components()`、`potential()`、`residual()`、`tangent()` 和 `solve()`。

### `material/`（材料参数）

- `material/isotropic.py`（各向同性材料定义）：`IsotropicMaterial` 保存 `young`（Young 模量）、`poisson`（Poisson 比）和 `thickness`（厚度），并提供 `lame_alpha`、`lame_mu`；`stvk_density()`（StVK 密度函数）实现能量中的各向同性本构收缩。

### `solver/`（约束、残量和非线性求解）

- `solver/boundary_conditions.py`（固定自由度）：`BoundaryConditions.fixed_nodes()`（固定节点构造函数）支持按节点和 `(x, y, z)` 分量固定；`fixed_boundary()`（只固定一条边界线）和 `clamped_boundary()`（固定边界及相邻节点线）提供结构化网格快捷方式。
- `solver/newton.py`（Newton 求解器）：`NewtonSettings` 是 `SolverSettings` 的公共别名；`value_gradient_hessian()`（自动微分梯度/Hessian 函数）和 `solve_load_level()`（单一载荷水平求解函数）来自当前论文非线性求解实现；`NewtonSolver`（对象封装）供 `CheckerboardShell.solve()` 调用。求解过程使用 Hessian 特征分解、下降方向、Armijo/残量线搜索，不改变能量表达式。
- `solver/result.py`（求解结果）：`SimulationResult` 保存参考节点、变形节点、自由坐标、能量、迭代历史、状态和残量，并提供 `displacement`、`max_displacement`、`plot()` 和 `save()`。

### `visualization/`（可视化）

`visualization/plotting.py`（变形网格绘图模块）的 `plot_deformation()`（三维参考/变形 quad 网格绘图函数）只负责显示和保存图像，不参与几何、能量或 Newton 计算。绘图依赖 Matplotlib，并在调用时延迟导入。

### `examples/`、`tests/` 和 `docs/`

- `examples/plate_demo.py`（最小矩形板 demo）展示一条完整 forward 调用链，并生成图像和 NumPy 结果文件。
- `tests/test_forward.py`（forward 回归测试）检查参考状态零能量、零载荷收敛、非零载荷变形，以及 residual/tangent 的形状、有限性和对称性。
- `docs/CORE_PIPELINE.md`（核心依赖链文档）记录抽取模块与 B/W、`X_M`、`X_W_rel`、能量、Newton 阶段的对应关系。

## 4. Quick Start

最简单的运行方式是：

```bash
python examples/plate_demo.py
```

该 demo（示例脚本）执行以下操作：

1. 用 `CheckerboardMesh.rectangle(width=1.0, height=0.6, nx=6, ny=4)` 创建规则矩形 Checkerboard geometry；
2. 用 `IsotropicMaterial(young=2.0e5, poisson=0.0, thickness=0.04)` 定义各向同性薄壳材料；
3. 用 `BoundaryConditions.clamped_boundary(mesh, "left")` 固定左边界及其相邻节点线；
4. 用 `NodalLoad.total_force()` 把 `(0.0, 0.0, -1.0e-1)` 的向下合力均分到右边界节点；
5. 创建 `CheckerboardShell(mesh, material, boundary, load, settings)` 并调用 `model.solve()`；
6. 打印 `status`、Newton 迭代数、`residual_inf`、最大位移和前五个变形节点；
7. 调用 `result.plot()` 写出 `outputs/plate_demo.png`（参考/变形网格图），调用 `result.save()` 写出 `outputs/plate_demo.npz`（参考、变形和位移数组）。

对应的核心代码如下（与 `examples/plate_demo.py`（最小矩形板示例脚本）使用同一套 API）：

```python
from checkerboard_shell import (
    BoundaryConditions, CheckerboardMesh, CheckerboardShell,
    IsotropicMaterial, NodalLoad, NewtonSettings,
)

mesh = CheckerboardMesh.rectangle(width=1.0, height=0.6, nx=6, ny=4)
material = IsotropicMaterial(young=2.0e5, poisson=0.0, thickness=0.04)
boundary = BoundaryConditions.clamped_boundary(mesh, "left")
load = NodalLoad.total_force(
    mesh, force=(0.0, 0.0, -1.0e-1),
    node_ids=mesh.boundaries["right"],
    name="downward_tip_resultant",
)
settings = NewtonSettings(
    residual_tolerance=1.0e-9,
    maximum_newton_iterations=20,
)

model = CheckerboardShell(mesh, material, boundary, load, settings)
result = model.solve()
print(result.status, result.max_displacement)
```

## 5. Forward Simulation Workflow

一次标准 forward simulation（正向静态仿真）在代码中的顺序是：

```text
CheckerboardMesh / custom structured vertices
                 ↓
       make_targets(reference, hx, hy)
                 ↓
   IsotropicMaterial + BoundaryConditions + NodalLoad
                 ↓
       CheckerboardShell(mesh, material, ...)
                 ↓
 model.potential(q) → residual(q) / tangent(q)
                 ↓
             NewtonSolver
                 ↓
           SimulationResult
                 ↓
        result.deformed_vertices / result.plot()
```

具体来说，`CheckerboardShell` 首先把自由位移坐标 `q` 还原成完整节点位置；`mechanics/energy.py`（能量模块）的 `paper_components()`（能量分块函数）计算 B/W metric、曲率、`X_M`、`X_W_rel` 和总内能；外载荷通过 `load.values * positions` 形成死载荷势；`residual()` 是总势能对自由坐标的梯度，`tangent()` 是其 Jacobian；最后 `NewtonSolver` 在同一势能上迭代并返回结果。

## 6. Creating New Geometry

### 使用内置矩形网格

推荐从 `CheckerboardMesh.rectangle()` 开始：

```python
from checkerboard_shell import CheckerboardMesh, CheckerboardShell, IsotropicMaterial

mesh = CheckerboardMesh.rectangle(
    width=2.0, height=1.0, nx=12, ny=6,
    origin=(-1.0, -0.5, 0.0),
)
model = CheckerboardShell(mesh, IsotropicMaterial(young=2e5, thickness=0.04))
result = model.solve()
```

如果没有显式传入 boundary condition，`CheckerboardShell` 会使用左边界的默认 `clamped_boundary`；如果没有传入 load，则使用 `NodalLoad.zero()`。

### 使用自定义曲面节点

自定义曲面仍需保持结构化 raw 网格。可以先生成矩形拓扑，再替换节点坐标：

```python
import torch
from checkerboard_shell import CheckerboardMesh, CheckerboardShell, IsotropicMaterial

flat = CheckerboardMesh.rectangle(width=1.0, height=0.6, nx=6, ny=4)
x = flat.reference_vertices[..., 0]
y = flat.reference_vertices[..., 1]
z = 0.05 * torch.sin(torch.pi * x) * torch.sin(torch.pi * y / 0.6)
curved_vertices = torch.stack((x, y, z), dim=-1)

mesh = CheckerboardMesh(
    reference_vertices=curved_vertices,
    hx=flat.hx,
    hy=flat.hy,
    quads=flat.quads,
    boundaries=flat.boundaries,
)
model = CheckerboardShell(mesh, IsotropicMaterial(young=2e5, thickness=0.04))
result = model.solve()
```

这里的 `quads`（quad connectivity）和 `boundaries`（边界节点编号）必须与 `[ny + 1, nx + 1]` 的 row-major 索引一致。当前 B/W 算子直接使用结构化张量切片，因此任意非结构化 connectivity、缺失节点或少于 3 个单元的方向都不受支持。

## 7. Adding New Loads

### 当前支持：节点死载荷

最直接的方法是使用 `NodalLoad.total_force()`：

```python
load = NodalLoad.total_force(
    mesh,
    force=(0.0, 0.0, -0.1),
    node_ids=mesh.boundaries["right"],
    name="right_edge_force",
)
```

如果每个节点的力不同，可以直接构造与 mesh 同形状的 PyTorch 张量：

```python
import torch
from checkerboard_shell import NodalLoad

values = torch.zeros_like(mesh.reference_vertices)
values[..., 2] = -1.0e-3
values[:, -1, 2] = -2.0e-3       # 仅示例：右边界额外向下力
load = NodalLoad(values=values, name="custom_nodal_load")
```

随后把 `load` 传给 `CheckerboardShell(mesh, material, boundary, load)`。`NodalLoad.scaled(factor)`（载荷缩放函数）可构造一个比例载荷。

### 当前不直接支持的载荷

代码中没有 `PressureLoad`、面 quadrature load 或 distributed-load 专用类。若要加入与变形无关的压力或连续分布载荷，需要在模型外先依据面积/积分规则转换成节点力 `[ny + 1, nx + 1, 3]`，再用 `NodalLoad(values=...)` 传入。未来如果要把这类载荷做成正式接口，建议新增 `mechanics/loads.py`（节点载荷模块）中的载荷构造器，而不要把载荷公式写进 `geometry/`（几何目录）或 `solver/`（求解器目录）。随当前法向变化的 follower pressure 不符合现有固定 `load.values` 接口，还需要扩展 `CheckerboardShell.potential()` 的外势组装。

## 8. Boundary Conditions

### 固定节点

`BoundaryConditions.fixed_nodes()` 接受 mesh 的扁平节点编号和布尔分量 `(x, y, z)`：

```python
boundary = BoundaryConditions.fixed_nodes(
    mesh,
    node_ids=mesh.boundaries["left"],
    components=(True, True, True),
)
```

只固定某些分量时，例如只固定 z 位移：

```python
boundary = BoundaryConditions.fixed_nodes(
    mesh,
    node_ids=mesh.boundaries["bottom"],
    components=(False, False, True),
)
```

### 固定边界线与固支边界

```python
fixed = BoundaryConditions.fixed_boundary(mesh, name="left")
clamped = BoundaryConditions.clamped_boundary(mesh, name="left")
```

`fixed_boundary()`（单边界位移约束）只固定命名边界上的节点；`clamped_boundary()`（离散固支约束）会同时固定边界线和相邻节点线，以约束 position-only 离散模型中的边界切向。可用的 `name` 是 `"left"`、`"right"`、`"bottom"`、`"top"`。

也可以直接创建布尔 `fixed_mask`（形状与 `mesh.reference_vertices` 相同），再调用 `BoundaryConditions(fixed_mask=mask)`。自由坐标数量由 `free_mask` 决定，不能在求解中改变。当前约束值均为相对于参考构型的零位移；没有内置非零 prescribed displacement 接口。

## 9. Material Definition

默认材料是 `IsotropicMaterial`：

```python
from checkerboard_shell import IsotropicMaterial

material = IsotropicMaterial(
    young=2.0e5,   # Young modulus
    poisson=0.30,  # Poisson ratio, must be in (-1, 1)
    thickness=0.04,
)
print(material.lame_alpha, material.lame_mu)
```

`CheckerboardShell` 当前的能量接口读取 `material.young`、`material.poisson` 和 `material.thickness`，并通过 `lame_surface_coefficients()` 和 `stvk_density()` 使用各向同性 StVK 收缩。替换材料参数只需传入另一个 `IsotropicMaterial` 实例；如果要增加各向异性或新的本构形式，需要扩展 `material/isotropic.py`（材料参数/本构模块）和 `mechanics/energy.py`（能量组装模块），同时保持 `paper_components()` 的几何输入输出契约。

## 10. Solver Output

`model.solve()` 返回 `SimulationResult`（正向求解结果对象）：

```python
result = model.solve()

result.status              # "CONVERGED", "MAXIMUM_ITERATIONS" 或 "LINE_SEARCH_FAILED"
result.reference_vertices  # [ny + 1, nx + 1, 3] 参考坐标
result.deformed_vertices   # [ny + 1, nx + 1, 3] 变形坐标
result.coordinates         # 仅自由 DOF 的位移坐标
result.displacement        # deformed_vertices - reference_vertices
result.max_displacement    # 最大节点位移（float）
result.energies            # 标量结果：能量分块、total 和 X_M
result.history             # 每次 Newton 迭代的 residual/potential 记录
result.residual_inf        # 最后一次记录的残量无穷范数
result.iterations          # Newton 迭代次数
```

能量和兼容量可以从模型再次读取：

```python
parts = model.energy_components(result.coordinates)
print(parts["membrane_b"])
print(parts["membrane_xm"])
print(parts["bending"])
print(parts["bending_compatibility"])
print(parts["X_M"])
print(parts["X_W_rel"])
```

可视化和保存：

```python
result.plot(save_path="outputs/custom_deformation.png")
result.save("outputs/custom_deformation.npz")
```

`plot()`（三维参考/变形网格绘图方法）返回 Matplotlib figure；`save()`（结果数组保存方法）保存 reference、deformed、displacement、status 和 residual_inf。求解失败时仍会返回当前最好的迭代状态，调用方应检查 `result.status` 和 `result.residual_inf`。

## 11. Future Extension: FPC Inverse Design

当前不实现 inverse design 优化器，但 API 已把 forward 所需输入拆开，未来可以采用：

```text
design parameters
        ↓
mesh / material / boundary / load
        ↓
CheckerboardShell.solve()
        ↓
deformation response
        ↓
optimization / inverse-design objective
```

例如，设计参数可以映射到网格尺寸或曲面高度、`IsotropicMaterial` 的厚度和模量、`NodalLoad` 的节点力，外层优化器再读取 `SimulationResult.deformed_vertices`、`displacement` 或能量目标。未来 FPC actuator design（FPC 执行器设计）可以在这层 forward engine 之上增加专门的电磁/驱动模型，而不需要把设计变量写死在 Checkerboard 几何或 Newton 代码中。

## 12. Developer Notes

### 分离原则

- `geometry/` 只计算结构化节点位置派生的 B/W 几何和参考量，不知道材料、载荷和求解策略；
- `mechanics/` 负责把几何量和材料系数组合成能量，并定义节点死载荷；
- `solver/` 只处理自由度约束、势能梯度/Hessian 和 Newton 迭代，不重新实现壳体数学；
- `visualization/` 只读取结果，不参与求解；
- `CheckerboardShell` 是组装层，适合连接未来的外部 design parameters。

### 推荐扩展位置

- 新的结构化曲面：扩展 `geometry/mesh.py` 或构造 `CheckerboardMesh` 实例；
- 新的几何诊断：放在 `geometry/core.py` 或 `geometry/boundary.py`；
- 新的节点载荷构造器：放在 `mechanics/loads.py`（节点载荷模块）；
- 新材料/本构：扩展 `material/`（材料目录）并在 `mechanics/energy.py`（能量组装模块）中显式接入；
- 新求解策略（例如载荷延续）：放在 `solver/`（求解器目录），复用现有 `potential()`、`residual()` 和 `tangent()`；
- 新绘图：放在 `visualization/`（可视化目录），保持 mechanics/solver 不依赖 Matplotlib。

### 不建议直接修改的部分

`geometry/boundary.py`（有限边界闭合模块）、`geometry/intrinsic.py`（`X_M` 模块）、`mechanics/energy.py`（论文冻结能量）和 `solver/newton.py`（Newton 核心）共同构成当前数学基线。若确实要改动，应同步更新 `tests/test_forward.py`（forward 回归测试），并检查参考状态零能量、residual/tangent 一致性和 demo 收敛。

### API/范围边界速查

| 需求 | 当前接口/状态 |
| --- | --- |
| 规则矩形 shell | `CheckerboardMesh.rectangle()` |
| 结构化自定义曲面 | 直接构造 `CheckerboardMesh(reference_vertices, hx, hy, quads, boundaries)` |
| 节点合力 | `NodalLoad.total_force()` |
| 任意节点力 | `NodalLoad(values=...)` |
| 压力/连续分布载荷 | 当前需外部转换为节点力 |
| follower pressure | 当前不支持，需要扩展外势接口 |
| 固定节点/分量 | `BoundaryConditions.fixed_nodes()` |
| 单边界位移固定 | `BoundaryConditions.fixed_boundary()` |
| 离散固支 | `BoundaryConditions.clamped_boundary()` |
| 非零 prescribed displacement | 当前不支持 |
| 静态正向求解 | `CheckerboardShell.solve()` |
| 载荷延续 | 当前不支持；`solve(load_factor=...)` 只求单一载荷水平 |
| residual/tangent | `CheckerboardShell.residual()` / `.tangent()` |
| 逆向设计 | 尚未实现，预留外层调用 |
| 非结构化网格 | 当前不支持 |
