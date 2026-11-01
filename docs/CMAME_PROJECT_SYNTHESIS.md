# Checkerboard Shell：工作全景、数值证据与 CMAME 投稿前审计

> 审计日期：2026-10-26  
> 审计范围：`/home/song/projects/checkerboard_shell`（当前工作区）与 `/home/song/projects/CB-shell`（arXiv 稿对应仓库，tag `v1.0-arxiv`，commit `093ad89`）  
> 目的：在正式撰写 CMAME 版本之前，把“已经做了什么、证据在哪里、哪些结论能说、哪些还不能说”一次梳理清楚。

---

## 0. 一页速览

### 0.1 这个工作是什么

一句话：**把 Checkerboard/Varignon 边中点几何做成一个只有节点位置（3V 个标量）的 Koiter 型薄壳离散，并用“商空间盲模分析 + 尺度匹配补全”解决中点几何在格子尺度上缺失的膜/曲率可观测性。**

它不只是“用边中点面做壳”。真正的新意有三层：

1. **几何层**：任意空间四边形的四条边中点构成严格共面的 Varignon 平行四边形，因此局部切平面、法向、第一基本形式可以只由节点位置唯一确定，不需要 director、rotation 或对角三角化。
2. **变分层**：把“raw checkerboard gauge”（表示冗余）与“physical blind modes”（真实缺失的物理方向）严格区分；证明膜侧缺 1 个方向、曲率侧缺 2 个方向，并构造 `X_M`（秩一、膜、O(t)）与 `X_W^rel`（客观、参考相对、O(t³)）做定向补全。
3. **实现层**：精确非线性切线可写成 `K_local (稀疏、有界模板) + α g gᵀ (全局秩一)`，用 bordered system 精确求解；边界重构从四个稠密矩阵改成每 W 点一个 3×3 局部 stencil。

### 0.2 论文与代码状态

| 项目 | 状态 |
|---|---|
| arXiv 初稿 | **arXiv:2608.30195**，*Checkerboard Shells: A Position-Only Thin-Shell Discretization with Scale-Compatible Completion*，39 页、8 图，2026-08-31 提交 |
| 作者 | Junjie Song, Xuanyu Wu, Zhifeng Zhang, Bingtao Hu, Zhaoxi Hong, Xiuju Song, Yixiong Feng, Jianrong Tan |
| 公式/能量 | 已“冻结”，`paper_components()` 与 arXiv 稿一致 |
| 当前包测试 | `19 passed`（allatom 环境，Python 3.10 / torch 2.3.1 / numpy 1.26.4 / scipy 1.15.3） |
| plate demo | `CONVERGED`，7 次 Newton，`residual_inf = 8.7e-13` |
| CMAME 加固 | Round 1（审计，已废弃）→ Round 2（驱动恢复、MidedgeTan、稀疏/秩一、scaling、畸变）→ Round 3（β 复核、最终对比表、图集、稿件映射） |
| 代码版本状态 | **不干净**：6 个已跟踪文件被修改，`experiments/` 及多个新模块未提交（审计时 10 个新条目，加上本文件后更多）；`README.md` 仍只描述 dense Newton，未写 sparse/rank-one 与 experiments |

### 0.3 最强的三个卖点

1. **严格平面 Varignon 几何 + 3V position-only 状态**：对任意 warped quad 都有唯一局部 frame，膜和弯共用同一个中点面。
2. **盲模分类与定向补全**：raw gauge 与物理盲模分开；`X_M` 只补 1 个膜缺陷，`X_W^rel` 补 2 个曲率缺陷且保持有限转动客观性。这是论文最“可发表”的理论贡献。
3. **精确稀疏 + 秩一实现**：非线性切线保留 `2c X_M D²X_M` 的局部项，全局项只存一个向量；n=128 可行，augmented nnz 近线性。这是 CMAME 很看重的“方法 + 实现”组合。

### 0.4 最大的五个风险

1. **非线性 benchmark 尚未在当前 production 代码里重跑**：arXiv 的 cantilever/hemisphere/open-cylinder 结果来自原 `CB-shell` 仓库；当前 workspace 只有一个小 cantilever tip-load 路径和“归档 hemisphere 末态切线”检查。
2. **Plate benchmark 是标量横向分支**：只用竖向 DOF（n=16 时 active=225），`X_M` 完全不激活；用它和 MidedgeTan 的 1025 个 active states 比“per-state”并不公平。
3. **M3 的 MidedgeTan 对比不是完全同口径**：Checkerboard 用物理中点矩重构观测，MidedgeTan 用 raw vertex；BC 的离散执行方式也不同；MidedgeTan 的 M3 误差序列 3.74% → 5.01% 非单调，需要解释/重跑。
4. **Roof 参考值有 0.3024 vs 0.30059 的歧义**：Checkerboard 在 n=32 仍为 0.31008（对 modern 参考约 3.16%），MidedgeTan 已到 0.29981；用 0.3024 会“美化”Checkerboard 的 roof 误差。
5. **鲁棒性与性能不能过度声称**：M3 在 δ=.30 表面保形重分配下误差 20.12%（smooth）/28.89%（random）；跨 backend 的 CPU/RSS 比较不公平；没有 MITC4。

---

## 1. 科学工作：从几何到壳能

### 1.1 问题背景

薄壳的膜能 ~ `t`、弯能 ~ `t³`。当 `t → 0` 时，任何虚假的膜约束都会被 `1/t²` 放大，造成 membrane locking。因此，一个“只用位置”的壳离散要成立，必须同时满足：

- 几何上：从当前位置能唯一、客观地得到第一、第二基本形式；
- 变分上：离散膜/弯观察的核空间不能包含非物理的软模，也不能把真实的近等距弯曲模误杀；
- 尺度上：补全项必须进入正确的 `O(t)` / `O(t³)` 扇区，不能把弯曲模变成膜模。

Checkerboard 边中点几何恰好提供了第一条的一个强性质：**Varignon 平面性**。论文的贡献是把后两条也做完。

### 1.2 离散几何

raw 状态是结构化网格上的 `q_ij ∈ R³`。物理中点面由水平/垂直边中点组成：

```
m^x_{i+1/2,j} = (q_ij + q_{i+1,j}) / 2
m^y_{i,j+1/2} = (q_ij + q_{i,j+1}) / 2
```

每个 raw quad 的四个边中点构成 B face（Varignon 平行四边形），相邻 B face 共享中点节点。B frame：

```
t1_B = (q_{i+1,j+1} - q_ij) / ℓ
t2_B = (q_{i,j+1} - q_{i+1,j}) / ℓ
a_B  = T_Bᵀ T_B
```

W 点位于 raw vertex，第二基本形式由相邻 B 法向的中心差与 W 点平均切向配对：

```
b̃_W,αβ = -D_α^h N_B · t_β^W
b_W     = sym(b̃_W)
```

边界用 moment-exact 单边重构：cell-center 样本在 `1/2, 3/2, 5/2`，得到值权重 `(15/8, -5/4, 3/8)`、导数权重 `(-2, 3, -1)/h`，二维边/角用张量积；dual fraction 为 `1, 1/2, 1/4`。

### 1.3 表示冗余 vs 物理盲模

这是论文最核心的区分：

- **raw gauge**：`g_ij = (-1)^{i+j} A`，`A ∈ R³`。它不改变任何中点坐标，只是 raw 表示冗余。对连通二分网格，midpoint map 的核恰好是这 3 维（每个分量一个 alternating 模式）。
- **physical blind modes**：在商掉 gauge 和刚体运动后，B metric 还缺 1 个膜方向；对称 W 曲率在 flat transverse sector 还缺 2 个曲率方向（`ip`, `jp`，`p = (-1)^{i+j}`）。

论文给出 flat 情形的完整分类：

```
ker J_B = R_h + G_h + span{z_M}
ker b_W = span{1, i, j, p, ip, jp}
```

其中 `p` 是 gauge，`ip, jp` 是真正的曲率盲模。这个分类解释了“为什么几何二阶收敛还不够”：光滑一致性不能排除格子尺度的盲模。

### 1.4 补全与能量

**膜补全 `X_M`**：由 B frame 的 polar normalization、当前/参考相对转动、相邻 B frame 的 intrinsic connection、plaquette mixed connection moment 与 M-character 投影构成。它只补那 1 个膜方向，归一化为 `β_M = 8μ`（trace-free shear 通道；不是 benchmark 拟合，但也不是经典 StVK 对 skew fiber 的唯一推论）。

**曲率补全 `X_W^rel`**：

```
X_W       = N_NE - N_NW + N_SW - N_SE
X_W^rel   = X_W(q) - R_W(q,q₀) X_W(q₀)
```

`R_W` 是 W frame 的 polar factor（surface-polar transport），保证有限转动下的客观性；系数 `β_W⁰ = μ γ_W⁰`，由参考 metric 上的 auxiliary-Q1 法向梯度能量精确积分得到。

**最终能量**：

```
E_m   = t/4 [ Σ_B A_B⁰ Q_{a_B⁰}(a_B - a_B⁰) + 8μ X_M² ]
E_b   = t³/12 Σ_W θ_W A_W⁰ Q_{a_W⁰}(b_W - b_W⁰)
E_XW  = t³ μ/12 Σ_{W∈int} γ_W⁰ ||X_W^rel||²
Π_h   = E_m + E_b + E_XW - W_ext
```

**切线的精确结构**（`c = 2tμ`）：

```
E_XM   = c X_M²
R_XM   = 2c X_M g,              g = D X_M
K_XM   = 2c [ g gᵀ + X_M D²X_M ]
```

`D²X_M` 是局部的（每个耦合在逻辑方向两步以内）；只有 `g gᵀ` 是全局的。production 实现把它存成 `K_local + α g gᵀ`，用 bordered system 精确求解，不显式组装外积。

### 1.5 解析结果

论文中可当作“严格定理/解析命题”的：

- 中点映射的 exact checkerboard gauge 核与商空间解释；
- Varignon B face 的严格平面性、B/W 几何的 gauge 不变性；
- `a_B`, `b_W` 在光滑族上的二阶一致性；边界 moment 重构的二阶局部一致性；
- flat first-form 的一维 quotient defect 分类，`X_M` 的 quotient-minimal 补全；
- flat transverse `b_W` 的核分类，`X_W` 的曲率补全；
- surface-polar `X_W^rel` 的客观性；
- aligned generalized cylinder 上的 `O(h²)` near-isometry approximation（**只对 arc-length、grid-aligned、row-independent isometry 这一族**）；
- `K(t) = t K_m + t³ K_b` 下的模态衰减关系 `t²/(λ_j + t²)`，解释 fixed-h ultra-thin residual locking。

不能当作定理的：一般曲面的 uniform locking-free 收敛、任意拓扑/多 patch、variable-step、fixed-h 极薄极限消失。

### 1.6 与现有方法的关系

论文的定位是“几何驱动的 position-only 壳离散”，与以下路线区分：

- **MITC / assumed strain / mixed**：通过应变重构或混合变量缓解 locking；Checkerboard 不引入这些变量，而是从中点几何的盲模分析出发做定向补全。
- **IGA / subdivision / DDG shells**：需要高连续性、细分或独立 edge/director 状态；Checkerboard 只用 3V。
- **Midedge shells**（Grinspun/Chen 等）：也从中点/边量构造弯曲，但通常保留 edge 状态；Checkerboard 的 B/W staggered 结构把膜和弯放在同一个中点载体上。
- **Dellinger 2024 的离散基本形式**：同属 Checkerboard 几何谱系；论文把 first form 对齐后做了比较，second form 的构造不同（B-face normal vs Gauss image）。

CMAME 版本的 Introduction 需要把这张“方法谱系表”做得更清楚，否则审稿人容易把它归类为“graphics/DDG 的壳应用”，低估变分分析的贡献。

---

## 2. 代码资产地图

### 2.1 目录结构

```
checkerboard_shell/
├── src/checkerboard_shell/
│   ├── geometry/
│   │   ├── mesh.py          结构化矩形网格
│   │   ├── core.py          B/W 几何、第一/第二基本形式
│   │   ├── boundary.py      有限边界 moment 闭合、局部 stencil、X_W^rel transport
│   │   └── intrinsic.py     polar B triad、intrinsic connection、X_M
│   ├── mechanics/
│   │   ├── energy.py        冻结的 paper 能量分块
│   │   ├── loads.py         节点死载荷
│   │   ├── low_rank.py      精确 sparse-local + rank-one 切线
│   │   └── model.py         CheckerboardShell 组装、dense/sparse 两条 solve 路径
│   ├── solver/
│   │   ├── newton.py        dense Newton 参考路径
│   │   ├── nonlinear.py     通用 NonlinearProblem + Newton-GMRES/CG
│   │   ├── continuation.py  自适应标量 continuation
│   │   ├── config.py        通用 SolverConfig
│   │   ├── convergence.py   残量/收敛监控
│   │   ├── boundary_conditions.py  固定 DOF 约束
│   │   └── result.py        SimulationResult
│   ├── material/isotropic.py
│   └── visualization/plotting.py
├── examples/plate_demo.py
├── tests/                   19 tests
├── docs/                    核心管线 + solver 架构/迁移文档
├── experiments/
│   ├── cmame_additional/    Round 1（审计，已废弃）
│   ├── cmame_round2/        Round 2（驱动恢复、稀疏/秩一、scaling、畸变）
│   └── cmame_round3/        Round 3（β、最终对比、图集、稿件映射）
├── outputs/plate_demo.*
└── README.md
```

### 2.2 核心模块与数学对应

| 数学对象 | 代码位置 | 说明 |
|---|---|---|
| B 对角切向、B 法向、`a_B` | `geometry/core.py:47-84` | 冻结的 B 几何 |
| W 平均切向、`b̃_W`, `b_W`, `χ` | `geometry/core.py:64-84` | 内部 W 几何 |
| 边界 moment 闭合 | `geometry/boundary.py:95-219` | `_axis_stencil`, `_tensor_reconstruct`, `boundary_geometry` |
| 局部 3×3 stencil | `geometry/boundary.py:46-87` | `LocalBoundaryStencil`，O(V) 存储 |
| `fast_boundary_geometry` | `geometry/boundary.py:319-362` | 与慢路径代数等价 |
| `γ_W`, `A_B`, `A_W`, dual fraction | `geometry/boundary.py:235-273` | 参考量缓存 |
| `X_W^rel` transport | `geometry/boundary.py:276-289` | polar Newton + surface-polar transport |
| `X_M` | `geometry/intrinsic.py:56-89` | polar triad、connection、M-character |
| 能量分块 | `mechanics/energy.py:28-76` | `membrane_b`, `membrane_xm`, `bending`, `bending_compatibility` |
| 精确 sparse+rank-one | `mechanics/low_rank.py:214-311` | `assemble_paper_low_rank_tangent` |
| 彩色 AD Hessian | `mechanics/low_rank.py:142-211` | `assemble_coloured_local_hessian` |
| bordered solve | `mechanics/low_rank.py:59-114` | `LowRankTangent.augmented_matrix/solve` |
| 模型组装 | `mechanics/model.py:21-186` | `CheckerboardShell` |
| dense Newton | `solver/newton.py:31-124` | 参考路径 |
| 通用 Newton-GMRES/CG | `solver/nonlinear.py:24-351` | 与壳模型解耦 |
| 自适应 continuation | `solver/continuation.py:26-83` | 与壳模型解耦 |

### 2.3 三条求解路径（要分清）

1. **Dense reference path**：`CheckerboardShell.solve(backend="dense_reference")` → `solver/newton.py`。用 `torch.autograd.functional.jacobian` 形成完整 Hessian，特征分解 + 线搜索。适合小网格/公式回归，不适合 scaling。
2. **Production sparse low-rank path**：`CheckerboardShell.solve(backend="sparse_low_rank")`（默认）→ `mechanics/model.py:_solve_sparse_low_rank` → `low_rank.py`。残量由 AD，切线由彩色 HVP 组装局部稀疏矩阵 + 一个秩一向量；`LowRankTangent.solve` 解 bordered system。适合中大规模，但 Newton 循环是固定载荷、无 continuation。
3. **Generic Newton-Krylov path**：`solver/nonlinear.py` + `solver/continuation.py`。这是从 FPC 项目提取的通用基础设施，**目前没有接入 `CheckerboardShell`**。它提供 `NonlinearProblem`、GMRES/CG、线搜索、`adaptive_continuation`，但 Shell 的 `_solve_sparse_low_rank` 没有用它。若要让 production 代码跑完整非线性载荷路径，最自然的做法是写一个 `NonlinearProblem` 适配器，让 `adaptive_continuation` 驱动现有的 sparse tangent。

### 2.4 测试与运行状态

- `tests/test_forward.py`：参考状态零能量、零载荷收敛、非零载荷变形、residual/tangent 形状/有限性/对称性。
- `tests/test_sparse_production.py`：局部 stencil 与慢路径一致（`rtol=2e-13`）、精确非线性切线 vs dense Hessian（`rtol=3e-10`）、`xm_multiplier` 只缩放秩一项。
- `tests/test_nonlinear_solver.py`：通用二次问题、GMRES/CG、预条件、continuation、配置校验。
- 运行：`PYTHONPATH=src /home/song/miniconda3/envs/allatom/bin/python -m pytest -q` → **19 passed**。
- plate demo：`CONVERGED`，7 iterations，`residual_inf=8.7e-13`，最大位移 `3.60e-2`。

### 2.5 代码层面的问题（投稿前要处理）

1. **未提交/未打 tag**：`git status` 显示 6 个修改 + 10 个新增；`experiments/` 全部 untracked。CMAME 的可复现性要求“论文图表 ↔ commit/tag ↔ 脚本”闭合。
2. **README 过期**：只写 dense Newton，没有 `low_rank.py`、`sparse_low_rank`、`adaptive_continuation`；读者会以为 production 还是 dense。
3. **无通用 gauge quotient API**：`BoundaryConditions` 只做固定 DOF；benchmark 里的 gauge 行是手工加的。论文说“raw gauge 用 quotient reduction 去掉”，但包本身没有暴露这个能力。至少要在文档/示例里给出 `rigid+gauge basis + orthogonal complement` 的用法。
4. **非线性驱动缺失**：当前 workspace 没有 cantilever/hemisphere/open-cylinder 的 production 驱动；`experiments/cmame_round2/run_exactness.py` 只验证了一个小 cantilever tip-load 路径和归档 hemisphere 末态切线。
5. **`solver/nonlinear.py` 与 Shell 解耦**：通用求解器没有被 `CheckerboardShell` 使用，存在“两套 Newton”的维护风险。CMAME 版本最好统一到一条 production 路径。

---

## 3. 数值验证资产

### 3.1 论文原有验证（arXiv 稿）

论文的数值部分按“几何 → 机理 → 线性 → 非线性”组织：

#### (a) 几何收敛（4 个光滑曲面）

- cylinder、sphere、twist、skew；
- `b_W` 与 Dellinger 的 Checkerboard second form 对比；
- 四个曲面均约二阶收敛；twist/skew 两条曲线几乎重合，cylinder/sphere 误差常数不同但同阶。

**证据位置**：`/home/song/projects/CB-shell/paper_draft_cmame/figures/figure_04_geometry_convergence.*`，数据在 `CB-shell` 的 paper benchmark outputs。

#### (b) 补全机理

- flat patch 的 exact nullity：B-only 6，B+X_M 5，full B/W 5；理论 rigid+gauge 维度 5 → `X_M` 恰好去掉 1 个方向。
- M3 active positive spectrum：`B+X_M` 基本保留 B 侧最小正奇异值序列；full B/W 把软空间抬高几个量级（over-complete 对照）。
- M3 载荷响应：B 与 B+X_M 的最大相对位移差 `1.19e-10`；`X_M` 最大膜能占比 `5.41e-12`。
- `X_W^rel` 对 plate/roof/M3 的有限网格位移影响：`0.0136% / 0.0210% / 0.2983%`，随加密下降。
- 附录 curved free patch nullity：B-only 51，B+X_M 50，full B/W 9；M3 active n=16：10 / 9 / 0。

**证据位置**：`figure_05_minimal_completion.*`；`CB-shell/outputs/checkerboard_benchmark`、`outputs/checkerboard_mechanism*` 等。

#### (c) 线性标准 benchmark

| benchmark | 网格 | Checkerboard 误差 | MidedgeTan 误差 | 备注 |
|---|---:|---:|---:|---|
| simply supported plate | n=8/16/32/64 | 5.596 / 1.488 / 0.399 / 0.104% | 4.953 / 1.248 / 0.313 / 0.078% | 参考 4.4360891e-4 |
| Scordelis–Lo roof | n=12/16/24/32 | 11.45 / 6.89 / 3.66 / 2.54% | 9.10 / 5.50 / 2.20 / 0.856% | 参考 0.3024；modern 0.30059 |
| MacNeal–Harder M3 | n=6/8/12/16 | 42.23 / 15.85 / 2.20 / 0.1105% | 43.12 / 13.21 / 3.741 / 5.009% | 参考 0.09352155；MT 非单调 |

- plate 的 Checkerboard 驱动是**标量横向分支**：只解竖向 DOF（n=16 active 225），不是完整 3V 求解；因为平板横向与面内解耦，横向响应本身是可信的，但 active-state 对比不公平。
- roof 的 Checkerboard 收敛到 0.310 左右，对 modern 0.30059 在 n=32 仍约 3.16%；MidedgeTan 收敛到 0.2998。用 0.3024 作参考会把 Checkerboard 的 roof 误差“美化”约 0.6 个百分点。
- M3 的 Checkerboard 在 n=16 达到 0.093418（0.11%），明显好于 MidedgeTan 的 0.098206（5.0%）；但 MT 的 3.741% → 5.009% 非单调，需要核查网格/BC/observable 是否完全同口径。

#### (d) 非线性 benchmark（arXiv 稿的核心非线性证据）

| case | 网格 | 关键结果 | 备注 |
|---|---|---|---|
| end-moment cantilever | 16×4, 24×4, 32×4 | 最终转角 **359.841°**；fine RMS 2.442%，max 2.792%；膜能占比 0.432% | 解析常曲率圆弧参考；position-only 端弯矩势 |
| alternating-force hemisphere | n=6,8,12,16 | n=16：A 终点误差 **2.759%**，B **0.149%**，joint RMS **1.736%**，max 1.859% | M3 几何的有限变形全路径；`X_W^rel` 最大 0.322% |
| open-ended cylinder pull-out | 8×8, 12×12, 16×16 | n=16：A **8.094%**，B **5.316%**，C **1.861%**，combined RMS **4.483%** | 45/45/46 load steps；C 点位移反转 |
| pinched cylinder | — | **未纳入正文**：`NONLINEAR_SOLVER_LIMITED`，path loop detected | 只在讨论里提一句 |

**证据位置**：`CB-shell/outputs/checkerboard_nonlinear_benchmarks/`（summary、raw paths、figures、manifest），驱动在 `CB-shell/experiments/checkerboard_nonlinear_benchmarks/`。

### 3.2 CMAME Round 1（`experiments/cmame_additional/`）

**性质：审计 + 初步补充实验，已被 Round 2/3 取代，不能作为正文证据。**

做了 6 组：

1. **mesh distortion**：75 个 Checkerboard 条件（72 收敛，3 个随机 δ=.40 网格反转）；plate 在 spacing ratio 8 时误差 49.25%；roof/M3 结果混合。
2. **gauge invariance**：对不同 raw gauge、线性变形 M3、非线性 cantilever 末态，物理中点坐标变化 `7.73e-15`，切线相对变化 `1.23e-15`，能量相对变化 `7.34e-13`，残量相对 `1.03e-11`。
3. **scaling**：n=16/32 的 dense tangent、rank-one audit、reference-linear 稀疏原型；n=64/128 被内存 guard 停止，只有投影值。
4. **β sensitivity**：c_beta=0…8；M3/roof 的载荷响应变化很小；但用的是重构驱动。
5. **curved spectrum**：sphere/saddle；full-tangent 最小特征值正，未发现新的近零模；但 BC 是强 clamp。
6. **boundary modes**：flat/cylinder strip；未观察到持续的 free-edge checkerboard 零模；但 2×2 corner quotient 比最小物理 BC 强。

**Round 1 的结论**：方法可跑、gauge 不变性强、但“三方法对比”不成立（MidedgeTan/MITC4 不可用、驱动是重构的、dense Hessian/boundary 不可扩展、spacing-ratio 超出 constant-h 形式范围）。报告中明确写了 **submission premature**。

### 3.3 CMAME Round 2（`experiments/cmame_round2/`）

**性质：解决 Round 1 的两个实现 blocker，并恢复真正的 manuscript drivers。这是当前最可信的一轮。**

#### (a) Driver recovery

- 从 `/home/song/projects/CB-shell` tag `v1.0-arxiv`（commit `093ad89`）恢复 benchmark 驱动、输出、LibShell wrapper。
- 复现论文数值：
  - plate n=8/16/32/64：4.6843e-4 / 4.5021e-4 / 4.4538e-4 / 4.4407e-4；误差 5.596 / 1.488 / 0.399 / 0.104%。
  - roof n=12/16/24/32：0.33704 / 0.32322 / 0.31348 / 0.31008；误差 11.45 / 6.89 / 3.66 / 2.54%。
  - M3 n=6/8/12/16：0.05403 / 0.07870 / 0.09146 / 0.09342；误差 42.23 / 15.85 / 2.20 / 0.1105%。
- 与历史 dense 结果一致：plate 精确一致；M3 n=8/16 相对差 1.31e-9 / 7.64e-10。
- 关键修正：Round 1 的 plate/roof/M3 观测、BC、载荷、点载荷符号、三阶点算子都有错；Round 1 的 displacement 表**不能引用**。

#### (b) 稀疏边界与精确秩一

- 四个稠密 `[V × F_B]` 重构矩阵 → 每 W 点一个 3×3 局部 stencil。
- n=64：528.1 MiB → **1.45 MiB**；n=128：8.13 GiB（old 为按精确形状投影）→ **5.71 MiB**（new measured）。
- 旧/新能量差 0，残量相对差 3.27e-17，切线相对差 6.80e-18。
- `X_M` 的精确切线 `K = K_local + 2c g gᵀ`；`K_local` 保留 `2c X_M D²X_M`。
- 验证：
  - nonlinear cantilever full tangent 1.47e-16；tangent action 7.72e-16；Newton increment 2.54e-12；
  - archived hemisphere final-state tangent 1.37e-16；tangent action 6.36e-16；
  - `X_M D²X_M` 的 Frobenius 范数：cantilever 0.01179，hemisphere 73.38（非零，说明保留了非线性项）。
- `g gᵀ` 满矩阵但秩 1；`H_X` 在 M3 n=6 上 6,669/21,609 非零，耦合不超过两个 raw node step。

#### (c) Scaling（n=16…128，measured）

M3 n=128：

| version | raw / active | local nnz | augmented nnz | assembly median | factor median | peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| CB-0 | 49,923 / 49,151 | 3,721,041 | 3,861,337 | 7.12 s | 18.60 s | 2,078 MiB |
| CB-M | 49,923 / 49,151 | 3,721,041 | 3,961,184 | 7.55 s | 17.38 s | 2,081 MiB |
| CB-Full | 49,923 / 49,151 | 3,721,041 | 3,961,184 | 11.49 s | 16.39 s | 2,073 MiB |
| MidedgeTan | 99,331 / 98,300 | not exposed | not exposed | different backend | different backend | 1,183 MiB |

- CB-M 的 border 只加 `2N+1 = 99,847` nnz，比 CB-0 多 2.59%。
- `g` 的计算 0.094 s，占 CB-M assembly 的 1.25%。
- augmented nnz 对 scalar states 的 log-log 斜率：plate ~1.09，M3 ~1.01。
- **但不能声称 CPU/RSS 优于 MidedgeTan**：Python/AD worker 2.07 GiB / ~28 s vs C++ baseline 1.18 GiB / ~14.5 s；瓶颈是 sparse factorization + AD local Hessian，不是秩一向量。
- 3V 状态数优势是真实的：CB-Full 的 raw states 是 MT 的 50.3%，M3 active states 是 50.0%。

#### (d) 表面保形畸变（revalidated）

- plate n=16、roof n=16、M3 n=12；smooth/random tangential redistribution，δ=0…0.30；边界固定；所有节点仍在同一平面/圆柱/球面；所有单元有效。
- 同一套顶点/载荷/BC 喂给 MidedgeTan（provenance-recorded adapter），uniform gate 复现历史值（最大相对差 2.0e-10）。

| case | method | uniform | smooth δ=.30 | random δ=.30 |
|---|---|---:|---:|---:|
| plate n=16 | CB-Full | 1.488% | 1.601% | 1.543% |
|  | MidedgeTan | 1.248% | 1.354% | 1.267% |
| roof n=16 | CB-Full | 6.886% | 7.502% | 6.635% |
|  | MidedgeTan | 5.504% | 5.852% | 5.956% |
| M3 n=12 | CB-Full | 2.201% | **20.124%** | **28.895%** |
|  | MidedgeTan | 3.741% | **30.219%** | **43.742%** |

- 结论：plate/roof 稳定、响应变化温和；M3 对球面样本重分配敏感，两个方法都退化，Checkerboard 退化较小但仍很大。
- **不能声称** Varignon 平面性带来普遍鲁棒性优势。

#### (e) 回归

- 参考能量 8.93e-29；参考残量 inf 1.37e-13；
- 刚体平移/转动能量相对变化 1.49e-13 / 1.24e-13；
- 切线对称性到存储精度；残量方向导数/切线作用误差 3.91e-9（FD step 1e-7）；
- dense/sparse 能量、残量、切线作用一致；cantilever 两条路径均 10 次 Newton 到同一末态；
- gauge α/L=0…1：物理中点最大变化 8.27e-15，切线相对变化 1.25e-15，变形状态能量相对变化 1.65e-12，参考残量 scaled 3.57e-17。

### 3.4 CMAME Round 3（`experiments/cmame_round3/`）

**性质：稿件收尾，β 复核 + 最终对比表 + 图集 + 章节映射。没有新 benchmark、没有改公式。**

#### (a) β_M 复核

- `β_M = c_beta (8μ)`，c_beta = 0, 0.25, 0.5, 1, 2, 4, 8。
- 相对 c_beta=1 的最大响应变化（0.25–4）：
  - plate：0（标量横向分支，X_M 不激活）；
  - roof n=16：1.63e-5；
  - M3 n=12：4.21e-10。
- 含 c_beta=0 和 8 的最大变化：roof 1.36e-3，M3 5.04e-10。
- completion-added `DX_M` stiffness 严格线性于 c_beta；roof 的 `E_XM/E_int` 在 c_beta=1 为 3.15e-7，M3 为 4.20e-24。
- c_beta=0 在 curved roof/M3 的 benchmark BC 下**没有**新增 full-tangent 零模；只有在 flat plate 的面内 membrane operator 上，nullity 从 5 变 6（多出的方向 base stiffness 9.62e-20，completion stiffness 8.27e-5）。
- 结论：c_beta=1 不是 sharp benchmark-tuned optimum；但不能说 β 不重要。

#### (b) 最终对比表

`02_final_comparison/MANUSCRIPT_COMPARISON_TABLE.md`，12 行：

| Benchmark | n | Checkerboard response (error) | MidedgeTan response (error) | CB raw/active | MT raw/active |
|---|---:|---:|---:|---:|---:|
| plate | 8 | 4.68434e-4 (5.596%) | 4.65579e-4 (4.953%) | 243/49 | 451/257 |
| plate | 16 | 4.50210e-4 (1.488%) | 4.49145e-4 (1.248%) | 867/225 | 1667/1025 |
| plate | 32 | 4.45377e-4 (0.3986%) | 4.44996e-4 (0.3126%) | 3267/961 | 6403/4097 |
| plate | 64 | 4.44070e-4 (0.1039%) | 4.43956e-4 (0.0782%) | 12675/3969 | 25091/16385 |
| roof | 12 | 0.337039 (11.45%) | 0.274884 (9.099%) | 507/468 | 963/912 |
| roof | 16 | 0.323223 (6.886%) | 0.285756 (5.504%) | 867/816 | 1667/1600 |
| roof | 24 | 0.313479 (3.664%) | 0.295748 (2.2%) | 1875/1800 | 3651/3552 |
| roof | 32 | 0.310083 (2.541%) | 0.299811 (0.8562%) | 3267/3168 | 6403/6272 |
| M3 | 6 | 0.0540279 (42.23%) | 0.0531934 (43.12%) | 147/107 | 267/212 |
| M3 | 8 | 0.0787009 (15.85%) | 0.0811630 (13.21%) | 243/191 | 451/380 |
| M3 | 12 | 0.0914627 (2.201%) | 0.0970198 (3.741%) | 507/431 | 963/860 |
| M3 | 16 | 0.0934182 (0.1105%) | 0.0982063 (5.009%) | 867/767 | 1667/1532 |

**这张表的数据来源需要标注清楚**：12 行的 MidedgeTan 数值全部来自 `driver_recovery.json` 中的 recovered historical rerun；新的 MidedgeTan distortion adapter 只在 plate n=16、roof n=16、M3 n=12 三个点做了 uniform gate。M3 的 MT raw states 是按 `6n² + 8n + 3` 反推的（历史 wrapper 只暴露 reduced dimension），不是逐项 measured。CMAME 版本应把 12 行全部用同一个 adapter 重跑，并尽量直接读出 raw/active states。

另：表中 roof 误差按 0.3024 计算；对 modern 参考 0.3005924566，Checkerboard n=32 为 3.16%，MidedgeTan n=32 为 −0.26%。正文应同时给出两个参考或明确选择理由。

#### (c) 图集与稿件映射

- `03_final_figures/`：Figure A1（gauge 物理商）、A2（gauge 误差）、B（accuracy per raw state）、C（sparse/rank-one scaling）、D（surface-preserving distortion）、E（β sensitivity）。
- `04_manuscript_map/CMAME_ROUND3_MANUSCRIPT_MAP.md`：逐节给出 retain/revise/delete/add，以及 claim-by-claim disposition。
- 这是写 CMAME 版本时最直接的“改稿 checklist”。

### 3.5 证据等级速查

| 主张 | 等级 | 证据 |
|---|---|---|
| Varignon B face 严格平面 | 定理 | 几何定义 |
| raw gauge 核 = `{(-1)^{i+j}A}` | 定理 | 中点映射线性代数 |
| flat B metric 缺 1 个膜方向 | 定理 | 附录 kernel proof |
| flat symmetric `b_W` 缺 2 个曲率方向 | 定理 | 附录 kernel proof |
| `X_M` quotient-minimal | 定理（flat, linear observation rank） | 附录 |
| `X_W^rel` 客观、reference-consistent | 定理 | 附录 objectivity |
| aligned generalized cylinder `O(h²)` | 定理（窄条件） | 附录 |
| 二阶几何收敛 | 数值 | 4 曲面 |
| X_M 去掉 1 个 flat nullity | 数值 | flat nullity 6→5 |
| X_W^rel 位移影响 <0.3% | 数值 | plate/roof/M3 |
| 线性 benchmark 竞争性 | 数值 | 12 行对比表 |
| 非线性路径能力 | 数值（旧 production 驱动） | cantilever/hemisphere/cylinder |
| 稀疏+秩一等价 | 数值 + 代数 | exactness 1e-16 |
| n=128 可行 | 数值 | scaling |
| 畸变鲁棒性 | **未证明/混合** | M3 20–29% |
| CPU/RSS 优势 | **未证明** | 跨 backend |
| 一般曲面 locking-free | **未证明** | 只有窄定理 |
| MITC4 对比 | **缺失** | 无 provenance-closed 实现 |

---

## 4. 关键结果总表（写论文时可直接取用）

### 4.1 状态数

- Checkerboard：`3V` 个 raw scalar states（V 为 raw 节点数）。
- MidedgeTan：`3V + E`（E 为 edge states）。
- 在 12 个 plate/roof/M3 网格对上，Checkerboard 使用 MidedgeTan raw states 的 **50.5–55.1%**。
- 但注意：plate 的 Checkerboard active states 只含竖向 DOF（225 @ n=16），不能与 MT 的 1025 直接比；M3/roof 的 active 对比更有意义（507/431 vs 963/860 @ n=12；867/816 vs 1667/1600 @ n=16）。

### 4.2 实现

- 边界存储：n=64 528.1 MiB → 1.45 MiB；n=128 8.13 GiB → 5.71 MiB。
- M3 n=128：49,923 raw / 49,151 active；local nnz 3,721,041；augmented 3,961,184；assembly 11.49 s；factor 16.39 s；RSS 2,073 MiB。
- CB-M border：+99,847 nnz（2.59%）。
- augmented nnz 斜率：plate 1.09，M3 1.01。
- 秩一存储：1 个向量 + 1 个标量；`rank_one_storage = N + 1`。

### 4.3 Gauge / 客观性

- 物理中点最大变化 8.27e-15；
- 切线相对变化 1.25e-15；
- 变形状态能量相对变化 1.65e-12；
- 参考残量 scaled 3.57e-17；
- 刚体平移/转动能量相对变化 1.49e-13 / 1.24e-13。

### 4.4 精确切线

- cantilever full tangent 1.47e-16；tangent action 7.72e-16；Newton increment 2.54e-12；
- hemisphere archived final-state tangent 1.37e-16；tangent action 6.36e-16；
- `X_M D²X_M` Frobenius：0.01179（cantilever）、73.38（hemisphere）。

### 4.5 畸变

- plate/roof：温和（≤ ~1% 绝对误差变化）；
- M3 δ=.30：CB 20.12% / 28.89%；MT 30.22% / 43.74%；
- 所有 sampled cells 有效；最小 scaled Jacobian 0.889（M3 random .30）。

### 4.6 β

- c_beta=1 在 0.25–4 范围内的响应变化：plate 0，roof 1.63e-5，M3 4.21e-10；
- completion stiffness 线性于 c_beta；
- flat membrane nullity 6→5。

### 4.7 谱

- sphere/saddle full-tangent 最小特征值正：sphere 7.64e-4 / 2.05e-4 / 5.01e-5（n=8/16/32）；saddle 1.52e-3 / 3.64e-4 / 8.74e-5；
- flat/cylinder strip 最小特征值正：flat 1.44e-7 / 2.67e-8 / 5.12e-9；cylinder 1.49e-7 / 2.41e-8 / 4.32e-9；
- 没有发现新的 persistent full-tangent 零模；但 BC 强、无 baseline spectra，只能作为附录诊断。

---

## 5. 证据边界：能说 / 不能说

### 5.1 可以安全写进 Abstract 的

- position-only `3V` state；
- exact planar Varignon B faces；
- raw gauge 与 physical blind modes 的区分；
- 1 个膜 + 2 个曲率盲模的分类；
- `X_M` / `X_W^rel` 的尺度匹配补全；
- gauge invariance 到 roundoff；
- 在 recovered plate/roof/M3 序列上，以约一半的 raw states 达到同阶精度；
- exact nonlinear tangent = sparse local + one rank-one border；
- n=128 feasible，near-linear nnz。

### 5.2 必须限定条件的

- “minimal completion”：只能写 **quotient-minimal in the linear observation-rank sense for the one-dimensional flat membrane defect**。
- “competitive accuracy”：写 **on the recovered plate/roof/M3 sequences**，并说明 roof 对 modern reference 仍偏硬、M3 的 MT 对比非单调。
- “exact planarity”：写 **guarantees a unique local frame, not distortion-insensitive accuracy**。
- “state economy”：写 **raw state count**，不要写成 CPU/memory economy。
- “locking”：写 **aligned generalized-cylinder O(h²) near-isometry approximation + observed fixed-h residual locking**，不要写 locking-free。
- “distortion”：写 **stable on sampled valid meshes; accuracy is problem-dependent; M3 is sensitive**。

### 5.3 不能说的

- “more robust on distorted/nonuniform meshes”（M3 20–29% 反例；spacing-ratio 超出 constant-h 形式范围）。
- “faster / lower memory than standard shells”（跨 backend，Python/AD 更慢、RSS 更高）。
- “MITC4 comparison”（没有 provenance-closed 实现）。
- “arbitrary unstructured / variable spacing”（当前公式要求 structured + constant hx,hy）。
- “general curved-shell convergence / uniformly locking-free”（只有窄定理）。
- “exact zero relative gauge error at stress-free state”（0/0 未定义；应给 absolute/scaled）。
- “X_M removes a curved-shell null mode”（只在 flat membrane operator 上严格；curved benchmark BC 下 nullity 不变）。

---

## 6. CMAME 版本建议

### 6.1 定位

CMAME 读者关心的是：**一个新的壳离散是否可靠、是否比现有方法有可解释的优势、实现是否可扩展、验证是否可复现**。因此 CMAME 版本不应写成“DDG 几何 + 一堆 benchmark”，而应写成：

> 一个由中点几何驱动的 position-only 壳离散；它的变分盲模被严格分类并定向补全；它的精确切线具有 sparse-local + rank-one 结构；它在标准线性/非线性 benchmark 上以约一半的 raw states 达到同阶精度，并明确披露畸变与薄极限的边界。

与 arXiv 稿相比，CMAME 版本要**新增**：
- 实现/算法 section（局部 stencil、彩色 AD Hessian、bordered solve、复杂度）；
- 更公平的 baseline 对比（同 mesh/load/BC/observable）；
- 更明确的 limitations（M3 畸变、roof 参考、plate 标量分支）；
- 可复现性声明（tag、脚本、环境、manifest）。

### 6.2 标题建议

1. 保留原题，加副标题：*Checkerboard Shells: A Position-Only Thin-Shell Discretization with Blind-Mode Completion and an Exact Sparse Rank-One Tangent*
2. 更 CMAME 风格：*A position-only thin-shell discretization from checkerboard midpoint geometry: quotient analysis, compatible completion, and scalable implementation*
3. 更短：*Blind-mode-completed checkerboard shells: geometry, analysis, and scalable computation*

### 6.3 建议结构（双栏，约 28–32 页）

| Section | 内容 | 页数 |
|---|---|---:|
| 1 | Introduction：locking、position-only 动机、方法谱系表、贡献 | 4 |
| 2 | Checkerboard geometry and midpoint quotient：raw/midpoint、B/W、boundary closure | 5 |
| 3 | Blind modes：gauge vs physical、flat kernel classification、curved observations | 4 |
| 4 | Compatible completion and energy：X_M、X_W^rel、coefficients、energy、scale separation | 5 |
| 5 | Analytical properties：objectivity、reference consistency、cylinder theorem、modal attenuation | 3 |
| 6 | Numerical implementation：local stencil、colored AD Hessian、rank-one border、complexity | 3 |
| 7 | Verification and benchmarks：geometry convergence、mechanism、linear、nonlinear、scaling、distortion | 8 |
| 8 | Discussion and limitations | 2 |
| 9 | Conclusions | 1 |
| Appendix/Supplement | proofs、additional spectra、timing、mesh quality、baseline provenance | — |

### 6.4 主图/表建议（8 图 + 3 表）

**Fig. 1** Geometry and quotient：raw grid、edge midpoints、B faces、W site、gauge orbit。  
**Fig. 2** Blind modes：gauge、`z_M`、`ip`、`jp` 的示意 + flat kernel 维数。  
**Fig. 3** Completion：`X_M` 的 polar/connection/moment 流程；`X_W^rel` 的 normal compatibility + polar transport；`t` vs `t³` 标度。  
**Fig. 4** Geometric convergence：4 曲面 `b_W` vs Dellinger，二阶斜率。  
**Fig. 5** Mechanism：flat nullity、M3 spectrum、B vs B+X_M 载荷响应、X_W 小修正。  
**Fig. 6** Linear benchmarks + state count：plate/roof/M3 的 setup、变形、error vs raw states；同时标出现代 roof 参考。  
**Fig. 7** Nonlinear benchmarks：cantilever 360°、hemisphere 全路径、open cylinder pull-out；每列 setup + 变形 + 曲线。  
**Fig. 8** Implementation + distortion：boundary storage、augmented nnz、n=128；distortion 的 plate/roof/M3，突出 M3 limitation。  
**Table 1** 方法对比：state variables、continuity、locking treatment、BC、reference。  
**Table 2** 最终 12 行 benchmark 对比：response、error（含 modern roof）、raw/active states。  
**Table 3** Scaling：n=16…128 的 raw/active、local/augmented nnz、assembly/factor、RSS；明确 backend caveat。

### 6.5 必须补的实验/工程（按优先级）

#### P0：不做会被审稿人直接问倒

1. **Clean release**
   - commit 所有修改；给 `checkerboard_shell` 打 tag（如 `v0.2-cmame`）；
   - 把 `experiments/cmame_round2/3` 纳入版本控制；
   - 更新 README：sparse/rank-one、`adaptive_continuation`、benchmark 复现命令、环境。
   - 产出 `paper/` 或 `reproduce/` 目录：一条命令生成所有图表；artifact manifest（脚本 ↔ 图 ↔ 数据 ↔ commit）。

2. **Full 3V plate**
   - 用 production `paper_components` 跑完整 3V simply supported plate（w=0，面内约束只去掉刚体）；
   - 验证 center deflection 与标量横向分支一致（线性解耦应一致）；
   - 报告 raw 3V 与 active all-DOF；这样和 MidedgeTan 的 active states 才是同口径。
   - 如果完整 3V 的 in-plane 约束引入边界层，明确说明并给对比。

3. **统一 MidedgeTan 对比**
   - 用一个 adapter、同一套 vertices/load/BC/observable，重跑 plate n=8/16/32/64、roof n=12/16/24/32、M3 n=6/8/12/16；
   - 解决 M3 的 observable 差异（Checkerboard 物理中点矩重构 vs MT raw vertex）和 BC 执行差异；
   - 解释 MT M3 的 3.741% → 5.009% 非单调：是 MT 自身收敛行为、网格差异、还是 BC/observable 差异；
   - 如果无法统一，就在表格里明确标注“different discrete observables”，不要写“identical”。
   - 现在最终表的 12 行 MT 数值来自 recovered historical rerun，只有三个网格点被新 adapter gate 过；P0 要全部由 adapter 重跑，并直接读出 raw/active states。

4. **把非线性 benchmark 移植到 production 代码**
   - cantilever（端弯矩）、hemisphere、open cylinder 至少各跑一个网格；
   - 用 `adaptive_continuation` + `LowRankTangent`，验证与 arXiv 存档结果一致（位移路径、误差、能量分块）；
   - 这一步同时验证“精确 sparse+rank-one 切线在完整非线性路径上可用”，是 CMAME 实现 section 的核心证据。
   - 如果时间有限，先做 cantilever + hemisphere；open cylinder 可作为补充。

5. **Roof 参考值透明化**
   - 表里同时给 0.3024（historical）和 0.3005924566（modern）；
   - 讨论 Checkerboard 的 roof 收敛偏慢（n=32 仍 3.16% vs modern）；
   - 不要只用 0.3024 计算误差。

6. **Gauge API**
   - 暴露 `rigid_and_gauge_basis` / `orthogonal_complement` 或等价的 quotient reduction；
   - 文档说明：物理 BC 不一定去掉 raw gauge；benchmark 里的 gauge slice 是必要的；
   - 在示例里演示 gauge 不变性测试。

#### P1：显著增强说服力

7. **MITC4 baseline**
   - 用一个 provenance-closed 的 Simo-type MITC4（商业代码或验证过的开源 FE）在同一 mesh/load/BC/observable 上跑 plate/roof/M3；
   - 如果做不到，明确写“only MidedgeTan is a verified external baseline; MITC4 is future work”，不要用文献数字假装自己跑的。

8. **同 backend 性能**
   - 在一个简单的 Python/NumPy 稀疏 FE baseline 上比 assembly/factor/solve；
   - 或者至少把 Python/AD 的 assembly、factorization、rank-one update 分项计时，避免跨语言比较；
   - 明确说“state-count economy is established; wall-clock superiority is not”。

9. **畸变机理**
   - 分析 M3 δ=.30 退化的原因：面积比 1.55、最大物理 aspect ratio 2.98、B/W frame 条件数、边界层；
   - 给 per-case 的 scaled Jacobian/aspect ratio/area ratio；
   - 如果可能，提出一个 local-step 或 weighting 修正；如果不行，明确 scope。

10. **Variable-step / unstructured**
    - 当前公式固定 `hx,hy`；如果审稿人问非均匀网格，要么推导 local weights，要么在 Introduction/Conclusion 明确“structured constant-step scope, extension is future work”。

#### P2：加分项

11. **Pinched cylinder**：如果能用 arc-length 跑通，加入；否则明确排除并说明 solver limitation。
12. **更多谱/BC**：free-edge、curved patch、baseline spectra；放 supplement。
13. **Thickness locking study**：更多 `t/h` 比，和 MITC4/IGA 对比；这是 CMAME 审稿人很可能会问的。
14. **Multi-patch / extraordinary vertex**：未来工作，不必在本文完成。

### 6.6 写作策略

- **先写“证据表”再写正文**：把每个 figure/table 对应到 claim，避免先写漂亮句子再找数据。
- **Introduction 用“方法谱系表”**：position-only vs director-based vs MITC vs IGA vs midedge；一眼看出差异。
- **Theory section 要克制**：只把真正证明的写成 theorem；其余写成 “numerical observation” 或 “design principle”。
- **Numerical section 要“先验证、后比较”**：几何收敛 → 机理 → 线性 → 非线性 → scaling → distortion。先证明方法正确，再谈优势。
- **Limitations 要主动写**：M3 畸变、roof 参考、plate 标量分支、跨 backend、无 MITC4。主动写比被审稿人挖出来好。
- **Abstract 只写四件事**：position-only 3V；blind-mode 分类与补全；sparse+rank-one；在标准 benchmark 上以约一半 raw states 达到同阶精度。不要写“robust”“efficient”“locking-free”这种没有限定词的词。

### 6.7 建议的 Abstract 骨架（英文）

> We present a position-only thin-shell discretization built on checkerboard edge-midpoint geometry. Because the four edge midpoints of any spatial quadrilateral form an exactly planar Varignon parallelogram, a local tangent frame and normal are determined by nodal positions alone, without independent director or rotation variables; the shell state is only `3V` scalar positions. Starting from the connected midpoint surface as the physical midsurface, we evaluate the first fundamental form on planar B faces and construct a W-centered second fundamental form from neighboring B-face normals. A quotient-space analysis separates the exact checkerboard gauge from genuine physical blind modes: the B metric misses one membrane direction, and the symmetric W curvature misses two curvature directions. We introduce a quotient-minimal membrane coordinate `X_M` and an objective reference-relative curvature coordinate `X_W^rel`, placed at the `O(t)` and `O(t³)` scales. The exact nonlinear tangent is a bounded-stencil sparse local operator plus one global rank-one update, solved through a bordered system without forming a dense outer product. On the recovered plate, Scordelis–Lo, and M3 benchmark sequences, the method uses about half the raw scalar states of a verified midedge baseline while achieving errors of comparable order; three finite-deformation benchmarks exercise large rotation, curved-shell pinching, and open-section ovalization. We also report the limitations: surface-preserving redistribution can materially degrade the M3 response, and no general curved-shell locking-free convergence is claimed.

（这段是骨架，具体数字和措辞按最终实验更新。）

---

## 7. 建议的下一步

### 7.1 如果目标是“尽快有一版可投的 CMAME 稿”

推荐顺序：

1. **本周**：clean release + 更新 README + artifact manifest；把现有 Round 3 的最终表/图固定下来。
2. **1–2 周**：P0-2（full 3V plate）、P0-3（统一 MidedgeTan 对比）、P0-5（roof 双参考）。
3. **2–4 周**：P0-4（非线性 benchmark 移植）；这一步最重，但也是 CMAME 实现 section 的关键。
4. **4–6 周**：写正文；同时做 P1-7（MITC4）如果可行。
5. **投稿前**：内部 review，按 Round 3 manuscript map 的 claim disposition 逐条核对。

### 7.2 如果目标是“先写一版完整草稿，再补实验”

可以先按第 6.3 的结构写，把 P0 实验留成 `TODO` 占位；但**不要**把旧 Round 1 的重构驱动结果写进去。写作时以 Round 3 的 12 行表和 arXiv 非线性结果为基础，标注哪些是“archived production data”、哪些是“new production rerun”。

### 7.3 现在最值得先做的一件事

**把非线性 benchmark 移植到当前 production 代码，并跑通 cantilever + hemisphere。**

理由：
- 它是当前最大的证据缺口（arXiv 非线性结果不在新代码里）；
- 它能同时验证 `LowRankTangent` 在完整载荷路径上的可用性；
- 它是 CMAME “implementation + validation” 故事的核心；
- 一旦跑通，scaling、distortion、gauge 的证据链就完整了。

---

## 附录 A：关键文件索引

### arXiv 稿与原始仓库

- `[/home/song/projects/CB-shell/arxiv_submission/Checkerboard_Shells.tex](/home/song/projects/CB-shell/arxiv_submission/Checkerboard_Shells.tex)`：arXiv 英文稿（1097 行）。
- `/home/song/projects/CB-shell/paper_draft_cmame/planning/manuscript_outline.md`：中文成稿骨架。
- `/home/song/projects/CB-shell/paper_draft_cmame/planning/material_audit.md`：公式/证据等级审计。
- `/home/song/projects/CB-shell/outputs/checkerboard_nonlinear_benchmarks/`：非线性 benchmark 原始输出。
- `/home/song/projects/CB-shell/experiments/checkerboard_nonlinear_benchmarks/`：非线性驱动。

### 当前工作区代码

- `src/checkerboard_shell/geometry/core.py`：B/W 几何。
- `src/checkerboard_shell/geometry/boundary.py`：边界闭合、`X_W^rel`、局部 stencil。
- `src/checkerboard_shell/geometry/intrinsic.py`：`X_M`。
- `src/checkerboard_shell/mechanics/energy.py`：能量分块。
- `src/checkerboard_shell/mechanics/low_rank.py`：精确 sparse+rank-one。
- `src/checkerboard_shell/mechanics/model.py`：`CheckerboardShell`。
- `src/checkerboard_shell/solver/newton.py`：dense Newton。
- `src/checkerboard_shell/solver/nonlinear.py`：通用 Newton-GMRES/CG。
- `src/checkerboard_shell/solver/continuation.py`：自适应 continuation。
- `tests/`：19 tests。

### CMAME 实验

- `experiments/cmame_additional/CMAME_NUMERICAL_REPORT.md`：Round 1 审计（**历史，不作为正文证据**）。
- `experiments/cmame_round2/CMAME_ROUND2_REPORT.md`：驱动恢复、稀疏/秩一、scaling、畸变。
- `experiments/cmame_round2/reports/formulation_scope.md`：constant-h 适用范围。
- `experiments/cmame_round2/reports/xm_derivative_structure.md`：`X_M` 导数结构。
- `experiments/cmame_round3/CMAME_ROUND3_REPORT.md`：β 复核、最终对比、最强 claim。
- `experiments/cmame_round3/04_manuscript_map/CMAME_ROUND3_MANUSCRIPT_MAP.md`：逐节改稿 checklist。
- `experiments/cmame_round3/02_final_comparison/MANUSCRIPT_COMPARISON_TABLE.md`：12 行最终表。
- `experiments/cmame_round3/03_final_figures/FIGURE_GUIDE.md`：每张图能支持/不能支持什么。

---

## 附录 B：复现命令

```bash
# 环境
cd /home/song/projects/checkerboard_shell
PY=/home/song/miniconda3/envs/allatom/bin/python

# 测试
PYTHONPATH=src $PY -m pytest -q

# plate demo
PYTHONPATH=src MPLBACKEND=Agg $PY examples/plate_demo.py

# Round 2（按需）
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_driver_recovery
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_structure_audits
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_boundary_equivalence
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_exactness
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_scaling
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_distortion
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_midedgetan_distortion
PYTHONPATH=src:. $PY -m experiments.cmame_round2.run_regressions
PYTHONPATH=src:. $PY -m experiments.cmame_round2.make_figures
PYTHONPATH=src:. $PY -m experiments.cmame_round2.finalize_artifacts

# Round 3
PYTHONPATH=src:. $PY -m experiments.cmame_round3.run_beta_revalidated
PYTHONPATH=src:. $PY -m experiments.cmame_round3.build_final_comparison
PYTHONPATH=src:. $PY -m experiments.cmame_round3.collect_final_figures
PYTHONPATH=src:. $PY -m experiments.cmame_round3.finalize_round3
```

注意：MidedgeTan 的 distortion 构建需要脚本里记录的 LibShell/Eigen 源码根目录；可用 `LIBSHELL_SOURCE_ROOT` / `EIGEN_SOURCE_ROOT` 覆盖。

---

## 8. Round 4 补实验进展（2026-10-26）

新增目录 `experiments/cmame_round4/`，用于补齐本文档第 6.5 节的 P0 实验；不修改任何公式、系数或边界闭合。完整报告见 [CMAME_ROUND4_REPORT.md](../experiments/cmame_round4/CMAME_ROUND4_REPORT.md)。

已完成：

- **完整 3V plate**：`w=0` 四边、面内仅三个刚体 pin，production `paper_components` + `sparse_low_rank`；n=8/16/32/64 的结果、active states 和与 scalar 分支的 cross-check 在 `experiments/cmame_round4/CMAME_ROUND4_REPORT.md`。
- **统一 12 行 MidedgeTan 对比**：从 pinned LibShell commit 重建 adapter；所有 MT 行用同一套 vertices/load/BC 重跑，raw/active states 为 measured；roof 同时给 0.3024 与 0.3005924566 两个参考；M3 主值统一为 raw equator vertex，physical-midpoint 重构值作为第二观测（n=16: raw 0.4795%，midpoint 0.1105%，MT 5.009%）。
- **cantilever**：16x4/24x4/32x4 全部 VALIDATED；32x4 为 359.841°、RMS 2.442%、max 2.792%，与 arXiv 一致。
- **hemisphere**：n=6/8/12/16 全部重跑；n=16 为 A 2.759%、B 0.149%、joint RMS 1.736%，与 arXiv 一致。
- **open cylinder**：8x8/12x12/16x16 全部重跑；16x16 combined RMS 4.483%，与 arXiv 一致。
- **production rank-one 切线**：在 cantilever 32x4、hemisphere n=12、open cylinder 8x8 的 archived 末态上，reduced `K_local + αggᵀ` 与 dense AD Hessian 的相对差分别为 1.64e-16、1.49e-15、4.95e-15。
- **M3 双观测**：主值统一为 raw equator vertex（与 MT 相同），physical-midpoint 重构值作为第二观测。n=16: raw 0.4795%，midpoint 0.1105%，MT 5.009%；结论不变，但 headline 不再口径错位。
- **MITC4 baseline**：OpenSeesPy 3.7.1 `ShellMITC4` + `ElasticMembranePlateSection`（6 DOF/node），同网格/载荷/物理 BC。plate 在 n=16 为 0.073%，但 t/a=0.01 有 ~0.5% 的 Mindlin 剪切 floor（t=0.001 时降到 0.001%）；roof 收敛到 0.3009（modern 参考 0.30059，误差 0.08–0.12%）；M3 因 raw corner point load 收敛慢，n=32 为 0.47%。**MITC4 对比是 mixed：不能声称全面优于 MITC4，但可声称 3V 状态数和特定 case 的优势。**

仍需处理：

- 通用 gauge quotient API、clean release/tag、README 更新。
- M3 的 load smoothing（Checkerboard 的物理点点算子与 MT/MITC4 的 raw 点载荷仍不同；当前 MITC4 raw/smoothed 结果都显示 coarse mesh 收敛慢）。

数据与命令见 `experiments/cmame_round4/README.md`。

---

## 9. Round 5 设计应用算例（2026-10-26）

新增 `experiments/cmame_round5/`，四个 bounded design demos，不修改 formulation：

1. **Adjoint 梯度验证**：9 个节点力设计变量、297 个自由 DOF、99 个观测；adjoint 与中心差分最大相对差 `1.5e-12`。
2. **驱动形状控制**：9 个节点力，目标位移相对误差 `6.5e-5`；力分布不唯一（相对差 `0.174`），因此 claim 是 shape control，不是 force identification；13 次 L-BFGS。
3. **参考形状优化**：12 个高度控制点，目标位移相对误差 `1.6e-4`；高度场非唯一（diagnostic `0.97`），因此主指标是位移，不是高度恢复；20 次迭代（iteration cap）。
4. **大转动形状设计**：2 个初始高度、46 个自由 DOF，`M/Mmax=0.25`；目标初始曲率 `h_true=[0,0.375,1.5]` 恢复到 `1e-6`，tip error `1.9e-6`；8 次迭代；无转动设计变量。

这些 demo 支撑：

- 3V 自由状态 → 3V adjoint；
- exact sparse+rank-one tangent → adjoint 梯度到 `1e-12`；
- 设计变量就是节点位置/力，不需要转动设计变量；
- 大转动下位置变量设计可精确恢复目标初始曲率。

不支撑：

- 比 MITC4 更准/更快；
- 驱动力或参考高度的唯一性；
- FPC/contact/batch 应用。

文件：`experiments/cmame_round5/DESIGN_REPORT.md`、`DESIGN_DEMO_TABLE.md`、`ADJOINT_CHECK.md`。

### 9.1 同目标设计循环对比（应用层）

同一个 9 作动器目标形状问题，同一网格/BC/观测，三个层次：

| 方法 | 设计变量 | 自由 DOF | 迭代 | 工作量 | 总时间 | 目标误差 |
|---|---:|---:|---:|---:|---:|---:|
| Checkerboard + adjoint | 9 | 297 | 13 | 79 次线性求解 | **1.05 s** | **0.0065%** |
| Checkerboard + finite difference | 9 | 297 | 12 | 450 次线性求解 | 1.96 s | 0.0066% |
| MITC4 + finite difference | 9 | 858 | 13 | 374 次 forward 分析 | 12.83 s | 0.0208% |

- **层 1（同方法 adjoint vs FD）**：79 vs 450 次求解，少 `5.7x`；优化时间 0.20 vs 1.11 s，快 `5.6x`；总时间 1.05 vs 1.96 s，因为 0.86 s 的切线组装是共同成本。
- **层 2（状态规模）**：297 vs 858 自由 DOF，约 `1/3`。
- **层 3（black-box MITC4 FD）**：12.83 s vs 1.05 s，约 `12.2x`。这是实现层对比：MITC4 每次分析重建模型且用 FD；如果给 MITC4 写 adjoint，它的工作量会降下来。不能据此说 MITC4 本身更慢。
