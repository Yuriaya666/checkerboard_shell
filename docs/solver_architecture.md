# Checkerboard Shell 求解器架构

安装包现在由原有 Checkerboard DDG formulation 和模型无关的非线性平衡求解框架组成。求解器不会直接计算 Checkerboard energy；具体应用通过实现 `NonlinearProblem`（通用非线性问题接口）提供所需算子。

## 1. 安装包结构

```text
checkerboard_shell/src/checkerboard_shell/solver/
    __init__.py       公共 API 导出
    config.py         SolverConfig 配置
    convergence.py    残量范数和 ConvergenceMonitor
    nonlinear.py      NonlinearProblem 和 Newton--Krylov corrector
    continuation.py   自适应标量 continuation
    newton.py         原有 dense Newton 参考 API
    result.py         原有正向仿真结果对象
```

`newton.py`（dense Newton 参考求解器）继续作为现有正向 API 使用的小规模 dense reference path。新提取的核心位于 `nonlinear.py`（矩阵无关非线性求解器）和 `continuation.py`（自适应续接驱动器）。

## 2. Problem 接口

```python
from checkerboard_shell.solver import NonlinearProblem

class MyProblem(NonlinearProblem):
    def energy(self, x): ...
    def residual(self, x, *, create_graph=False): ...  # 可选覆写
    def hvp(self, x, v): ...                            # 可选覆写
    def update_parameter(self, value): ...              # 可选的续接参数更新
```

基类默认从 `energy` 自动微分得到 `residual`，并将 `gradient` 作为 `residual` 的别名，同时提供自动微分 HVP。拥有解析切线作用或缓存切线作用的模型可以覆写 `hvp`。`update_parameter` 默认不执行操作；`residual_scale` 是供相对残量判据使用的可选接口。

该接口不会假定载荷类型、材料、几何、单位或参数名称。FPC 适配层可以把自身的载荷因子映射到 `update_parameter`，但适配层必须留在安装包之外。

## 3. 非线性 corrector

```python
from checkerboard_shell.solver import SolverConfig, solve_newton_gmres

config = SolverConfig(
    method="newton-gmres",
    max_iterations=20,
    max_krylov_iterations=50,
    tolerance=1e-8,
    line_search=True,
    line_search_merit="mixed",
)
result = solve_newton_gmres(problem, initial, config)
```

底层入口 `solve_level` 始终只求解一个参数值。模型无关的 `solve`/`solve_nonlinear` dispatcher 会识别 `SolverConfig(continuation=True)` 并转入 `adaptive_continuation`；如果只需要单次 corrector，应显式调用 `solve_newton_gmres` 或 `solve_newton_cg`。

每个 Newton 状态只建立一次残量计算图，切线仅通过 HVP 作用于向量。由于壳模型的切线可能不定，默认方法是 GMRES。对于 SPD 区间的诊断，可以调用截断 Newton--CG：

```python
result = solve_newton_cg(problem, initial, config)
```

通用 `solve` 函数可以接收可选的 `preconditioner(vector)` 回调，但安装包不会构造任何 DDG-aware 预条件器。

`NonlinearResult`（非线性求解结果）记录：

- 求解状态和返回坐标；
- nonlinear 与 Krylov 迭代次数；
- HVP 和梯度计算次数；
- 线搜索拒绝次数；
- 残量范数与最终能量；
- 最小接受步长；
- 每次非线性迭代和线性求解历史。

## 4. 参数接口

`SolverConfig`（统一求解器配置）包含以下算法参数：

- `method`：`"newton-gmres"` 或 `"newton-cg"`；
- `max_iterations`、`max_krylov_iterations`：非线性和 Krylov 最大迭代次数；
- `tolerance`、`relative_tolerance`：绝对和可选相对残量容差；
- `gmres_tolerance`、`forcing_floor`：非精确 Krylov forcing 控制；
- `damping`：可选切线阻尼；
- `line_search`：是否启用线搜索；
- `line_search_merit`：`"energy"`、`"residual"` 或 `"mixed"`；
- `max_backtracks`、`armijo`：回溯次数和 Armijo 参数；
- `continuation`：是否由通用入口启用自适应 continuation；
- `initial_load_step`、`min_load_step`、`max_load_step`：初始、最小和最大参数步长；
- `step_growth`、`step_shrink`：步长增长和缩小系数；
- `easy_newton_iterations`、`easy_alpha`：判断 corrector 是否容易求解的阈值；
- `max_continuation_attempts`：最大 continuation 尝试次数；
- `trust_region`：预留开关。目前设置为 `True` 时会明确抛出 `NotImplementedError`，不会假装启用尚未验证的核心算法。

配置对象会在构造时拒绝无效的容差、负 damping、不合法的线搜索参数以及相互矛盾的 continuation 参数。

## 5. 收敛监控与 continuation

`convergence.py`（收敛判据模块）提供：

- 残量无穷范数和 L2 范数；
- 绝对或相对残量收敛判据；
- 只追加记录的 `ConvergenceMonitor`（收敛监控器）。

`adaptive_continuation` 将归一化标量参数从 0 推进到 1：

```python
from checkerboard_shell.solver import SolverConfig, adaptive_continuation

config = SolverConfig(
    continuation=True,
    initial_load_step=0.1,
    min_load_step=0.01,
    max_load_step=0.2,
)
path = adaptive_continuation(problem, initial, config)
```

corrector 失败时，驱动器保持上一个已接受坐标不变，缩小步长后重试。连续出现容易求解的 corrector 时，步长会在配置上限内增大。如果最后一次尝试失败，驱动器会恢复最后一个已接受参数，并返回全部 attempt 和 corrector 历史。

## 6. 不同问题的调用方式

### 6.1 直接能量问题

如果问题只有能量表达式，只需实现 `energy(x)`；基类会自动生成残量和 HVP：

```python
class EnergyProblem(NonlinearProblem):
    def energy(self, x):
        return model_energy(x)

result = solve_newton_gmres(
    EnergyProblem(), initial,
    SolverConfig(tolerance=1e-8),
)
```

### 6.2 自定义残量或 HVP 的问题

如果应用已有经过优化的 residual 或 tangent action，可以覆写对应方法。求解器会使用应用提供的残量图或 HVP，而不需要知道这些算子的物理来源。

### 6.3 参数化载荷问题

参数化问题通过 `update_parameter(value)` 接收 continuation 参数。该参数可以由应用解释成载荷系数、边界控制量或其他标量，但核心求解器不会对其赋予物理含义。

## 7. 依赖边界

核心 solver 只依赖 Python、PyTorch 以及自身的配置和收敛模块。它不会导入 `FPC`、`experiments`、Lorentz force、conductor、magnetic field、material 或 Checkerboard energy 函数。

以下内容继续保留在应用或实验层：

- benchmark、plot 和 case study；
- neural initialization；
- pseudo arc-length；
- trust-region prototype；
- DDG-aware preconditioner；
- FPC 专用载荷、材料和诊断。

