# Checkerboard 非线性求解器迁移报告

## 1. 修改和新增文件

安装包代码：

- `src/checkerboard_shell/solver/config.py`（通用求解器参数与边界校验）；
- `src/checkerboard_shell/solver/convergence.py`（残量范数与收敛监控器）；
- `src/checkerboard_shell/solver/nonlinear.py`（抽象问题接口以及 HVP Newton--GMRES/CG）；
- `src/checkerboard_shell/solver/continuation.py`（自适应标量续接）；
- `src/checkerboard_shell/solver/__init__.py`（solver 公共 API 导出）；
- `src/checkerboard_shell/__init__.py`（安装包顶层公共 API 导出）；
- `tests/test_nonlinear_solver.py`（通用二次问题与 solver 回归测试）。

文档：

- `docs/solver_extraction_analysis.md`（提取前的现状、边界和迁移分析）；
- `docs/solver_architecture.md`（最终模块、参数和调用架构）；
- `docs/solver_migration_report.md`（本文档，记录迁移内容和验证结果）。

FPC 专用适配器继续保留在 `FPC/experiments/checkerboard_solver_adapter.py`（FPC 应用适配层）中，没有把任何 FPC 模块复制到安装包。Checkerboard geometry 和 energy 文件均未修改。

## 2. 已迁移内容

安装包现在提供以下通用能力：

- 最小 `NonlinearProblem` 接口；
- residual、gradient 和 HVP 的自动微分默认实现；
- 可配置的矩阵无关 HVP Newton--GMRES；
- 带负曲率截断的 Newton--CG；
- 基于能量、残量或混合 merit 的线搜索全局化；
- 残量与能量收敛历史和计算量计数；
- 带 rollback、步长缩小和步长增长的自适应 continuation；
- 通用的可调用预条件器接口。

本次还补齐了以下接口行为：

- 应用只覆写 `residual()` 而不覆写 `hvp()` 时，`linearize()` 会从应用提供的残量计算图构造 HVP，不会绕过自定义残量；
- Newton--CG 会实际使用传入的通用预条件器，并采用负曲率阈值进行截断；
- `solve_nonlinear` 会识别 `SolverConfig.continuation=True` 并进入自适应 continuation；
- `NonlinearResult.krylov_iterations` 提供模型无关的 Krylov 计数别名；
- 最终能量严格在返回坐标处重新计算；
- 无效的 Krylov、damping、线搜索和 continuation 参数会在配置构造时被拒绝。

## 3. 明确未迁移的内容

以下内容仍然留在实验或应用层：

- trust-region 实现；
- pseudo arc-length；
- neural initialization；
- DDG-aware preconditioner；
- FPC 载荷、材料和诊断；
- benchmark、plot、CSV/JSON 写出和 case study。

`trust_region=True` 目前只是预留配置，会明确报错，而不会静默改变求解算法。

## 4. 验证结果

验证使用现有 `allatom` Python 环境完成。

### 4.1 安装包测试

- `tests/`（安装包测试目录）共 16 项测试，全部通过；
- 原有 forward tests 继续通过；
- 新测试覆盖默认 HVP、自定义 HVP、自定义 residual graph、Newton--GMRES、预条件 Newton--CG、绝对/相对收敛判据、自适应 continuation 和无效配置拒绝。

### 4.2 FPC 原型回归

- `FPC/experiments/test_experiments.py`（FPC 原型与实验回归测试）共 15 项测试，全部通过。
- 已重新运行 `FPC/experiments/validate_solver_migration.py`（迁移前后回归脚本）。
- 回归汇总保存在 `FPC/results/solver_migration/regression.json`（求解器迁移回归结果）中。

主要对照结果如下：

| 案例 | 状态 | 主要对照结果 |
|---|---|---|
| FPC case B | 迁移前后均收敛 | 两条路径均为 7 次 continuation 尝试和 20 次 Newton；残量差 `1.82e-11 N`，能量差 `2.50e-16 J`，最大位移差 `7.40e-12 m`，坐标差 `1.12e-11 m` |
| nonlinear benchmark 完整路径 | 迁移前后均收敛 | 原型与提取路径分别为 16 和 18 次尝试、31 和 30 次 Newton；坐标差 `1.97e-3 m`，能量差 `7.13e-9 J` |

FPC case B 为提取前后等价性提供了有力证据。另一次单载荷级 nonlinear corrector 重放使用完全相同的载荷级、容差、初始状态和 continuation 状态，得到：

- 两条路径均为 6 次 Newton；
- 两条路径均为 133 次 Krylov；
- 坐标差 `3.75e-12 m`；
- 能量差 `1.36e-17 J`。

完整 nonlinear benchmark 不是严格的数学一致性对照，因为其 legacy 和 extracted 路径使用了不同的 geometry-aware 与 generic continuation 策略。两条路径都到达单位载荷并满足指定残量容差；最终状态差异来自验证脚本中的载荷路径控制差异，不能据此认为 DDG formulation 发生了变化。

## 5. 当前限制

- 核心 continuation 驱动器只跟踪从 0 到 1 的归一化标量参数，具体应用负责将它映射为实际物理载荷。
- trust-region 和 pseudo arc-length 目前只是未来扩展方向，不是生产级核心方法。
- Newton--CG 主要用于诊断和参考；遇到不定 shell tangent 时可能提前截断。
- 若要严格比较完整 nonlinear benchmark 路径，还需让迁移前后实现使用完全相同的 continuation 策略进行重放。
- 原有 dense `newton.py`（小规模 dense Newton 参考求解器）继续用于向后兼容和小系统验证，它没有被矩阵无关实现替换。

## 6. 结论

本次迁移只提取了已经验证、具有通用意义的 nonlinear solver infrastructure。最终安装包具备：

```text
Checkerboard DDG formulation
        +
通用 nonlinear solver framework
```

安装包没有演变成 FPC simulation package。Checkerboard 数学模型、position-only 状态、能量定义和原有 forward API 均保持不变。

