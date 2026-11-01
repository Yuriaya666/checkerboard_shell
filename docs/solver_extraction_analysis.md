# Checkerboard DDG 非线性求解器提取分析

## 1. 分析范围与审阅资料

本次分析的目标，是在不改变已冻结的 Checkerboard position-only 能量

\[
E(q)=E_M(a_B(q))+X_M(q)+E_B(b_W(q))+X_W(q)
\]

的前提下，将可复用的非线性算法与 FPC 应用分离。

本次审阅覆盖了当前 `checkerboard_shell` 安装包、`FPC/experiments` 下的求解器原型与适配层，以及 `FPC/docs` 下的求解器报告。重点对照了：

- `hvp_newton.py`（矩阵无关 corrector、Krylov 方法、线搜索和 DDG 预条件器实验）；
- `adaptive_continuation.py`（FPC 载荷驱动以及回滚、增步策略）；
- `trust_region_newton.py`（坐标范数和几何范数 trust-region 研究）；
- `solvers.py`（早期 dense、CG 和 GMRES 实现）；
- `nonlinear_benchmark.py`、`fpc_cases.py` 和 `validate_solver_migration.py`（benchmark、应用调用路径和回归判据）；
- 最终求解器报告、设计空间报告和提取计划。

安装包原本已经在 `solver/newton.py`（dense Newton 参考求解器）中提供 dense reference Newton。审阅开始时，工作区还存在一份尚未提交的提取草稿，包括 `config.py`（求解器配置）、`convergence.py`（收敛监控）、`nonlinear.py`（非线性求解器）、`continuation.py`（自适应续接）和对应的单元测试。本次工作将这些内容视为需要审查和完善的已有工作，而不是进行大规模重写的依据。

## 2. 当前已有的求解器组件

### 2.1 安装包原有求解器

`checkerboard_shell.solver.newton` 提供了一条显式构造 dense Hessian 的参考求解路径，包括：

- 通过自动微分计算能量、梯度和完整 Hessian；
- 对称特征分解以及经过特征值过滤的 Newton 方向；
- 基于能量或残量的线搜索接受判据；
- `CheckerboardShell.solve()` 中的固定载荷级序列；
- 面向现有正向仿真公共 API 的紧凑结果对象。

这套实现适合小规模问题和 formulation 回归验证，但其 dense Hessian 不适合作为可扩展的通用非线性求解框架。本次提取保留其现有公共 API，不对其进行重写。

### 2.2 FPC 原型中的组件

FPC 研究代码中包含以下求解器能力：

1. 矩阵无关 Newton corrector：每次非线性迭代只建立一次残量计算图，再通过二次 VJP 得到 Hessian-vector product（HVP）。
2. 确定性、无重启、左预条件 GMRES。
3. 带负曲率处理的截断 HVP-CG。
4. 基于残量、能量或二者混合的线搜索 merit，并在必要时回退到残量下降方向。
5. 类 Eisenstat--Walker 的非精确 Krylov forcing 策略。
6. 绝对残量收敛判据，以及 nonlinear、Krylov、HVP、梯度评估和线搜索拒绝次数等诊断计数。
7. 载荷续接：失败时回滚到上一个已接受状态并缩小步长，连续出现容易求解的 corrector 后增大步长。
8. 面向 DDG 的材料缩放、block-Jacobi、曲率加权和 local-Krylov 预条件器。
9. trust-region Newton--GMRES 原型：支持坐标范数或 Checkerboard 几何范数、基于模型下降量的接受判据和 Cauchy fallback。
10. 实验性的几何感知 continuation、`X_M` 低秩修正、神经网络初值以及尚在规划中的 pseudo arc-length 机制。

FPC 原型中没有单一、模型无关的 convergence monitor 类。其收敛数据主要以求解器或案例专用的列表和 CSV 行累积。真正可复用的是残量范数、收敛判据和记录行为，而不是具体 CSV 字段格式。

## 3. 可通用化的非线性求解器基础设施

以下组件只依赖抽象的可微非线性问题，适合进入安装包：

- **HVP Newton--GMRES**：只需要残量和切线作用，不依赖 Checkerboard 几何或材料模型，并且无需组装全局 Hessian。
- **Newton--CG**：适合局部 SPD 切线，也可作为参考方法。由于壳模型切线可能不定，CG 不应取代 GMRES 成为默认方法。
- **线搜索全局化**：能量、残量或混合接受判据都可以完全通过抽象问题接口表达。
- **收敛检查与监控**：绝对/相对残量范数、迭代历史和操作计数与具体问题无关。
- **自适应 continuation**：标量参数更新、已接受状态回滚和步长调整都具有通用意义。安装包中的参数不应被命名或解释为电流、压力、磁场或 Lorentz 载荷。
- **可选预条件器回调**：非线性核心可以接收一个可调用的线性预条件器，但其构造方法和物理含义必须留在核心之外。

## 4. FPC 专用或验证尚不充分的组件

以下内容应继续保留在 `FPC/experiments` 或其他应用、研究层：

- `FPCProblem`、导体路径、Lorentz 积分、磁场、电流、FPC 材料常数、边界案例以及 FPC 专用结果格式；
- 膜/弯曲坐标缩放，以及所有需要读取网格步长、材料厚度、曲率、B/W stencil 或 `X_M` 的预条件器；
- 根据曲率变化、力的单位或案例名称作决定的 FPC 专用 continuation 策略；
- 几何范数 trust-region，以及配套的 Checkerboard metric/curvature 安全诊断；
- benchmark 调度、绘图、CSV/JSON 输出和 case study 脚本；
- 神经网络初值和训练数据生成。

普通坐标范数的 trust-region 算法比几何感知版本更通用，但现有报告显示：在已经验证、没有 fold 的载荷路径上，它增加了成本且没有带来收益。因此，本次提取仍将其保留为 experimental。`SolverConfig.trust_region` 可以作为预留开关，但遇到尚未支持的配置时必须明确报错，不能静默切换到未经充分验证的算法。

Pseudo arc-length 目前也不足以进入核心。代码不应提供一个看似可用、实际却没有完成分支追踪的占位实现；只需保证未来可以扩展 continuation API，而不是现在就把它暴露为可用方法。

## 5. 本次应迁移的组件

最小安装包提取应包含：

1. `NonlinearProblem`（通用非线性问题接口），提供 `energy(x)`、`residual(x)`、`gradient(x)`、`hvp(x, v)` 和 `update_parameter(parameter)`。基类可以提供自动微分默认实现，具体问题也可以覆写切线作用。
2. `SolverConfig`（统一求解器配置），管理非线性方法、迭代与容差、Krylov forcing、damping、线搜索、标量 continuation 和预留的 trust-region 开关。
3. 共用同一套非线性循环与线搜索策略的矩阵无关 Newton--GMRES 和截断 Newton--CG。
4. 与模型无关的残量范数、收敛判据和只追加不修改的监控记录。
5. 带回滚、步长缩小/增长限制以及完整 attempt/corrector 诊断的标量自适应 continuation。
6. 位于安装包之外的 FPC 专用适配器，用于把 `FPCProblem.potential` 及其载荷因子映射到 `NonlinearProblem`。

原有 dense 求解器继续保留，以维持向后兼容，但不重写，也不把它描述为新矩阵无关核心的一部分。

## 6. 本次暂不迁移的组件

- DDG-aware、material-aware 和 geometry-aware 预条件器构造；
- trust-region corrector，包括几何范数版本；
- pseudo arc-length 和 fold 检测；
- neural initialization；
- experiment、benchmark、plot、报告生成和 case study 代码；
- `checkerboard_shell` 中的 FPC 专用诊断或结果类型；
- 对 Checkerboard 几何、能量项或其权重的任何修改。

## 7. 目标设计与依赖边界

目标安装包结构如下：

```text
checkerboard_shell/solver/
    __init__.py       稳定的公共求解器导出
    config.py         只包含算法参数的配置模块
    convergence.py    残量范数、收敛判据和监控器
    nonlinear.py      Problem 接口和 Newton--Krylov corrector
    continuation.py   标量自适应 continuation
    newton.py         原有 dense 参考求解器，API 不变
```

依赖方向如下：

```text
应用/模型适配层 -> NonlinearProblem -> 非线性求解器
                                      |
                                      +-> continuation 驱动器
```

提取后的 solver 模块不得导入 FPC、experiments、Lorentz force、conductor、magnetic field、Checkerboard energy 函数或 material 类。求解器只能调用抽象问题接口。

## 8. 对已有提取草稿的审查结论

已有草稿基本符合目标模块划分，也将应用依赖隔离在安装包之外，并保留了原有 dense Newton API。不过，在完成提取之前仍需核查或小范围修正以下问题：

- 配置校验不能只检查方法名，还必须覆盖步长、forcing、damping、线搜索和 continuation 参数边界；
- 结果字段应尽量采用模型无关的命名，同时可以保留必要的兼容别名；
- corrector 被拒绝后，continuation 必须同时恢复最后一个已接受状态和对应参数；
- 即使线搜索失败或达到最大迭代次数，corrector 也应报告返回坐标处的最终能量和残量；
- 自定义 residual/HVP 的自动微分计算图语义必须明确；
- 可选预条件器接口必须保持通用，不能让使用者误以为实验性的 DDG 预条件器已经迁入；
- `continuation` 配置开关必须有明确公共语义：要么由顶层 solve dispatcher 执行，要么在文档中明确要求调用 `adaptive_continuation`。

## 9. 验证状态与验收标准

在现有 `allatom` Python 环境中，安装包测试共 16 项，全部通过；覆盖了非线性问题接口、求解器、continuation 和配置校验。已有 FPC 回归结果显示：

- FPC case B：迁移前后 continuation 尝试次数和迭代数一致；残量差约 `8.1e-11`，能量差约 `1.12e-15 J`，最大位移差约 `3.31e-11 m`，坐标差约 `5.01e-11 m`；
- nonlinear benchmark：两条路径均收敛，但由于 continuation 控制规则不同，分别进行了 16 次和 18 次尝试，最大坐标差约 `1.97e-3 m`。

FPC 结果为迁移等价性提供了有力证据。完整 nonlinear benchmark 不能作为严格的一一对照，因为其中比较了不同的步长控制策略。为此，又使用完全相同的初始状态、载荷、容差和 corrector 配置进行了单载荷级重放：迁移前后均为 6 次 Newton 和 133 次 Krylov 迭代，坐标差为 `3.75e-12 m`，能量差为 `1.36e-17 J`。以后若比较完整载荷路径，也应首先让两套实现重放完全相同的 continuation 调度或策略，避免把载荷路径差异误判为求解器差异。

