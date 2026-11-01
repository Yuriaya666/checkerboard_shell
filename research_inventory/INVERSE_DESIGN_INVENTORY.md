# 全部逆向设计实验专项清查

## 口径与结论

按**独立设计问题**计数，本项目有 **3 个非工程基础逆向案例**（I01–I03）和 **4 个工程或工程前驱案例**（I04–I07）。伴随梯度核对、FD 比较、不同网格、Round 6→7 的飞机物理参数更新，均为上述案例的实验系列，**不另算新逆问题**。外部 `4d_printing_inverse_design` 是独立双层材料/打印研究支线 X04；未逐项审计其数千个输出，不计入本篇 CB-Full 的三个非工程逆向案例。当前仓库中未找到第四个已有结果、代码和物理目标均完整的非工程逆向问题。

| ID | 非工程问题 | 设计变量与平衡 | 实际结果 | 证据状态 |
|---|---|---|---|---|
| I01 | 九作动力平板**形状控制** | 9 个 z 力；12×8，297 个 CB 自由状态、99 个观测 | 目标误差 0.0065%；力恢复误差 0.174；一致梯度与中心差分一致 | 已验证梯度和目标形状，未验证力唯一识别 |
| I02 | 12 参数参考形状 | 预弯板 4×4 控制格中 12 个高度；参考度量/曲率影响正向 | 20 次迭代达上限；位移相对误差 1.60e-4，参数误差 0.968 | 目标拟合有效，参数非唯一且优化终止受限 |
| I03 | 大转动参考几何反演 | 2 个初始高度；端矩；非线性平衡、约 90° 旋转 | 8 次迭代，tip 误差 1.86e-6，恢复高度 [0,0.374998,1.499998] | 小规模合成问题完成；最直接耦合几何参数与大转动平衡 |

### I01：平板九力逆向与梯度比较

- **研究目的：** 检验精确切线下的力变量隐式/伴随梯度、目标形状控制与完整优化成本。目标由已知力生成，因此是合成形状控制。
- **原始 Round 5 数据：** CB + adjoint 13 次优化迭代、79 次线性解、1.05 s、目标相对误差 0.0065%；同 CB 的中心 FD 为 12 次、438 次线性解、2.08 s、0.0066%；MITC4 + FD 为 17 次、533 次正向分析、11.85 s、0.0026%。MITC4 最终目标误差更低。早期比较混入 MITC4 12×12 与 CB 12×8 网格，后续等网格 FD/FD 脚本已补充；原三法表只能当**工作流**比较。
- **梯度与可辨识性：** 原报告指出伴随相对核对约 1e-12 量级；`raw_data/adjoint_check.json` 保留精确数字。9 力恢复相对误差 0.174，表明单一平滑位移目标不能唯一确定力。后期重跑给出约 68/391 次线性解以及约 1.076/1.908 s，属于同案例不同实现/计时，不应与原始 79/438 行混写。
- **程序：** [`adjoint_check.py`](../experiments/cmame_round5/adjoint_check.py)、[`actuation_design.py`](../experiments/cmame_round5/actuation_design.py)、[`fair_fd_adjoint.py`](../experiments/cmame_round5/fair_fd_adjoint.py)、[`fd_vs_fd.py`](../experiments/cmame_round5/fd_vs_fd.py)。**原始记录：** [`raw_data/`](../experiments/cmame_round5/raw_data)、[`DESIGN_REPORT.md`](../experiments/cmame_round5/DESIGN_REPORT.md)。**图：** [`Figure_R5_actuation_design.pdf`](../experiments/cmame_round5/figures/Figure_R5_actuation_design.pdf)。
- **可支持：** 同一 CB 模型上的一致梯度与 FD 数值一致，且降低线性求解工作量。**不可支持：** CB 单元固有 CPU 加速、唯一力反演、MITC4 精度更差。

### I02：参考形状的 12 参数逆向

- **目的与设定：** 通过改变预弯板的 12 个内部高度控制量影响参考度量、参考曲率和 297 自由状态的平衡；99 个观测。初始/最终目标分别为 8.58e-7/5.66e-10，位移相对误差 1.60e-4，但高度场误差 0.968；20 次迭代到上限。
- **程序/数据/图：** [`shape_design.py`](../experiments/cmame_round5/shape_design.py)、[`shape_design.json`](../experiments/cmame_round5/raw_data/shape_design.json)、[`shape_design_history.csv`](../experiments/cmame_round5/raw_data/shape_design_history.csv)、[`Figure_R5_shape_design.pdf`](../experiments/cmame_round5/figures/Figure_R5_shape_design.pdf)。
- **判断：** 它检验了与 I01 不同的**初始几何设计通路**，不只是另一个力逆向。单一位移目标下参数严重非唯一，不能宣称参考几何已正确识别；没有独立同目标传统壳反演效率基线。

### I03：约 90° 条带大转动参考几何逆向

- **目的与设定：** 两个初始中心线高度作为设计变量；端矩 (M/M_{max}=0.25)；46 个自由状态。已知目标 tip (6.342622, 8.435393)，优化 tip (6.342623, 8.435391)。真高度 [0,0.375,1.5]、优化 [0,0.374998,1.499998]，8 次迭代、tip 误差 1.86e-6。
- **程序/数据/图：** [`large_rotation_design.py`](../experiments/cmame_round5/large_rotation_design.py)、[`large_rotation_design.json`](../experiments/cmame_round5/raw_data/large_rotation_design.json)、[`large_rotation_design_history.csv`](../experiments/cmame_round5/raw_data/large_rotation_design_history.csv)、[`Figure_R5_large_rotation_design.pdf`](../experiments/cmame_round5/figures/Figure_R5_large_rotation_design.pdf)。
- **判断：** 三个基础逆问题中，它**最直接**验证几何参数化、参考曲率、有限转动、非线性平衡和目标拟合的耦合。目标仍由已知参数生成，只有两个变量；不能推断高维形状设计一般可辨识或优于成熟壳方法。`DESIGN_REPORT.md` 顶部“demo 4 running”是旧状态，数值文件和报告后文已存在。

## 工程逆向及其前驱

| ID | 阶段 | 设计任务与结果 | 最终可信范围 |
|---|---|---|---|
| I04 | Round 6/7 飞机前驱 | 9 个蒙皮节点力；Round 7 NACA 2412 后缘、0.6×0.3 m、2 mm Al、5 kPa、外给 +5 mm 目标；CB 0.19 s、3.63% 相对目标误差，MITC4 8.98 s、3.55% | 能算简化力反演；节点力不是实际作动器，旧时间混合不同梯度/重建策略 |
| I05 | 工程飞机 | 5 候选安装区域、31 个线性布局、两工况；三处 [0,2,4] 非线性 CB RMS/max 为 0.146/0.412、0.168/0.599 mm；加密 ASD 为 0.181/0.412、0.206/0.644 mm | 经验证的**可行三处方案**；单/双点的局部审计不能证明三处全局最少 |
| I06 | Round 6 汽车前驱 | 1×0.6 简化曲板、9 集中力、自生目标；CB 1.00 s、误差 1.28e-4；MITC4 7.91 s、2.55e-4 | 一般曲板力反演；不是安装点调整或制造偏差装配 |
| I07 | 工程汽车 | 0.8×0.45 m 代表板、0.7 mm 钢、4 有限支架、8 行程；最终 CB 一致梯度 35 次迭代、50 平衡、200.3 s；共同 ASD 24×20 参考最大边缘误差 0.2569 mm、最大支架反力 88.55 N | 同一确定边界下满足研究容限；不能解释为真实车型已标定 |

**I04 溯源：** [`Round 6`](../experiments/cmame_round6_applications/APPLICATIONS_REPORT.md)、[`Round 7`](../experiments/cmame_round7_aircraft/AIRCRAFT_ACTUAL_REPORT.md)。I04 是 I05 的早期不同物理模型，不能把旧九力计时与新布局结果连成一条同一算法加速曲线。

**I05 溯源：** [`aircraft_skin.py`](../experiments/engineering_inverse/aircraft_skin.py)、[`aircraft_nonlinear_reference.py`](../experiments/engineering_inverse/aircraft_nonlinear_reference.py)、[`aircraft_layout_summary.csv`](../experiments/engineering_inverse/results/aircraft_layout_summary.csv)、[`aircraft_corotational_reference.json`](../experiments/engineering_inverse/results/aircraft_corotational_reference.json)。共享算子 31 布局三网格最终 CB 1.688/2.310/3.890 s，MITC4 0.497/0.873/2.033 s；**仅线性筛选**，非完整非线性布局优化。位置式加强线与传统转动梁仍有运动学差异。

**I06 溯源：** [`automotive_panel.py`](../experiments/cmame_round6_applications/automotive_panel.py)、[`automotive_panel.json`](../experiments/cmame_round6_applications/raw_data/automotive_panel.json)。两法按各自生成目标运行；不可用误差差值排名壳单元。

**I07 溯源与修订：** [`automotive_spring_design.py`](../experiments/engineering_inverse/automotive_spring_design.py)、[`tangent_reliability_report.md`](../experiments/engineering_inverse/tangent_reliability_report.md)、[`final optimization CSV`](../experiments/engineering_inverse/results/tangent_automotive_optimization_final.csv)。旧“硬安装区”在同一输入下加密 MITC4 反力 1245 N，旧 CB 缩放设计 264 N，跨模型失败；随后采用有限支架。最终修正第二子格 x 位置基准，旧貌似“多平衡分支”被证明为零能 gauge，统一了 CB/ASD 的物理定位。最终 CB-FD 2509.4 s、766 平衡，CB 一致梯度 200.3 s、50 平衡、相近目标；同模型完整优化收益 12.53×。ASD-FD 3676.5 s、708 平衡，得到**更低的**加密参考目标 0.022173（CB 0.037786），ASD 终点驻点未按 CB 精度独立认证。CB 最终投影 KKT 梯度 4.73e-7 目标单位/mm；这个证据只支持确定冷启动路径的局部驻点。两种壳的完整时间差混合梯度策略、非线性路径、软件栈，不能当作单元速度比。

## 逆向方法与证据矩阵

| 能力或效率问题 | 最直接案例 | 已有证据 | 不能声称 |
|---|---|---|---|
| 力设计梯度 | I01 | 伴随/中心差分核对及同模型优化解数 | 力唯一识别、壳单元固有加速 |
| 参考几何参数化 | I02 | 参考形状改变平衡并拟合位移 | 12 参数真实几何恢复 |
| 大转动+非线性逆向 | I03 | 约 90° 条带的两参数合成反演 | 一般高维可靠性 |
| 多布局重复线性计算 | I05/P03 | 两法共享算子后 MITC4 更快 | CB 完整布局搜索加速 |
| 非线性约束优化与反力梯度 | I07/P05 | 共同参考可行、CB KKT、CB 一致 vs FD | 相同最优目标下跨壳速度优势 |

## 独立研究支线与未发现项

- 外部 [`4d_printing_inverse_design`](</home/song/projects/inverse_design/4d_printing_inverse_design/README.md>) 研究双层打印自然形式/路径，内部 vendored CB；它有自己的材料可达性和非唯一性问题。已登记 X04，**不计本文的基础逆向案例**；未逐文件审计其输出。
- 旧 [`CB-shell/PROJECT_HISTORY.md`](</home/song/projects/CB-shell/PROJECT_HISTORY.md>) 的 GPU、FPC 和候选膜形式是独立支线。FPC 洛伦兹力/接触并非已通过门槛的 CB-Full 逆向实验。
- 在本仓库的 123 个实验 Python 文件、Round 5/6/7 和工程结果目录中，未找到已完成的自然曲率/打印路径反演、第三个非工程大变形逆问题或其他可独立计数的形状匹配案例。若在未入库的私人运行日志中存在，应另补来源；本清查不以关键词臆造案例。
