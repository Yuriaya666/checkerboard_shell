# 全项目独立数值实验总清单

本清单以**研究问题/物理案例**而非网格、厚度、优化方法或轮次计数。每项下的系列见 [CASE_DETAILS.md](CASE_DETAILS.md)，逐系列机器索引见 [EXPERIMENT_MASTER.csv](EXPERIMENT_MASTER.csv)。

检索范围：当前仓库 `experiments/`、`examples/`、`outputs/`、`tests/`、`docs/`、`manuscript/`；原始 arXiv 历史仓库 `/home/song/projects/CB-shell` 的 `v1.0-arxiv` 标签及输出；外部 4D 打印仓库仅作独立支线登记。当前仓库 Git 仅有 `0699638` 基线提交，绝大多数实验目录未被该 Git 追踪，不能从当前 Git 精确恢复每轮运行版本。

**计数：32 项**＝机制 9＋正向 7＋逆向 7＋性能 5＋历史/外部 4。主仓库或原始论文相关内容为 28 项（含历史失败 F07）；X01–X04 是不同程度的独立支线，**不计作本篇壳单元论文的四个新验证**。其中非工程逆向为 I01–I03 共 3 项，工程/工程前驱为 I04–I07 共 4 项。

## 几何与理论机制

| ID | 案例 | 来源 | 系列 | 可信状态 | 证据作用 |
|---|---|---|---|---|---|
| [M01](CASE_DETAILS.md#m01) | B/W 中点基本形式与曲率收敛 | 原始论文；后续复核 | 四曲面 b_W 收敛；Dellinger 对照；后期同 Pq 曲面观测 | 已验证，几何范围有限 | A,B |
| [M02](CASE_DETAILS.md#m02) | 物理商空间、gauge 与刚体客观性 | 原始论文；Round 1/2；定位实验 | P 秩/核；参考和有限构形能量、残差、切线不变性；n=8/12 Pq 诊断 | 多轮复核 | A |
| [M03](CASE_DETAILS.md#m03) | 膜盲模态与 X_M 补全 | 原始论文；Round 3；相对价值研究 | B-only/B+XM 核空间；目标方向刚度；M3 载荷干扰 | 已验证于所测平板/曲面 | A |
| [M04](CASE_DETAILS.md#m04) | W 曲率观测缺失、X_W^rel 与曲面 transport | 原始论文；后续消融 | 曲率秩；有限构形 gauge/刚体；CB-0/M/W/Full 消融 | 有条件验证 | A,D |
| [M05](CASE_DETAILS.md#m05) | 近等距几何、厚度—网格与残余膜锁定 | 原始论文；相对价值研究 | 柱面第一形式/能量；M3/屋顶跨方法厚度扫描 | 部分验证 | A,C,D |
| [M06](CASE_DETAILS.md#m06) | 补全系数敏感性与四模型消融 | Round 1 新增；Round 3 修正；相对研究 | Round 1 reconstructed drivers；Round 3 recovered drivers；Round relative ablation | 已修正；结论有限 | A,D |
| [M07](CASE_DETAILS.md#m07) | 曲面及自由边低端谱 | Round 1 后新增 | 强双层夹持曲面谱；2×2 角块自由边商空间谱 | 有条件验证 | A,D |
| [M08](CASE_DETAILS.md#m08) | 网格畸变与结构适用范围 | Round 1 新增；Round 2 重做 | Round 1 reconstructed/inverted/out-of-scope spacing；Round 2 recovered surface-preserving CB/MT | Round 2 可信；Round 1 部分失效 | B,D |
| [M09](CASE_DETAILS.md#m09) | 相同 Pq 设计几何的传统壳替代检验 | 定位研究新增 | P 秩/gauge；B/W 外部观测；CB/MITC4/ASD 位移与参数灵敏度 | 诊断可信；无认证连续体参考 | A,D |

## 正向壳力学

| ID | 案例 | 来源 | 系列 | 可信状态 | 证据作用 |
|---|---|---|---|---|---|
| [F01](CASE_DETAILS.md#f01) | Simply Supported Plate | 原始论文；Round 2/4 重算 | 原稿 scalar transverse；Round 4 full 3V；跨方法/状态；厚度控制 | 最新 full 3V 可信 | B,D |
| [F02](CASE_DETAILS.md#f02) | Scordelis–Lo Roof | 原始论文；Round 2/4 重算 | 原稿半长 recovered；Round 1 全长误设；双参考值；厚度/网格；三方法 | Recovered driver 可信 | B,D |
| [F03](CASE_DETAILS.md#f03) | MacNeal–Harder M3 Hemisphere | 原始论文；后续多轮 | 标准厚度；raw vs midpoint 观测；载荷修正；厚度×网格；畸变 | 选择性优势；参考与边界有限制 | B,C,D |
| [F04](CASE_DETAILS.md#f04) | 大转动端矩悬臂 | 原始论文；Round 4 重跑 | 16×4/24×4/32×4 路径；最终状态 tangent 复核 | 已复现全路径 | B,E |
| [F05](CASE_DETAILS.md#f05) | 非线性交替力半球 | 原始论文；Round 4 重跑 | n=6/8/12/16 全路径；末态 sparse tangent | 已复现；数字化参考有限精度 | B |
| [F06](CASE_DETAILS.md#f06) | 开端圆柱壳 pull-out | 原始论文；Round 4 重跑 | 8×8/12×12/16×16 全路径；C 点位移反转；末态 tangent | 已复现；局部误差保留 | B,D |
| [F07](CASE_DETAILS.md#f07) | 闭端 pinched cylinder 折叠分支探索 | 原始历史仓库；原稿未纳入 | 载荷控制 6/8/12 网格；SVD null-vector pseudo-arclength；分支诊断 | 失败/未纳入论文 | D |

## 逆向设计

| ID | 案例 | 来源 | 系列 | 可信状态 | 证据作用 |
|---|---|---|---|---|---|
| [I01](CASE_DETAILS.md#i01) | 九节点作动力平板形状控制 | Round 5 新增；工程轮次回归 | adjoint 中心差分；CB adjoint/FD；MITC4 FD own/same target；后期重跑 | 同模型梯度强验证；力识别非唯一 | E |
| [I02](CASE_DETAILS.md#i02) | 12 参数参考形状优化 | Round 5 新增 | 参考形状参数与目标位移拟合；优化历史 | 目标拟合有效；参数恢复失败 | E,D |
| [I03](CASE_DETAILS.md#i03) | 90°级大转动参考几何反演 | Round 5 新增 | 两参数曲率恢复；非线性平衡路径；目标 tip | 已完成；小参数合成逆问题 | E |
| [I04](CASE_DETAILS.md#i04) | 早期飞机蒙皮九力反演 | Round 6/7 后续新增 | R6 各自生成目标；R7 外给 +5mm 目标；adjoint/黑箱 FD | 历史/已被工程布局替代 | F,D |
| [I05](CASE_DETAILS.md#i05) | 飞机五候选作动布局与行程 | engineering_inverse 后续新增 | 31 线性布局；单/双点非线性本地审计；[0,2,4] 非线性；参考与复用性能 | 有可行布局；最少数量未证明 | F,D |
| [I06](CASE_DETAILS.md#i06) | 早期汽车曲板九力 gap/flush 反演 | Round 6 后续新增 | CB adjoint 与 MITC4 FD 各自生成目标 | 历史/被安装调整替代 | F,D |
| [I07](CASE_DETAILS.md#i07) | 汽车覆盖件非线性安装调整 | engineering_inverse 后续新增并多次修正 | 硬区域旧模型；有限支架；共同 locator/gauge；CB consistent/FD 与 ASD-FD；KKT | 研究容限通过；物理支架未标定 | E,F,D |

## 性能与实现

| ID | 案例 | 来源 | 系列 | 可信状态 | 证据作用 |
|---|---|---|---|---|---|
| [P01](CASE_DETAILS.md#p01) | 精确稀疏＋秩一切线与边界存储 | 理论原稿；Round 1/2 生产化 | Round 1 dense/reference-linear 原型；Round 2 exact nonlinear sparse border；内存/nnz/solve | 精确结构已验证 | A,D |
| [P02](CASE_DETAILS.md#p02) | 跨方法精度—状态数—时间对照 | 原稿 MT；Round 3/4 MITC；相对研究 | recovered CB/MT；MITC4；同目标精度阈值；CPU/RSS | 混合结果；限定标量/实现 | C,D |
| [P03](CASE_DETAILS.md#p03) | 飞机 31 布局共享响应算子对照 | engineering_inverse 后续新增 | CB direct/shared；MITC direct/shared；三次中位；同解一致性 | 公平共享筛选；仅线性 | D |
| [P04](CASE_DETAILS.md#p04) | 局部 Hessian 瓶颈与缓存实验 | engineering_inverse 最终审计 | 147 graph colours；legacy/cached AD batch；力学同解回归 | 实现等价；收益混合 | D |
| [P05](CASE_DETAILS.md#p05) | 汽车三种完整非线性优化成本与梯度归因 | engineering_inverse 最终审计 | CB consistent；CB 8变量中心 FD；ASD corotational FD；KKT/branch | 同 CB 梯度收益强；跨壳限于达容限 | E,D |

## 历史与独立支线

| ID | 案例 | 来源 | 系列 | 可信状态 | 证据作用 |
|---|---|---|---|---|---|
| [X01](CASE_DETAILS.md#x01) | 旧仓库替代膜/补全形式探索 | 外部历史仓库；未迁入正式内核 | full-energy 反例；locking-free extension；最终信息层级搜索 | 历史失败；与当前正式能量隔离 | D |
| [X02](CASE_DETAILS.md#x02) | 旧仓库 matrix-free/GPU 批量计算 | 外部历史仓库；非本篇当前 CPU 结论 | matrix-free CPU；GPU HVP；多场景与失败区间 | 独立历史结果 | D |
| [X03](CASE_DETAILS.md#x03) | 旧仓库 FPC 洛伦兹/接触柔性片 | 外部历史项目；非本论文实验 | 节点/积分点接触；后屈曲分支；三角壳替代；细网格门槛 | 独立且含失败结果 | D |
| [X04](CASE_DETAILS.md#x04) | 4D 打印双层自然形式逆向 | 外部独立仓库；不计本文逆向案例 | 几何可达性、打印参数非唯一、目标分支和前向验收 | 独立项目，未逐项审计 | E,D |
