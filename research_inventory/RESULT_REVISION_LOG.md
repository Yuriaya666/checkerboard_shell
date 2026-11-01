# 结果修订、重复与争议台账

此文件保留**旧结果→发现的问题→后续记录→当前引用口径**。未经重新计算。文中的“失效”指不能继续支持原主张，原始数据仍保留。

## 版本来源边界

- 当前 Git 仓库只见基线提交 `0699638944e35262dec343a13e6398da203b91f7`（2026-09-04）；`experiments/` 多数文件未被当前 Git 追踪。因此当前仓库无法独立提供 Round 1–7 的逐提交精确来源。**文件名/报告时间不等同 Git 运行版本。**
- 邻接旧仓库 `/home/song/projects/CB-shell` 有 `v1.0-arxiv`，指向 `093ad89aa0f56b513dd89e207cfb1d69ca4a21f0`（2026-08-29）。该标签中的 `paper_draft_arxiv_cn/figures/` 是原稿图的可信路径；旧仓库工作树中的大量 `outputs/` 是后续历史运行，不能仅因目录名为 `paper_benchmark` 就视为标签时发布的原始数据。
- `research_inventory/catalog.py` 对每项标“原始论文/后续 Round/外部支线”。`COVERAGE_AUDIT.json` 检查代表路径实际存在；逐文件图表是否已在原稿出现，以标签和 `FIGURE_AND_DATA_INDEX.csv` 路径核查。

| 编号/案例 | 历史数据或主张 | 发现与修正 | 当前可信引用口径 | 证据路径 |
|---|---|---|---|---|
| R01 F01–F03 | Round 1 重建板/屋顶/M3 驱动得出位移、畸变和 beta 结论 | Round 2 找到原稿真实边界、载荷、观测量；旧驱动的 BVP 不一致 | **禁用** Round 1 重建基准位移排名；采用 Round 2 recovered drivers、Round 4 同观测对照 | [`Round 2`](../experiments/cmame_round2/CMAME_ROUND2_REPORT.md)、[`run_driver_recovery.py`](../experiments/cmame_round2/run_driver_recovery.py) |
| R02 F01 | 早期板中心观测使用不同节点/物理中面约定 | Round 4 统一 full 3V 和中心观测；n=8/16/32/64 CB 误差为 5.5475/1.4479/0.3605/0.0664% | 写明是 full 3V 还是旧 scalar active；不用不同观测行拼接收敛率 | [`Round 4`](../experiments/cmame_round4/CMAME_ROUND4_REPORT.md) |
| R03 F02 | 屋顶历史参考 `0.3024` 给出 n=32 CB 2.541% | 现代参考 `0.3005924566` 给出同响应 `0.310083` 的误差 3.157% | 两参考都列或明确现代参考；MidedgeTan n=32 对现代参考误差 0.260% | [`Round 4`](../experiments/cmame_round4/CMAME_ROUND4_REPORT.md) |
| R04 F03 | M3 n=12 “CB 2.20%，MT 3.74%”用于精度/状态结论 | CB 2.20% 对应重建的物理中点观测；raw 测点 CB `0.0909082` 对标准 `0.09352155` 为 2.794%；不可跨观测排位 | 对同一原始观测比较：n=12 CB 2.794% 与 MT 3.741%；如展示 2.201% 必标为 midpoint-reconstructed | [`Round 4`](../experiments/cmame_round4/CMAME_ROUND4_REPORT.md)、[`relative value`](../experiments/cmame_relative_value/report.md) |
| R05 M05/F03 | 固定网格或仅标准厚度的 M3 较好响应被解读为普遍无膜锁定 | 后续 4 网格×5 厚度显示固定 h 仍有残余锁定；极薄 M3 n=16 CB/MT/MITC 响应差约 19.73/17.28/88.92% | CB 在此 MITC 实现前响应抑制较弱，但不优于 MidedgeTan；不宣称 locking-free | [`relative value`](../experiments/cmame_relative_value/report.md) |
| R06 M08 | Round 1 随机畸变 δ=0.40 的三个失败曾混入方法稳健性解释 | 原报告注明面片倒置/无效；Round 2 用 δ=0.30 surface-preserving 重新比较；M3 CB 平滑/随机误差 20.124/28.895% | 保留失效网格，引用 Round 2 有效物理曲面对照；CB 并无一般畸变稳健性优势 | [`additional`](../experiments/cmame_additional/CMAME_NUMERICAL_REPORT.md)、[`Round 2`](../experiments/cmame_round2/CMAME_ROUND2_REPORT.md) |
| R07 M08 | 非均匀逻辑间距扫描被解释为当前 CB 方法适用结果 | 现有几何与差分假设使用固定逻辑步长，非均匀间距不在已验证方法范围 | 只作为失效/适用范围提醒，不用于正式精度排名 | [`Round 2`](../experiments/cmame_round2/CMAME_ROUND2_REPORT.md) |
| R08 M06 | Round 1 beta 对 roof/M3 正常响应无影响 | 当时使用 reconstructed drivers；Round 3 用 recovered drivers 重跑：roof n=16 c=.25–4 最大位移差 1.63e-5，M3 n=12 为 4.21e-10 | 使用 Round 3 数据；只说明所测 BVP 的低敏感性 | [`Round 3`](../experiments/cmame_round3/CMAME_ROUND3_REPORT.md) |
| R09 P01 | Round 1 dense 边界矩阵/参考线性 prototype 内存估计 | Round 2 实现局部边界 stencil、精确 nonlinear sparse+rank-one 生产切线；M3 n=128 11.49s assembly/16.39s factor/2073MiB | 用实测生产内存/时长；旧 dense 投影仅作为修改历史 | [`Round 2`](../experiments/cmame_round2/CMAME_ROUND2_REPORT.md) |
| R10 F04–F06 | 原文非线性只展示终点或个别载荷 | Round 4 整理完整路径、加密及 tangent：悬臂 n=32×4 路径 RMS 2.442%；半球 n=16 联合 RMS 1.736%；开端圆柱 n=16×16 联合 RMS 4.483% | 路径和终点一起报告；外部半球数字化参考有有限读数不确定性 | [`Round 4`](../experiments/cmame_round4/CMAME_ROUND4_REPORT.md) |
| R11 F07 | 旧仓库闭端 pinched cylinder 重载荷/折叠分支算例 | 载荷控制和伪弧长尝试均未可靠重接 post-snap 分支；原观察值不能当作已验证解 | 保留失败与负切线/路径诊断；不纳入“已通过非线性基准” | [`历史审计`](</home/song/projects/CB-shell/docs/checkerboard_nonlinear_benchmark_audit.md>) |
| R12 I01 | 9 力优化被描述为“真实作动力识别” | 力恢复相对误差 0.174，但目标形状 0.0065%；逆映射非唯一 | 称“形状控制/合成逆设计”；不要称唯一力识别 | [`DESIGN_REPORT.md`](../experiments/cmame_round5/DESIGN_REPORT.md) |
| R13 I01 | CB-adjoint 1.05s 与 MITC4-FD 11.85s 被归因为 3V 优势 | 比较混合单元、梯度、模型建立、软件栈；历史 MITC own mesh 12×12 与 CB 12×8 不匹配 | 同 CB adjoint/FD 79/438 线性解证明梯度工作量；跨壳时间只作工作流参考，MITC own-target 误差还更小 | [`FD_VS_FD.md`](../experiments/cmame_round5/FD_VS_FD.md)、[`DESIGN_REPORT.md`](../experiments/cmame_round5/DESIGN_REPORT.md) |
| R14 I02/I03 | `DESIGN_REPORT.md` 顶部“demo 4 running” | 后文与 JSON/CSV 已有 I02/I03 完成结果 | 以原始 JSON/history 和后文数值为准，旧状态行仅为文档未更新 | [`Round 5 raw_data`](../experiments/cmame_round5/raw_data) |
| R15 I04/I06 | Round 6/7 九力曲板和蒙皮作为“工程作动/装配” | 物理装置、安装路径和部分目标来自自身求解；后来建立 I05/I07 的独立工程问题 | 旧模型留作历史过渡，不能当真正作动器数目/汽车安装工程约束结果 | [`Round 6`](../experiments/cmame_round6_applications/APPLICATIONS_REPORT.md)、[`Round 7`](../experiments/cmame_round7_aircraft/AIRCRAFT_ACTUAL_REPORT.md) |
| R16 I05 | 线性筛选认为三处 [0,2,4] 最少，个别 1/2 点 Newton 失败被判不可行 | 后续续接/本地审计为 5 单点+10 双点找到收敛状态，但没有完整非线性全局可行性排除；线性最优布局随网格变 | [0,2,4] 是两工况经独立加密核验的**可行**布局，三处最少未证明 | [`engineering report`](../experiments/engineering_inverse/report.md)、[`nonlinear audit data`](../experiments/engineering_inverse/results/aircraft_nonlinear_audit_12x24.json) |
| R17 I07 | 刚性安装区粗 CB 缩放行程后 264N、0.385mm | 相同行程在加密 MITC4 得 1245N，区域连接/反力网格依赖严重；后改有限刚度支架并做共同参考 | 旧结果失效为工程可行设计；引用有限支架 + 新边界的最终结果 | [`RESEARCH_REPORT.md`](../experiments/engineering_inverse/RESEARCH_REPORT.md)、[`report.md`](../experiments/engineering_inverse/report.md) |
| R18 I07 | 旧 CB 优化 328s/42 次，后 238s、1384s、3560s 三法表 | 连接、x locator、路径及缓存几轮变化；旧不同物理定位/优化终点不能同最终表比较 | 最终共同 locator：CB consistent 200.3s/50 平衡、CB-FD 2509.4s/766、ASD-FD 3676.5s/708 | [`tangent reliability`](../experiments/engineering_inverse/tangent_reliability_report.md) |
| R19 I07 | 同行程冷/热初值解似有不同稳定分支 | 差异为相反棋盘子格 x 向 0.956mm 的零能 gauge；加第二子格 x datum；ASD 共同物理约束也重跑 | 不再把此现象声称物理多稳态；CB 冷路径 KKT 4.73e-7；ASD 端点驻点精度仍未完全认证 | [`tangent reliability`](../experiments/engineering_inverse/tangent_reliability_report.md) |
| R20 P03/P04 | 旧飞机 31 布局 CB 1.24/2.15/4.67s 对 MITC 9.46/23.15/37.26s，约 8×“单元加速” | 旧 MITC 每布局重建；后实现同等共享算子，三网格最终 CB 1.688/2.310/3.890s，MITC 0.497/0.873/2.033s | MITC 对此线性筛选**更快**；CB 20×40 局部 Hessian 约 3.13s 主导 | [`efficiency_report.md`](../experiments/engineering_inverse/efficiency_report.md)、[`tangent reliability`](../experiments/engineering_inverse/tangent_reliability_report.md) |
| R21 P05 | 汽车一致梯度 vs ASD-FD 约 18×被表述为 CB 单元优势 | 最终 CB 同模型一致梯度/FD 为 200.3/2509.4s，完整优化加速 12.53×；ASD-FD 的共同参考目标 0.022173 比 CB 0.037786 更低，且不同软件/路径 | 12.53× 是可微工作流收益，不是 CB 几何独有；跨壳没有等最优目标加速比 | [`tangent reliability`](../experiments/engineering_inverse/tangent_reliability_report.md) |
| R22 M09 | “原生 Pq 几何必然免除传统壳转换误差”的候选卖点 | 同一 (Pq) 给 MITC4/ASD 建模可复用外部参数，几何 B/W 观测约 1e-15，设计响应可计算；CB 并无独占转换能力 | Pq 的兼容几何优点真实，但传统壳可共享；不称独有应用必要性。MITC 力矩账例外仍未闭合 | [`positioning report`](../experiments/cmame_positioning/report.md) |
| R23 X01–X03 | 旧候选膜、GPU、FPC 输出容易误并入当前 CB-Full | 旧膜报告最终 `NO_UNIFIED_CHECKERBOARD_MEMBRANE_FORMULATION_FOUND`；FPC Mechanics Gate 未通过；GPU 是不同硬件/批量 regime | 分别以 X 编号保留，不能混入当前能量的验证结果或当前 CPU 速度结论 | [`PROJECT_HISTORY.md`](</home/song/projects/CB-shell/PROJECT_HISTORY.md>) |

## 当前仍未解决、不可自行补写的缺口

1. 当前 `experiments/` 缺 Git 逐轮提交 provenance；多数原始 JSON/CSV 有文件路径但不能以当前 Git 证明它们在原稿标签时存在。
2. 图表多处只有目录级原始数据候选，逐图生成脚本/具体 CSV 需人工复核；索引明确标记，不臆造。
3. F07 的闭端折叠分支没有可信完整路径；I05 三处“最少数量”没有全局证明。
4. I07 ASDShellQ4 最优终点驻点精度及 OpenSees 反力矩账仍有保留；有限支架刚度是代表性设定，非实车标定。
5. 缺同条件、可复现 BSQ 等 3V rotation-free 基线；不能由现有 MidedgeTan/MITC4 对照推断相对所有位置式壳单元更优。
