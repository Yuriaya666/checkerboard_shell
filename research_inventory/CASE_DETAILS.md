# 独立案例详细记录

数值均引用已有记录；未重新运行求解器。`源程序`、`原始数据`、`代表图表`是已核实存在的**代表路径**，全部文件见 `ARTIFACT_MANIFEST.csv` 和 `FIGURE_AND_DATA_INDEX.md`。一项案例可以横跨多轮，系列名不另计案例。

<a id="m01"></a>

## M01 B/W 中点基本形式与曲率收敛

- **原始研究目的：** 检验 B 第一形式、W 第二形式的几何一致性。
- **物理对象与设定：** 柱面、球面、扭曲和 skew 曲面；光滑网格加密。
- **实验系列：** 四曲面 b_W 收敛；Dellinger 对照；后期同 Pq 曲面观测
- **方法：** CB B/W；Dellinger 几何构造
- **主要数值结果：** 原稿四曲面二阶左右收敛；新 Pq 诊断的 B/W 外部观测差约 1e-15。
- **参考与比较基准：** 连续基本形式；Dellinger 几何量
- **版本/轮次：** 原始论文；后续复核
- **验证状态：** 已验证，几何范围有限
- **已知问题：** 该几何谱系已有先行研究；不能计为本文独创，也不能仅由几何误差推断力学精度。
- **能够支持：** 统一几何量的平滑收敛。
- **不能支持：** 普遍更优壳力学或独占几何转换优势。
- **论文证据作用：** A,B
- **源程序：** `/home/song/projects/CB-shell/experiments/checkerboard_paper_benchmark/run_geometry.py`；[`experiments/cmame_positioning/diagnostic.py`](../experiments/cmame_positioning/diagnostic.py)
- **原始数据：** `/home/song/projects/CB-shell/outputs/checkerboard_paper_benchmark/geometry_convergence.csv`；[`experiments/cmame_positioning/results/diagnostic.json`](../experiments/cmame_positioning/results/diagnostic.json)
- **代表图表：** `/home/song/projects/CB-shell/paper_draft_arxiv_cn/figures/figure_04_geometry_convergence.pdf`

<a id="m02"></a>

## M02 物理商空间、gauge 与刚体客观性

- **原始研究目的：** 区分 raw checkerboard gauge、物理形变及刚体零模态。
- **物理对象与设定：** 连通二分网格；有限 gauge、刚体运动和变形状态。
- **实验系列：** P 秩/核；参考和有限构形能量、残差、切线不变性；n=8/12 Pq 诊断
- **方法：** CB-Full；独立 midpoint map
- **主要数值结果：** ker P 维数 3；Round 1 最大切线相对变化 1.23e-15；新 gauge 能量差 2.08e-17 J。
- **参考与比较基准：** 代数秩；相同物理中面
- **版本/轮次：** 原始论文；Round 1/2；定位实验
- **验证状态：** 多轮复核
- **已知问题：** 载荷自身必须 gauge 不变；raw 节点固定不能冒充物理边界。
- **能够支持：** 商空间与有限构形客观性。
- **不能支持：** 较其他壳单元更准确或更快。
- **论文证据作用：** A
- **源程序：** [`experiments/cmame_additional/run_primary.py`](../experiments/cmame_additional/run_primary.py)；[`experiments/cmame_positioning/gauge_check.py`](../experiments/cmame_positioning/gauge_check.py)
- **原始数据：** [`experiments/cmame_additional/02_gauge_invariance/raw_data/invariance_errors.csv`](../experiments/cmame_additional/02_gauge_invariance/raw_data/invariance_errors.csv)；[`experiments/cmame_positioning/results/gauge_check.json`](../experiments/cmame_positioning/results/gauge_check.json)
- **代表图表：** [`experiments/cmame_additional/02_gauge_invariance/figures/Figure_G1.pdf`](../experiments/cmame_additional/02_gauge_invariance/figures/Figure_G1.pdf)

<a id="m03"></a>

## M03 膜盲模态与 X_M 补全

- **原始研究目的：** 验证 B 度量的物理缺失方向和秩一膜兼容补全。
- **物理对象与设定：** 平板 in-plane 线性算子；M3 正常加载响应。
- **实验系列：** B-only/B+XM 核空间；目标方向刚度；M3 载荷干扰
- **方法：** CB-0；CB-M；CB-Full
- **主要数值结果：** 平板核维 6→5；目标方向刚度约 8.27e-5；M3 正常响应 B→B+XM 变化至多 1.19e-10。
- **参考与比较基准：** 解析秩计数；正谱与载荷响应
- **版本/轮次：** 原始论文；Round 3；相对价值研究
- **验证状态：** 已验证于所测平板/曲面
- **已知问题：** 曲面受约束 full tangent 去掉 XM 未必新增精确零模态；勿把 gauge 算作材料盲模。
- **能够支持：** XM 对目标缺失方向的针对性及低正常响应干扰。
- **不能支持：** 普遍无锁定。
- **论文证据作用：** A
- **源程序：** `/home/song/projects/CB-shell/experiments/checkerboard_paper_benchmark/run_mechanism.py`；[`experiments/cmame_round3/run_beta_revalidated.py`](../experiments/cmame_round3/run_beta_revalidated.py)
- **原始数据：** `/home/song/projects/CB-shell/outputs/checkerboard_paper_benchmark/mechanism_results.json`；[`experiments/cmame_round3/01_beta_revalidated/raw_data/beta_revalidated.csv`](../experiments/cmame_round3/01_beta_revalidated/raw_data/beta_revalidated.csv)
- **代表图表：** `/home/song/projects/CB-shell/paper_draft_arxiv_cn/figures/figure_05_minimal_completion.pdf`

<a id="m04"></a>

## M04 W 曲率观测缺失、X_W^rel 与曲面 transport

- **原始研究目的：** 验证 W-centered 曲率盲通道、参考相对补全和客观 transport。
- **物理对象与设定：** 平/曲参考、有限转动；四补全模型。
- **实验系列：** 曲率秩；有限构形 gauge/刚体；CB-0/M/W/Full 消融
- **方法：** CB variants；原稿曲面 transport
- **主要数值结果：** 原稿 fine-grid XW 引起板/屋顶/M3 位移变化 0.0136%/0.0210%/0.2983%；后续 roof/M3 正常模态轻微变化。
- **参考与比较基准：** AD 导数、秩和加载响应
- **版本/轮次：** 原始论文；后续消融
- **验证状态：** 有条件验证
- **已知问题：** 初始曲面和边界会改变谱；曲面 full-tangent 数值近零数不是一般核定理。
- **能够支持：** 曲率兼容补全和有限参考客观性。
- **不能支持：** 对所有物理弯曲模态无干扰。
- **论文证据作用：** A,D
- **源程序：** [`experiments/cmame_relative_value/run_completion.py`](../experiments/cmame_relative_value/run_completion.py)；[`experiments/cmame_additional/run_spectra.py`](../experiments/cmame_additional/run_spectra.py)
- **原始数据：** [`experiments/cmame_relative_value/results/completion_ablation.csv`](../experiments/cmame_relative_value/results/completion_ablation.csv)；[`experiments/cmame_relative_value/results/completion_prior_spectrum_summary.csv`](../experiments/cmame_relative_value/results/completion_prior_spectrum_summary.csv)
- **代表图表：** [`experiments/cmame_relative_value/figures/completion_ablation.pdf`](../experiments/cmame_relative_value/figures/completion_ablation.pdf)

<a id="m05"></a>

## M05 近等距几何、厚度—网格与残余膜锁定

- **原始研究目的：** 检验 t/h 耦合、柱面近等距和薄壳残余刚化。
- **物理对象与设定：** 平面至柱面规定等距；M3/屋顶 4 网格×5 厚度。
- **实验系列：** 柱面第一形式/能量；M3/屋顶跨方法厚度扫描
- **方法：** CB-Full；MidedgeTan；MITC4
- **主要数值结果：** 柱面 n=8→32 第一形式误差 4.091e-4→2.557e-5；M3 n=16 最薄 CB/MT/MITC4 响应差 19.73%/17.28%/88.92%。
- **参考与比较基准：** 解析等距；逐厚度加密外部参考
- **版本/轮次：** 原始论文；相对价值研究
- **验证状态：** 部分验证
- **已知问题：** CB 固定网格仍有残余锁定；最薄 M3 不优于 MT；屋顶并非纯弯曲。
- **能够支持：** 受限族近等距逼近与具体 M3 响应趋势。
- **不能支持：** locking-free 或一般薄壳优势。
- **论文证据作用：** A,C,D
- **源程序：** [`experiments/cmame_relative_value/run_scan.py`](../experiments/cmame_relative_value/run_scan.py)；[`experiments/cmame_relative_value/extend_thin_reference.py`](../experiments/cmame_relative_value/extend_thin_reference.py)
- **原始数据：** [`experiments/cmame_relative_value/results/cylinder_isometry.csv`](../experiments/cmame_relative_value/results/cylinder_isometry.csv)；[`experiments/cmame_relative_value/results/thickness_scan.csv`](../experiments/cmame_relative_value/results/thickness_scan.csv)
- **代表图表：** [`experiments/cmame_relative_value/figures/thickness_mesh.pdf`](../experiments/cmame_relative_value/figures/thickness_mesh.pdf)

<a id="m06"></a>

## M06 补全系数敏感性与四模型消融

- **原始研究目的：** 辨别补全是否是单基准调参、是否明显干扰正常响应。
- **物理对象与设定：** 理论 beta_M=8mu；c_beta=0,.25,.5,1,2,4,8；CB-0/M/W/Full。
- **实验系列：** Round 1 reconstructed drivers；Round 3 recovered drivers；Round relative ablation
- **方法：** CB variants
- **主要数值结果：** Round 3 roof n=16 c=.25–4 响应最大变化 1.63e-5；M3 n=12 为 4.21e-10；roof/M3 CB-0→Full 低端特征值约 +0.26%/+1.03%。
- **参考与比较基准：** 固定理论系数；同模型扰动
- **版本/轮次：** Round 1 新增；Round 3 修正；相对研究
- **验证状态：** 已修正；结论有限
- **已知问题：** Round 1 reconstructed benchmark 响应被 Round 3 取代；不能扩展到非线性所有载荷。
- **能够支持：** 所测正常响应不依赖狭窄系数窗口。
- **不能支持：** 可任意去掉补全或普遍参数无关。
- **论文证据作用：** A,D
- **源程序：** [`experiments/cmame_additional/run_primary.py`](../experiments/cmame_additional/run_primary.py)；[`experiments/cmame_round3/run_beta_revalidated.py`](../experiments/cmame_round3/run_beta_revalidated.py)；[`experiments/cmame_relative_value/run_completion.py`](../experiments/cmame_relative_value/run_completion.py)
- **原始数据：** [`experiments/cmame_additional/04_beta_sensitivity/raw_data/beta_sensitivity.csv`](../experiments/cmame_additional/04_beta_sensitivity/raw_data/beta_sensitivity.csv)；[`experiments/cmame_round3/01_beta_revalidated/raw_data/beta_revalidated.csv`](../experiments/cmame_round3/01_beta_revalidated/raw_data/beta_revalidated.csv)；[`experiments/cmame_relative_value/results/completion_ablation.csv`](../experiments/cmame_relative_value/results/completion_ablation.csv)
- **代表图表：** [`experiments/cmame_round3/01_beta_revalidated/figures/Figure_Beta_revalidated.pdf`](../experiments/cmame_round3/01_beta_revalidated/figures/Figure_Beta_revalidated.pdf)

<a id="m07"></a>

## M07 曲面及自由边低端谱

- **原始研究目的：** 查找曲率改变或自由边引入的额外软模态。
- **物理对象与设定：** 球/鞍面 n=8/16/32；平/柱自由条带。
- **实验系列：** 强双层夹持曲面谱；2×2 角块自由边商空间谱
- **方法：** CB-Full 观测算子与 full tangent
- **主要数值结果：** 所测球/鞍 full tangent 无额外近零模态；自由条带正最小特征值随细化下降。
- **参考与比较基准：** 受约束谱；保存特征向量
- **版本/轮次：** Round 1 后新增
- **验证状态：** 有条件验证
- **已知问题：** 强夹持/角块不是所有自然边界；初次 ARPACK 停滞，后用 dense 子集。
- **能够支持：** 特定边界下未见持久虚假零模。
- **不能支持：** 所有自由边稳定性定理。
- **论文证据作用：** A,D
- **源程序：** [`experiments/cmame_additional/run_spectra.py`](../experiments/cmame_additional/run_spectra.py)
- **原始数据：** [`experiments/cmame_additional/05_curved_spectrum/raw_data/spectra.csv`](../experiments/cmame_additional/05_curved_spectrum/raw_data/spectra.csv)；[`experiments/cmame_additional/06_boundary_modes/raw_data/spectra.csv`](../experiments/cmame_additional/06_boundary_modes/raw_data/spectra.csv)
- **代表图表：** [`experiments/cmame_additional/05_curved_spectrum/figures/Figure_K3.pdf`](../experiments/cmame_additional/05_curved_spectrum/figures/Figure_K3.pdf)

<a id="m08"></a>

## M08 网格畸变与结构适用范围

- **原始研究目的：** 检验非矩形嵌入和表面保持扰动下的响应。
- **物理对象与设定：** 板/屋顶/M3；平滑和随机切向扰动；固定逻辑步长。
- **实验系列：** Round 1 reconstructed/inverted/out-of-scope spacing；Round 2 recovered surface-preserving CB/MT
- **方法：** CB-Full；MidedgeTan
- **主要数值结果：** Round 2 δ=.30 M3 CB smooth/random 20.124%/28.895%，MT 30.219%/43.742%；板/屋顶 MT 更好。
- **参考与比较基准：** 有效网格 Jacobian、共同顶点/载荷
- **版本/轮次：** Round 1 新增；Round 2 重做
- **验证状态：** Round 2 可信；Round 1 部分失效
- **已知问题：** Round 1 非均匀步长超出当前公式范围；随机 δ=.40 三网格翻转；M3 对点重分布敏感。
- **能够支持：** 固定逻辑网格有效扰动可计算但精度混合。
- **不能支持：** 普遍畸变鲁棒性。
- **论文证据作用：** B,D
- **源程序：** [`experiments/cmame_round2/run_distortion.py`](../experiments/cmame_round2/run_distortion.py)；[`experiments/cmame_round2/run_midedgetan_distortion.py`](../experiments/cmame_round2/run_midedgetan_distortion.py)
- **原始数据：** [`experiments/cmame_round2/04_distortion_revalidated/raw_data/distortion.csv`](../experiments/cmame_round2/04_distortion_revalidated/raw_data/distortion.csv)；[`experiments/cmame_additional/01_mesh_distortion/raw_data/all_cases.csv`](../experiments/cmame_additional/01_mesh_distortion/raw_data/all_cases.csv)
- **代表图表：** [`experiments/cmame_round2/04_distortion_revalidated/figures/Figure_D1_distortion_comparison.pdf`](../experiments/cmame_round2/04_distortion_revalidated/figures/Figure_D1_distortion_comparison.pdf)

<a id="m09"></a>

## M09 相同 Pq 设计几何的传统壳替代检验

- **原始研究目的：** 验证传统壳是否必然损失 Checkerboard 参考几何/参数。
- **物理对象与设定：** 合成非均匀曲率 Pq；p=±.1,0；n=8–80 固定输入。
- **实验系列：** P 秩/gauge；B/W 外部观测；CB/MITC4/ASD 位移与参数灵敏度
- **方法：** CB-Full；ShellMITC4；ASDShellQ4
- **主要数值结果：** B/W 观测差约 1e-15；CB n=32 几何参数灵敏度 0.01780 mm/unit，外部 n=64 约 0.01965–0.01967。
- **参考与比较基准：** 共同 Pq 参考几何；两外部细网格
- **版本/轮次：** 定位研究新增
- **验证状态：** 诊断可信；无认证连续体参考
- **已知问题：** 合成对象非发表设计；MITC4 曲面力矩账不闭合；变形运动学不同。
- **能够支持：** 无必然几何转换代价；Pq 兼容维数。
- **不能支持：** CB 独占应用优势或精度排序。
- **论文证据作用：** A,D
- **源程序：** [`experiments/cmame_positioning/diagnostic.py`](../experiments/cmame_positioning/diagnostic.py)；[`experiments/cmame_positioning/run_reference.py`](../experiments/cmame_positioning/run_reference.py)
- **原始数据：** [`experiments/cmame_positioning/results/diagnostic.json`](../experiments/cmame_positioning/results/diagnostic.json)；[`experiments/cmame_positioning/results/summary.csv`](../experiments/cmame_positioning/results/summary.csv)
- **代表图表：** [`experiments/cmame_positioning/figures/design_sensitivity.pdf`](../experiments/cmame_positioning/figures/design_sensitivity.pdf)

<a id="f01"></a>

## F01 Simply Supported Plate

- **原始研究目的：** 基本平板弯曲精度和全 3V 状态验证。
- **物理对象与设定：** 单位方板，E=1, nu=.3, t=.01, 压力 1e-8，四边 w=0。
- **实验系列：** 原稿 scalar transverse；Round 4 full 3V；跨方法/状态；厚度控制
- **方法：** CB-Full；MidedgeTan；MITC4
- **主要数值结果：** full 3V n=8/16/32/64 误差 5.5475/1.4479/.3605/.0664%；n=64 MT .0782%；MITC4 粗网格更准确。
- **参考与比较基准：** Navier 4.4360891e-4
- **版本/轮次：** 原始论文；Round 2/4 重算
- **验证状态：** 最新 full 3V 可信
- **已知问题：** 原稿 scalar branch 225 vs MT 1025 的 active-DOF 口径不公平；Round 1 smoothed center 非原文观测。
- **能够支持：** 平板收敛与 full 3V 基本正确性。
- **不能支持：** 普遍板精度/CPU 优势。
- **论文证据作用：** B,D
- **源程序：** [`experiments/cmame_round4/plate_full.py`](../experiments/cmame_round4/plate_full.py)；[`experiments/cmame_round4/mitc4.py`](../experiments/cmame_round4/mitc4.py)
- **原始数据：** [`experiments/cmame_round4/raw_data/plate_full_crosscheck.csv`](../experiments/cmame_round4/raw_data/plate_full_crosscheck.csv)；[`experiments/cmame_relative_value/results/accuracy_states.csv`](../experiments/cmame_relative_value/results/accuracy_states.csv)
- **代表图表：** [`experiments/cmame_round4/figures/Figure_R4_full_plate.pdf`](../experiments/cmame_round4/figures/Figure_R4_full_plate.pdf)

<a id="f02"></a>

## F02 Scordelis–Lo Roof

- **原始研究目的：** 曲率/膜弯耦合及端隔板边界下的经典响应。
- **物理对象与设定：** R=25,L=50,t=.25,E=4.32e8,nu=0, global -z 压力 90。
- **实验系列：** 原稿半长 recovered；Round 1 全长误设；双参考值；厚度/网格；三方法
- **方法：** CB-Full；MidedgeTan；MITC4
- **主要数值结果：** n=32 CB 位移 .310083；相对现代 .30059246 误差 3.157%；MT .260%；MITC4 同网格更准。
- **参考与比较基准：** 历史 .3024；现代 .3005924566
- **版本/轮次：** 原始论文；Round 2/4 重算
- **验证状态：** Recovered driver 可信
- **已知问题：** Round 1 几何/边界/载荷/观测不同；现代参考比历史参考更不利于 CB。
- **能够支持：** 曲壳基本收敛。
- **不能支持：** 屋顶精度优势或纯弯曲锁定判据。
- **论文证据作用：** B,D
- **源程序：** [`experiments/cmame_round2/recovered_drivers.py`](../experiments/cmame_round2/recovered_drivers.py)；[`experiments/cmame_round4/unified_comparison.py`](../experiments/cmame_round4/unified_comparison.py)
- **原始数据：** [`experiments/cmame_round2/00_driver_recovery/raw_data/checkerboard_manuscript_rerun.csv`](../experiments/cmame_round2/00_driver_recovery/raw_data/checkerboard_manuscript_rerun.csv)；[`experiments/cmame_round4/raw_data/combined_three_methods.csv`](../experiments/cmame_round4/raw_data/combined_three_methods.csv)
- **代表图表：** [`experiments/cmame_round4/figures/Figure_R4_combined_three_methods.pdf`](../experiments/cmame_round4/figures/Figure_R4_combined_three_methods.pdf)

<a id="f03"></a>

## F03 MacNeal–Harder M3 Hemisphere

- **原始研究目的：** 锁定敏感曲面上的位移、状态数及厚度趋势。
- **物理对象与设定：** R=10,18° cutout,E=6.825e7,nu=.3,t0=.04, quarter symmetry/交替径向力。
- **实验系列：** 标准厚度；raw vs midpoint 观测；载荷修正；厚度×网格；畸变
- **方法：** CB-Full；MidedgeTan；MITC4
- **主要数值结果：** raw n=12 CB .0909082、2.794% 对 .09352155；midpoint 重构 2.201%；n=16 raw .479%，MT 5.009%；最薄 n=16 CB 19.73% 差。
- **参考与比较基准：** 文献标准 .09352155；逐厚度加密外部参考
- **版本/轮次：** 原始论文；后续多轮
- **验证状态：** 选择性优势；参考与边界有限制
- **已知问题：** 旧数据混用 raw/midpoint 或载荷；MT n=12→16 非单调；极薄 M3 MT 稍优；跨方法对称运动学不同。
- **能够支持：** 局部精度/状态数和较 MITC4 弱的固定网格抑制。
- **不能支持：** 全局最优、无锁定或跨 3V 壳优势。
- **论文证据作用：** B,C,D
- **源程序：** [`experiments/cmame_round4/unified_comparison.py`](../experiments/cmame_round4/unified_comparison.py)；[`experiments/cmame_relative_value/run_scan.py`](../experiments/cmame_relative_value/run_scan.py)
- **原始数据：** [`experiments/cmame_round4/raw_data/combined_three_methods.csv`](../experiments/cmame_round4/raw_data/combined_three_methods.csv)；[`experiments/cmame_relative_value/results/thickness_scan.csv`](../experiments/cmame_relative_value/results/thickness_scan.csv)
- **代表图表：** [`experiments/cmame_relative_value/figures/thickness_mesh.pdf`](../experiments/cmame_relative_value/figures/thickness_mesh.pdf)

<a id="f04"></a>

## F04 大转动端矩悬臂

- **原始研究目的：** 检验位置型端矩、大转动与完整载荷路径。
- **物理对象与设定：** L=12,width=1,t=.1,E=1.2e6,nu=0；最高 2πEI/L 端矩。
- **实验系列：** 16×4/24×4/32×4 路径；最终状态 tangent 复核
- **方法：** CB-Full；解析常曲率圆弧
- **主要数值结果：** 32×4 最大转角 359.841°；路径 RMS 2.442%、最大 2.792%；膜能占 0.432%；tangent 相对差 1.64e-16。
- **参考与比较基准：** 解析常曲率路径
- **版本/轮次：** 原始论文；Round 4 重跑
- **验证状态：** 已复现全路径
- **已知问题：** 载荷路径和几何单元与其他壳的完全同条件比较未做。
- **能够支持：** 大转动平衡和精确非线性切线。
- **不能支持：** 普遍非线性壳优势。
- **论文证据作用：** B,E
- **源程序：** [`experiments/cmame_round4/cantilever.py`](../experiments/cmame_round4/cantilever.py)；[`experiments/cmame_round4/tangent_final_states.py`](../experiments/cmame_round4/tangent_final_states.py)
- **原始数据：** [`experiments/cmame_round4/raw_data/cantilever/mesh_32x4_raw_steps.csv`](../experiments/cmame_round4/raw_data/cantilever/mesh_32x4_raw_steps.csv)；[`experiments/cmame_round4/raw_data/tangent_final_states.json`](../experiments/cmame_round4/raw_data/tangent_final_states.json)
- **代表图表：** [`experiments/cmame_round4/figures/Figure_R4_nonlinear_paths.pdf`](../experiments/cmame_round4/figures/Figure_R4_nonlinear_paths.pdf)

<a id="f05"></a>

## F05 非线性交替力半球

- **原始研究目的：** 检验曲面参考 transport、非线性膜弯及完整压力/力路径。
- **物理对象与设定：** M3 quarter 几何；四交替径向力总量 400。
- **实验系列：** n=6/8/12/16 全路径；末态 sparse tangent
- **方法：** CB-Full；文献曲线
- **主要数值结果：** n=16 A/B 末态误差 2.759%/.149%，joint 路径 RMS 1.736%、最大 1.859%；末态 tangent 相对差 1.49e-15。
- **参考与比较基准：** 文献数字化曲线；B 端表值 3.8796
- **版本/轮次：** 原始论文；Round 4 重跑
- **验证状态：** 已复现；数字化参考有限精度
- **已知问题：** 参考曲线数字化约 .03 位移不确定性；不能与线性 M3 合并为同一物理响应。
- **能够支持：** 曲面大变形路径与 tangent。
- **不能支持：** 与成熟非线性壳同精度成本优势。
- **论文证据作用：** B
- **源程序：** [`experiments/cmame_round4/hemisphere.py`](../experiments/cmame_round4/hemisphere.py)；[`experiments/cmame_round4/tangent_final_states.py`](../experiments/cmame_round4/tangent_final_states.py)
- **原始数据：** [`experiments/cmame_round4/raw_data/hemisphere/mesh_16_raw_steps.csv`](../experiments/cmame_round4/raw_data/hemisphere/mesh_16_raw_steps.csv)；[`experiments/cmame_round4/raw_data/tangent_final_states.json`](../experiments/cmame_round4/raw_data/tangent_final_states.json)
- **代表图表：** [`experiments/cmame_round4/figures/Figure_R4_nonlinear_paths.pdf`](../experiments/cmame_round4/figures/Figure_R4_nonlinear_paths.pdf)

<a id="f06"></a>

## F06 开端圆柱壳 pull-out

- **原始研究目的：** 检验局部载荷、开自由端及强曲率梯度。
- **物理对象与设定：** L=10.35,R=4.953,t=.094,E=10.5e6,nu=.3125；1/8 域。
- **实验系列：** 8×8/12×12/16×16 全路径；C 点位移反转；末态 tangent
- **方法：** CB-Full；文献数字化路径
- **主要数值结果：** 16×16 末态 A/B/C 误差 8.094%/5.316%/1.861%，合并路径 RMS 4.483%；末态 tangent 相对差 4.95e-15。
- **参考与比较基准：** 文献图形路径
- **版本/轮次：** 原始论文；Round 4 重跑
- **验证状态：** 已复现；局部误差保留
- **已知问题：** 参考路径数字化；A 末态仍有 8% 差异。
- **能够支持：** 能追踪开边局部非线性反转。
- **不能支持：** 高精度强局部曲率或全局分支捕捉。
- **论文证据作用：** B,D
- **源程序：** [`experiments/cmame_round4/open_cylinder.py`](../experiments/cmame_round4/open_cylinder.py)；[`experiments/cmame_round4/tangent_final_states.py`](../experiments/cmame_round4/tangent_final_states.py)
- **原始数据：** [`experiments/cmame_round4/raw_data/open_cylinder/open_cylinder_pullout/summary.csv`](../experiments/cmame_round4/raw_data/open_cylinder/open_cylinder_pullout/summary.csv)；[`experiments/cmame_round4/raw_data/open_cylinder/open_cylinder_pullout/mesh_16x16_raw_path.csv`](../experiments/cmame_round4/raw_data/open_cylinder/open_cylinder_pullout/mesh_16x16_raw_path.csv)
- **代表图表：** [`experiments/cmame_round4/figures/Figure_R4_nonlinear_paths.pdf`](../experiments/cmame_round4/figures/Figure_R4_nonlinear_paths.pdf)

<a id="f07"></a>

## F07 闭端 pinched cylinder 折叠分支探索

- **原始研究目的：** 检验严重 snap-through 和分支延拓能力。
- **物理对象与设定：** R=100,L=200,t=1,E=30000,nu=.3；刚性端隔板；最高 12000 双力。
- **实验系列：** 载荷控制 6/8/12 网格；SVD null-vector pseudo-arclength；分支诊断
- **方法：** 历史 CB-Full 研究驱动
- **主要数值结果：** 三网格载荷控制遇折点；弧长探针跨过局部载荷反转但未连接文献后屈曲支，状态 NONLINEAR_SOLVER_LIMITED。
- **参考与比较基准：** 文献数字化后屈曲曲线
- **版本/轮次：** 原始历史仓库；原稿未纳入
- **验证状态：** 失败/未纳入论文
- **已知问题：** 不得把末次已收敛状态当作全路径成功；历史输出未迁入当前 repo。
- **能够支持：** 明确非线性分支求解局限。
- **不能支持：** 完成严重 snap-through benchmark。
- **论文证据作用：** D
- **源程序：** `/home/song/projects/CB-shell/experiments/checkerboard_nonlinear_benchmarks/pinched_cylinder.py`；`/home/song/projects/CB-shell/experiments/checkerboard_nonlinear_benchmarks/pinched_solver_refinement.py`
- **原始数据：** `/home/song/projects/CB-shell/outputs/checkerboard_nonlinear_benchmarks/pinched_cylinder/results.json`；`/home/song/projects/CB-shell/outputs/checkerboard_nonlinear_benchmarks/pinched_cylinder/solver_refinement_results.json`；`/home/song/projects/CB-shell/outputs/checkerboard_nonlinear_benchmarks/pinched_cylinder/solver_refinement_branch_diagnostics.csv`
- **代表图表：** 未找到图表；请按数据/报告核实

<a id="i01"></a>

## I01 九节点作动力平板形状控制

- **原始研究目的：** 验证线性状态下力设计、隐式/伴随梯度和非唯一性。
- **物理对象与设定：** 12×8 平板；9 z 力；99 观测；297 自由 CB 状态。
- **实验系列：** adjoint 中心差分；CB adjoint/FD；MITC4 FD own/same target；后期重跑
- **方法：** CB-Full；MITC4
- **主要数值结果：** 原记录 CB adjoint 13 次迭代/79 解/1.05 s/目标相对误差 6.51e-5；FD 438 解/2.08 s；力恢复相对误差 .174；后期同方法重跑 68/391 解。
- **参考与比较基准：** 生成的目标位移；中心差分梯度
- **版本/轮次：** Round 5 新增；工程轮次回归
- **验证状态：** 同模型梯度强验证；力识别非唯一
- **已知问题：** MITC4 own-target误差 .0026% 比 CB 低；FD/adjoint与求解器/重建混合；历史 MITC 12×12 行失配已修正。
- **能够支持：** 一致梯度精度及同模型优化工作量收益。
- **不能支持：** CB 单元自身速度优势或力唯一反演。
- **论文证据作用：** E
- **源程序：** [`experiments/cmame_round5/actuation_design.py`](../experiments/cmame_round5/actuation_design.py)；[`experiments/cmame_round5/fair_fd_adjoint.py`](../experiments/cmame_round5/fair_fd_adjoint.py)；[`experiments/cmame_round5/fd_vs_fd.py`](../experiments/cmame_round5/fd_vs_fd.py)
- **原始数据：** [`experiments/cmame_round5/raw_data/actuation_design.json`](../experiments/cmame_round5/raw_data/actuation_design.json)；[`experiments/cmame_round5/raw_data/adjoint_check.json`](../experiments/cmame_round5/raw_data/adjoint_check.json)；[`experiments/cmame_round5/raw_data/fair_fd_adjoint.json`](../experiments/cmame_round5/raw_data/fair_fd_adjoint.json)
- **代表图表：** [`experiments/cmame_round5/figures/Figure_R5_actuation_design.pdf`](../experiments/cmame_round5/figures/Figure_R5_actuation_design.pdf)

<a id="i02"></a>

## I02 12 参数参考形状优化

- **原始研究目的：** 测试参考度量/曲率随初始几何变化的线性化设计通路。
- **物理对象与设定：** 预弯板 4×4 高度控制格内 12 参数；297 自由状态、99 观测。
- **实验系列：** 参考形状参数与目标位移拟合；优化历史
- **方法：** CB-Full 参考 tangent
- **主要数值结果：** 20 次迭代达上限；位移相对误差 1.60e-4，参数高度相对误差 .968，显示非唯一。
- **参考与比较基准：** 由已知参考几何生成的位移目标
- **版本/轮次：** Round 5 新增
- **验证状态：** 目标拟合有效；参数恢复失败
- **已知问题：** 非线性几何平衡与独立同目标传统壳比较未做；优化到 iteration cap。
- **能够支持：** 参考几何进入力学灵敏度与逆向目标。
- **不能支持：** 几何参数唯一可辨识。
- **论文证据作用：** E,D
- **源程序：** [`experiments/cmame_round5/shape_design.py`](../experiments/cmame_round5/shape_design.py)
- **原始数据：** [`experiments/cmame_round5/raw_data/shape_design.json`](../experiments/cmame_round5/raw_data/shape_design.json)；[`experiments/cmame_round5/raw_data/shape_design_history.csv`](../experiments/cmame_round5/raw_data/shape_design_history.csv)
- **代表图表：** [`experiments/cmame_round5/figures/Figure_R5_shape_design.pdf`](../experiments/cmame_round5/figures/Figure_R5_shape_design.pdf)

<a id="i03"></a>

## I03 90°级大转动参考几何反演

- **原始研究目的：** 将初始曲率参数与非线性平衡和目标位姿耦合。
- **物理对象与设定：** 端矩 M/Mmax=.25 的条带；2 个初始高度参数；46 自由状态。
- **实验系列：** 两参数曲率恢复；非线性平衡路径；目标 tip
- **方法：** CB-Full 非线性 solve
- **主要数值结果：** 8 次迭代；tip 误差 1.86e-6；真高度 [0,.375,1.5]，恢复 [0,.374998,1.499998]。
- **参考与比较基准：** 已知初始曲率生成的目标
- **版本/轮次：** Round 5 新增
- **验证状态：** 已完成；小参数合成逆问题
- **已知问题：** 历史 DESIGN_REPORT 顶部“demo 4 running”已过时；合成目标不证明真实结构识别唯一性。
- **能够支持：** 参考几何参数/大转动/非线性平衡耦合。
- **不能支持：** 一般大规模逆向最优或相对壳效率。
- **论文证据作用：** E
- **源程序：** [`experiments/cmame_round5/large_rotation_design.py`](../experiments/cmame_round5/large_rotation_design.py)
- **原始数据：** [`experiments/cmame_round5/raw_data/large_rotation_design.json`](../experiments/cmame_round5/raw_data/large_rotation_design.json)；[`experiments/cmame_round5/raw_data/large_rotation_design_history.csv`](../experiments/cmame_round5/raw_data/large_rotation_design_history.csv)
- **代表图表：** [`experiments/cmame_round5/figures/Figure_R5_large_rotation_design.pdf`](../experiments/cmame_round5/figures/Figure_R5_large_rotation_design.pdf)

<a id="i04"></a>

## I04 早期飞机蒙皮九力反演

- **原始研究目的：** 将九力平板流程迁移到弯曲蒙皮并采用实际单位目标。
- **物理对象与设定：** Round 6 无量纲曲面；Round 7 NACA 后缘 .6m×.3m、2mm 铝、5kPa。
- **实验系列：** R6 各自生成目标；R7 外给 +5mm 目标；adjoint/黑箱 FD
- **方法：** CB-Full；MITC4
- **主要数值结果：** R7 CB 10 迭代/44 解/.19s/3.63% 目标误差、最大 .250mm；MITC4 8.98s/3.55%；均非真实作动机构。
- **参考与比较基准：** 自身或外给目标；传统模型交叉
- **版本/轮次：** Round 6/7 后续新增
- **验证状态：** 历史/已被工程布局替代
- **已知问题：** 九个节点力非九个实际作动器；材料/边界/载荷代表性不足；黑箱 MITC 计时不公平。
- **能够支持：** 力反演扩展可计算。
- **不能支持：** 真实飞机作动布局或壳单元快 48×。
- **论文证据作用：** F,D
- **源程序：** [`experiments/cmame_round6_applications/aircraft_skin.py`](../experiments/cmame_round6_applications/aircraft_skin.py)；[`experiments/cmame_round7_aircraft/aircraft_actual.py`](../experiments/cmame_round7_aircraft/aircraft_actual.py)
- **原始数据：** [`experiments/cmame_round6_applications/raw_data/aircraft_skin.json`](../experiments/cmame_round6_applications/raw_data/aircraft_skin.json)；[`experiments/cmame_round7_aircraft/raw_data/aircraft_actual.json`](../experiments/cmame_round7_aircraft/raw_data/aircraft_actual.json)
- **代表图表：** [`experiments/cmame_round7_aircraft/figures/Figure_R7_aircraft_actual_profile.pdf`](../experiments/cmame_round7_aircraft/figures/Figure_R7_aircraft_actual_profile.pdf)

<a id="i05"></a>

## I05 飞机五候选作动布局与行程

- **原始研究目的：** 找两工况可行布局、反力/行程及 31 布局筛选成本。
- **物理对象与设定：** NACA 后缘 .6m 弦/1m 展、.5mm 铝、位置梁；5 候选；4mm/80N 限。
- **实验系列：** 31 线性布局；单/双点非线性本地审计；[0,2,4] 非线性；参考与复用性能
- **方法：** CB-Full；ShellMITC4；ASD corotational/frame
- **主要数值结果：** [0,2,4] 两工况 CB RMS/max .146/.412 与 .168/.599 mm；20×40 ASD .181/.412 与 .206/.644 mm；1/2 点未证全局不可行。
- **参考与比较基准：** 同输入加密 ASD；形状容限 .25/.8mm
- **版本/轮次：** engineering_inverse 后续新增
- **验证状态：** 有可行布局；最少数量未证明
- **已知问题：** 位置梁与转动梁有限网格运动学不同；线性布局不稳定；不是完整31布局非线性全局优化。
- **能够支持：** 两工况可行设计和局部竞争布局诊断。
- **不能支持：** 三处全局最少或完整非线性搜索加速。
- **论文证据作用：** F,D
- **源程序：** [`experiments/engineering_inverse/aircraft_skin.py`](../experiments/engineering_inverse/aircraft_skin.py)；[`experiments/engineering_inverse/nonlinear_audit.py`](../experiments/engineering_inverse/nonlinear_audit.py)；[`experiments/engineering_inverse/aircraft_nonlinear_reference.py`](../experiments/engineering_inverse/aircraft_nonlinear_reference.py)
- **原始数据：** [`experiments/engineering_inverse/results/aircraft_layout_summary.csv`](../experiments/engineering_inverse/results/aircraft_layout_summary.csv)；[`experiments/engineering_inverse/results/aircraft_nonlinear_audit_12x24.json`](../experiments/engineering_inverse/results/aircraft_nonlinear_audit_12x24.json)；[`experiments/engineering_inverse/results/aircraft_corotational_reference.json`](../experiments/engineering_inverse/results/aircraft_corotational_reference.json)
- **代表图表：** [`experiments/engineering_inverse/figures/aircraft_layout_curve.svg`](../experiments/engineering_inverse/figures/aircraft_layout_curve.svg)

<a id="i06"></a>

## I06 早期汽车曲板九力 gap/flush 反演

- **原始研究目的：** 在简化曲板中反演九个夹具力。
- **物理对象与设定：** 1×.6 曲板，厚 .02，9 集中力，11 测量点。
- **实验系列：** CB adjoint 与 MITC4 FD 各自生成目标
- **方法：** CB-Full；MITC4
- **主要数值结果：** CB 24 次迭代/89 解/1.00s/目标误差 1.28e-4；MITC4 7.91s/2.55e-4；力场不唯一。
- **参考与比较基准：** 各模型自身力场生成的目标
- **版本/轮次：** Round 6 后续新增
- **验证状态：** 历史/被安装调整替代
- **已知问题：** 不是制造偏差+安装点调整；两模型目标不同，时间混入 FD/重建差别。
- **能够支持：** 一般曲板力反演可计算。
- **不能支持：** 真实汽车弹性装配优化。
- **论文证据作用：** F,D
- **源程序：** [`experiments/cmame_round6_applications/automotive_panel.py`](../experiments/cmame_round6_applications/automotive_panel.py)
- **原始数据：** [`experiments/cmame_round6_applications/raw_data/automotive_panel.json`](../experiments/cmame_round6_applications/raw_data/automotive_panel.json)
- **代表图表：** [`experiments/cmame_round6_applications/figures/Figure_R6_automotive_panel.pdf`](../experiments/cmame_round6_applications/figures/Figure_R6_automotive_panel.pdf)

<a id="i07"></a>

## I07 汽车覆盖件非线性安装调整

- **原始研究目的：** 通过八个安装行程修正 gap/flush 并直接约束反力。
- **物理对象与设定：** .80×.45m 曲板、.7mm 钢、四有限支架、两 datum 加一 x locator；±1.5mm/300N/.6mm 限。
- **实验系列：** 硬区域旧模型；有限支架；共同 locator/gauge；CB consistent/FD 与 ASD-FD；KKT
- **方法：** CB-Full；ASDShellQ4 corotational；线性 MITC 诊断
- **主要数值结果：** 最终 CB-consistent 35 迭代/50 平衡/200.3s/设计目标 .031260202；共同 ASD 24×20 误差 .2569mm/最大支架 88.55N；ASD-FD 3676.5s但参考目标更低。
- **参考与比较基准：** 24×20 ASD 共同参考；独立 KKT/梯度复核
- **版本/轮次：** engineering_inverse 后续新增并多次修正
- **验证状态：** 研究容限通过；物理支架未标定
- **已知问题：** 旧硬区域 1245N 失败；单 x gauge 漏约束曾伪装多分支；ASD 最优驻点未等精度认证；反力矩账未闭合。
- **能够支持：** 条件化工程可用性和同 CB 梯度收益。
- **不能支持：** 真实车型预测、同最优目标的跨壳加速。
- **论文证据作用：** E,F,D
- **源程序：** [`experiments/engineering_inverse/automotive_spring_design.py`](../experiments/engineering_inverse/automotive_spring_design.py)；[`experiments/engineering_inverse/efficiency_automotive.py`](../experiments/engineering_inverse/efficiency_automotive.py)；[`experiments/engineering_inverse/validate_automotive_reliability.py`](../experiments/engineering_inverse/validate_automotive_reliability.py)
- **原始数据：** [`experiments/engineering_inverse/results/automotive_efficiency_cb_consistent_12x8_zero.json`](../experiments/engineering_inverse/results/automotive_efficiency_cb_consistent_12x8_zero.json)；[`experiments/engineering_inverse/results/automotive_efficiency_common_reference_zero.json`](../experiments/engineering_inverse/results/automotive_efficiency_common_reference_zero.json)；[`experiments/engineering_inverse/results/automotive_reliability_cb_consistent_zero.json`](../experiments/engineering_inverse/results/automotive_reliability_cb_consistent_zero.json)
- **代表图表：** [`experiments/engineering_inverse/figures/final_automotive_optimization_cost.svg`](../experiments/engineering_inverse/figures/final_automotive_optimization_cost.svg)

<a id="p01"></a>

## P01 精确稀疏＋秩一切线与边界存储

- **原始研究目的：** 确认全局 XM 仅秩一、非线性局部项/边界重构可稀疏存储。
- **物理对象与设定：** 板/M3 n=16–128；有限变形 tangent 与 dense AD 比较。
- **实验系列：** Round 1 dense/reference-linear 原型；Round 2 exact nonlinear sparse border；内存/nnz/solve
- **方法：** CB-0/M/Full；历史 dense 对照
- **主要数值结果：** M3 n=128 CB-Full 49923 raw states、3961184 augmented nnz、11.49s assembly、16.39s factor、2073MiB RSS；末态 tangent 差 1e-16–1e-15。
- **参考与比较基准：** 同能量 dense AD Hessian；历史 boundary 实现
- **版本/轮次：** 理论原稿；Round 1/2 生产化
- **验证状态：** 精确结构已验证
- **已知问题：** Round 1 n=64/128 内存投影非实测；CB Python/AD 内存不低于 C++ MT；factors expensive。
- **能够支持：** 稀疏局部+单秩一精确结构。
- **不能支持：** CPU 或 RSS 总体优于传统壳。
- **论文证据作用：** A,D
- **源程序：** [`experiments/cmame_additional/run_primary.py`](../experiments/cmame_additional/run_primary.py)；[`experiments/cmame_round2/run_exactness.py`](../experiments/cmame_round2/run_exactness.py)；[`experiments/cmame_round2/run_scaling.py`](../experiments/cmame_round2/run_scaling.py)；[`experiments/cmame_round2/run_boundary_equivalence.py`](../experiments/cmame_round2/run_boundary_equivalence.py)
- **原始数据：** [`experiments/cmame_additional/03_scaling/raw_data/rank_one_audit.json`](../experiments/cmame_additional/03_scaling/raw_data/rank_one_audit.json)；[`experiments/cmame_round2/02_low_rank_xm/raw_data/exactness.json`](../experiments/cmame_round2/02_low_rank_xm/raw_data/exactness.json)；[`experiments/cmame_round2/03_scaling/raw_data/scaling.csv`](../experiments/cmame_round2/03_scaling/raw_data/scaling.csv)；[`experiments/cmame_round2/01_sparse_boundary/raw_data/boundary_storage.json`](../experiments/cmame_round2/01_sparse_boundary/raw_data/boundary_storage.json)
- **代表图表：** [`experiments/cmame_round2/03_scaling/figures/Figure_C3_nnz_scaling.pdf`](../experiments/cmame_round2/03_scaling/figures/Figure_C3_nnz_scaling.pdf)

<a id="p02"></a>

## P02 跨方法精度—状态数—时间对照

- **原始研究目的：** 在共同观测和参考上区分状态效率与时间。
- **物理对象与设定：** 板/屋顶/M3，多网格和逐厚度参考。
- **实验系列：** recovered CB/MT；MITC4；同目标精度阈值；CPU/RSS
- **方法：** CB-Full；MidedgeTan；MITC4
- **主要数值结果：** M3 标准厚度 1% scalar error 首测 CB 867、MITC4 6534 raw states；屋顶/板 MITC 更好；roof n=24 总时间 CB/MT/MITC 1.50/.137/.230s。
- **参考与比较基准：** 解析/文献值；逐厚度外部加密 surrogate
- **版本/轮次：** 原稿 MT；Round 3/4 MITC；相对研究
- **验证状态：** 混合结果；限定标量/实现
- **已知问题：** 无可靠同条件 BSQ 3V 基线；部分 MT M3 非单调；CPU 跨不同软件栈。
- **能够支持：** M3 特定精度/状态数收益及无普遍 CPU 优势。
- **不能支持：** 全问题精度/状态优势。
- **论文证据作用：** C,D
- **源程序：** [`experiments/cmame_round3/build_final_comparison.py`](../experiments/cmame_round3/build_final_comparison.py)；[`experiments/cmame_round4/combined_comparison.py`](../experiments/cmame_round4/combined_comparison.py)；[`experiments/cmame_relative_value/analyze.py`](../experiments/cmame_relative_value/analyze.py)
- **原始数据：** [`experiments/cmame_round3/02_final_comparison/raw_data/final_comparative_evidence.csv`](../experiments/cmame_round3/02_final_comparison/raw_data/final_comparative_evidence.csv)；[`experiments/cmame_relative_value/results/target_states.csv`](../experiments/cmame_relative_value/results/target_states.csv)；[`experiments/cmame_relative_value/results/runtime_summary.csv`](../experiments/cmame_relative_value/results/runtime_summary.csv)
- **代表图表：** [`experiments/cmame_relative_value/figures/accuracy_states.pdf`](../experiments/cmame_relative_value/figures/accuracy_states.pdf)

<a id="p03"></a>

## P03 飞机 31 布局共享响应算子对照

- **原始研究目的：** 分离逐布局重建与通用矩阵复用的成本。
- **物理对象与设定：** 固定 5 候选、31 布局、两工况；12×24/16×32/20×40。
- **实验系列：** CB direct/shared；MITC direct/shared；三次中位；同解一致性
- **方法：** CB-Full；ShellMITC4
- **主要数值结果：** 最终 shared n=12×24/16×32/20×40 CB 1.688/2.310/3.890s，MITC .497/.873/2.033s；MITC 较快。
- **参考与比较基准：** 相同线性布局任务；同共享策略
- **版本/轮次：** engineering_inverse 后续新增
- **验证状态：** 公平共享筛选；仅线性
- **已知问题：** 较早 8× 声称混入 OpenSees 重建；加强构件运动学仍不同；不得称完整非线性搜索速度。
- **能够支持：** 共享策略属通用技术；当前 CB 总时间未胜。
- **不能支持：** CB 单元计算速度优越。
- **论文证据作用：** D
- **源程序：** [`experiments/engineering_inverse/efficiency_aircraft.py`](../experiments/engineering_inverse/efficiency_aircraft.py)；[`experiments/engineering_inverse/mitc4_shared.py`](../experiments/engineering_inverse/mitc4_shared.py)；[`experiments/engineering_inverse/validate_aircraft_efficiency.py`](../experiments/engineering_inverse/validate_aircraft_efficiency.py)
- **原始数据：** [`experiments/engineering_inverse/results/tangent_aircraft_final.csv`](../experiments/engineering_inverse/results/tangent_aircraft_final.csv)；[`experiments/engineering_inverse/results/aircraft_efficiency_four_methods.json`](../experiments/engineering_inverse/results/aircraft_efficiency_four_methods.json)
- **代表图表：** [`experiments/engineering_inverse/figures/final_aircraft_shared_screen.svg`](../experiments/engineering_inverse/figures/final_aircraft_shared_screen.svg)

<a id="p04"></a>

## P04 局部 Hessian 瓶颈与缓存实验

- **原始研究目的：** 测量 CB 精确切线为何慢并测试完全等价优化。
- **物理对象与设定：** 20×40 飞机 profile；三个布局网格；汽车固定输入和完整路径。
- **实验系列：** 147 graph colours；legacy/cached AD batch；力学同解回归
- **方法：** CB-Full 旧/新切线策略
- **主要数值结果：** 20×40 tangent 3.059s，其中 coloured HVP 2.906s；缓存后飞机 3.876→3.890s无稳定收益；汽车匹配路径 266.3→200.3s。
- **参考与比较基准：** 同能量/残差/Hessian 及布局输出
- **版本/轮次：** engineering_inverse 最终审计
- **验证状态：** 实现等价；收益混合
- **已知问题：** 汽车总时间改变也含求解路径；飞机无稳定速度收益。
- **能够支持：** 瓶颈来自 AD 局部 Hessian。
- **不能支持：** 3V 必然更快。
- **论文证据作用：** D
- **源程序：** [`experiments/engineering_inverse/profile_hessian.py`](../experiments/engineering_inverse/profile_hessian.py)；[`experiments/engineering_inverse/make_tangent_reliability_figures.py`](../experiments/engineering_inverse/make_tangent_reliability_figures.py)
- **原始数据：** [`experiments/engineering_inverse/results/hessian_profile_20x40_legacy.json`](../experiments/engineering_inverse/results/hessian_profile_20x40_legacy.json)；[`experiments/engineering_inverse/results/tangent_aircraft_final.csv`](../experiments/engineering_inverse/results/tangent_aircraft_final.csv)；[`experiments/engineering_inverse/results/automotive_tangent_before_after.json`](../experiments/engineering_inverse/results/automotive_tangent_before_after.json)
- **代表图表：** [`experiments/engineering_inverse/figures/final_aircraft_tangent_assembly.svg`](../experiments/engineering_inverse/figures/final_aircraft_tangent_assembly.svg)

<a id="p05"></a>

## P05 汽车三种完整非线性优化成本与梯度归因

- **原始研究目的：** 比较同问题 CB 一致/FD 与 ASD-FD 的工作量、终点和共同参考。
- **物理对象与设定：** 12×8 设计网格、零行程初值、共同 x locator；24×20 ASD 验收。
- **实验系列：** CB consistent；CB 8变量中心 FD；ASD corotational FD；KKT/branch
- **方法：** CB-Full；ASDShellQ4
- **主要数值结果：** 完整时间 200.3/2509.4/3676.5s；CB 同模型一致梯度快 12.53×；CB/ASD 最终参考目标 .037786/.022173，不是同最优目标。
- **参考与比较基准：** 同共同参考工程容限；CB 独立 KKT
- **版本/轮次：** engineering_inverse 最终审计
- **验证状态：** 同 CB 梯度收益强；跨壳限于达容限
- **已知问题：** ASD 终点更低参考目标且站点精度未等同；数值栈与非线性路径不同。
- **能够支持：** 一致灵敏度工作流价值。
- **不能支持：** 壳单元独有 18× 加速或同最优目标加速。
- **论文证据作用：** E,D
- **源程序：** [`experiments/engineering_inverse/efficiency_automotive.py`](../experiments/engineering_inverse/efficiency_automotive.py)；[`experiments/engineering_inverse/efficiency_automotive_reference.py`](../experiments/engineering_inverse/efficiency_automotive_reference.py)
- **原始数据：** [`experiments/engineering_inverse/results/tangent_automotive_optimization_final.csv`](../experiments/engineering_inverse/results/tangent_automotive_optimization_final.csv)；[`experiments/engineering_inverse/results/automotive_efficiency_common_reference_zero.json`](../experiments/engineering_inverse/results/automotive_efficiency_common_reference_zero.json)
- **代表图表：** [`experiments/engineering_inverse/figures/final_automotive_optimization_cost.svg`](../experiments/engineering_inverse/figures/final_automotive_optimization_cost.svg)

<a id="x01"></a>

## X01 旧仓库替代膜/补全形式探索

- **原始研究目的：** 寻找无盲模、低锁定的统一膜形式。
- **物理对象与设定：** B-only/W 观测、intrinsic connection、mixed/assumed strain 候选；M1/M2/M3。
- **实验系列：** full-energy 反例；locking-free extension；最终信息层级搜索
- **方法：** 历史候选 CB 能量，非当前 CB-Full
- **主要数值结果：** 最终 verdict NO_UNIFIED_CHECKERBOARD_MEMBRANE_FORMULATION_FOUND；connection nonlinear M3 n=8 5.1096、n=12 6.7254 对参考 3.8796。
- **参考与比较基准：** 候选 kernel/patch/非线性 M3
- **版本/轮次：** 外部历史仓库；未迁入正式内核
- **验证状态：** 历史失败；与当前正式能量隔离
- **已知问题：** 不能把不同试验能量混作当前 CB-Full；数个旧假设已被反例否定。
- **能够支持：** 记录方法边界及失败机制。
- **不能支持：** 当前生产能量已解决所有膜锁定。
- **论文证据作用：** D
- **源程序：** `/home/song/projects/CB-shell/experiments/checkerboard_membrane_final`；`/home/song/projects/CB-shell/experiments/checkerboard_locking_free_extension`
- **原始数据：** `/home/song/projects/CB-shell/outputs/checkerboard_membrane_final/decision_summary.json`；`/home/song/projects/CB-shell/outputs/checkerboard_locking_free_extension/prototype_summary.json`
- **代表图表：** `/home/song/projects/CB-shell/outputs/intrinsic_m_general_locking/m3_relative_error_vs_thickness.pdf`

<a id="x02"></a>

## X02 旧仓库 matrix-free/GPU 批量计算

- **原始研究目的：** 验证冻结物理的 HVP、CPU Krylov、GPU 多同拓扑场景吞吐。
- **物理对象与设定：** M1/M2/M3 多场景；GPU 同拓扑 batch。
- **实验系列：** matrix-free CPU；GPU HVP；多场景与失败区间
- **方法：** 历史 CB 后端；独立硬件
- **主要数值结果：** M2 n=30 CPU 24.681→17.461s；n=32/B=256 GPU 717.4 solves/s，对 CPU matrix-free sequential 169.5×，限定同拓扑批量。
- **参考与比较基准：** 同物理 HVP 与显式切线；硬件环境
- **版本/轮次：** 外部历史仓库；非本篇当前 CPU 结论
- **验证状态：** 独立历史结果
- **已知问题：** GPU 批量吞吐不证明单次壳求解快；高精度 M3 大网格受预条件限制。
- **能够支持：** 单独的并行计算研究线索。
- **不能支持：** 当前 CPU 元素级优势或普遍 GPU 优势。
- **论文证据作用：** D
- **源程序：** `/home/song/projects/CB-shell/experiments/checkerboard_matrix_free_gpu`；`/home/song/projects/CB-shell/experiments/checkerboard_multiscenario_gpu`
- **原始数据：** `/home/song/projects/CB-shell/outputs/checkerboard_multiscenario_gpu/14_regime_map.csv`；`/home/song/projects/CB-shell/outputs/checkerboard_multiscenario_gpu/13_dataset_generation.csv`
- **代表图表：** `/home/song/projects/CB-shell/outputs/checkerboard_matrix_free_gpu/15_paper_figures/figure_05_batch_throughput.pdf`

<a id="x03"></a>

## X03 旧仓库 FPC 洛伦兹/接触柔性片

- **原始研究目的：** 研究柔性印刷电路、导线洛伦兹力和接触。
- **物理对象与设定：** 100mm 柔性板及两端口四支路细条；地面接触。
- **实验系列：** 节点/积分点接触；后屈曲分支；三角壳替代；细网格门槛
- **方法：** 旧 DDG/FPC 耦合；三角壳候选
- **主要数值结果：** 窄条旧 CB 奇偶节点零能机制导致高载荷形貌撤回；替代三角 hinge 弯曲绝对系数不对，Mechanics Gate NOT VALIDATED。
- **参考与比较基准：** 力/矩守恒和网格对照
- **版本/轮次：** 外部历史项目；非本论文实验
- **验证状态：** 独立且含失败结果
- **已知问题：** 不能将其高载荷形状或三角替代后端作为当前 CB-Full 方法证据。
- **能够支持：** 历史独立支线中的离散风险提醒。
- **不能支持：** 本篇材料壳精度或工程优势。
- **论文证据作用：** D
- **源程序：** `/home/song/projects/CB-shell/PROJECT_HISTORY.md`；`/home/song/projects/CB-shell/experiments/archive`
- **原始数据：** `/home/song/projects/CB-shell/outputs/two_port_parallel_strip/summary.csv`；`/home/song/projects/CB-shell/outputs/triangular_shell_validation`
- **代表图表：** 未找到图表；请按数据/报告核实

<a id="x04"></a>

## X04 4D 打印双层自然形式逆向

- **原始研究目的：** 从目标自然形式到可打印相位、双层本征应变和平衡形状。
- **物理对象与设定：** 双层打印材料；相位/路径/自然曲率；vendored CB commit 0699638。
- **实验系列：** 几何可达性、打印参数非唯一、目标分支和前向验收
- **方法：** inverse4d + vendored CB Shell
- **主要数值结果：** 同层厚度/刚度及方向控制下 tr(kappa0)=0，圆柱/球目标一般不能精确匹配；多种逆向记录在独立 docs/outputs。
- **参考与比较基准：** 独立材料可达性与前向模型
- **版本/轮次：** 外部独立仓库；不计本文逆向案例
- **验证状态：** 独立项目，未逐项审计
- **已知问题：** 外部仓库有自己的 AGENTS/研究范围；不能自动并入本篇 CB 壳验证。
- **能够支持：** 仅提示潜在后续耦合应用。
- **不能支持：** 本项目已完成的非工程逆向案例。
- **论文证据作用：** E,D
- **源程序：** `/home/song/projects/inverse_design/4d_printing_inverse_design/README.md`；`/home/song/projects/inverse_design/4d_printing_inverse_design/examples/run_benchmarks.py`
- **原始数据：** `/home/song/projects/inverse_design/4d_printing_inverse_design/outputs/benchmarks`；`/home/song/projects/inverse_design/4d_printing_inverse_design/docs/INVERSE_4D_EXPLORATION.md`
- **代表图表：** 未找到图表；请按数据/报告核实
