# 旧仓库实验目录逐项归属

来源为邻接 `/home/song/projects/CB-shell/experiments/` 工作树；`v1.0-arxiv` 标签是否包含**具体文件**请查 [`HISTORICAL_ARTIFACT_MANIFEST.csv`](HISTORICAL_ARTIFACT_MANIFEST.csv)，不能由目录名推断。以下按目录合并到本清查案例，既防遗漏，也不把同一个 M3/屋顶重复算成多个独立案例。

| 历史实验目录 | 对应案例或支线 | 归属理由与限制 |
|---|---|---|
| `archive` | X03/待逐文件归属 | 旧 FPC 与原型存档，内容混合；不得整体计作正式 CB-Full |
| `checkerboard_benchmark` | F01–F03/M01–M04 | 早期经典基准和机制；以原稿与 Round 2/4 复核为准 |
| `checkerboard_bform` | M01 | 第二基本形式几何研究 |
| `checkerboard_boundary_closure` | P01/M04 | W/边界重构与存储；后期生产实现另在 Round 2 |
| `checkerboard_cross_method_gpu` | X02 | 跨后端/GPU 独立性能线 |
| `checkerboard_curved_reference` | M04 | 曲面参考/transport 变体 |
| `checkerboard_curved_transport` | M04 | 曲面补全客观 transport |
| `checkerboard_final_benchmark` | F01–F06/P02 | 原稿前后经典与非线性基准汇总，不另算新物理案例 |
| `checkerboard_full_energy` | M03/M04 | Koiter 膜弯能及补全研究 |
| `checkerboard_locking_free_extension` | X01/M05 | 历史试验形式及薄极限；非当前正式能量 |
| `checkerboard_locking_free_membrane` | X01 | 历史候选膜形式；失败也保留 |
| `checkerboard_matrix_free_gpu` | X02 | HVP/matrix-free/GPU 独立性能线 |
| `checkerboard_membrane_completion` | M03/X01 | 膜盲模态及候选补全的历史阶段 |
| `checkerboard_membrane_final` | X01 | 后续替代膜形式信息层级搜索，不能混入正式 CB-Full |
| `checkerboard_membrane_locking` | X01/M05 | 薄极限与候选膜形式筛选 |
| `checkerboard_mesh_thickness` | M05 | 厚度—网格系列，不另算一个新的结构案例 |
| `checkerboard_multiscenario_gpu` | X02 | 多同拓扑场景吞吐，独立硬件/计算域 |
| `checkerboard_nonlinear_benchmarks` | F04–F07 | 悬臂、半球、开端圆柱及失败闭端 pinched cylinder |
| `checkerboard_open_cylinder_pullout` | F06 | 开端圆柱更早入口；同物理案例 |
| `checkerboard_operator_audit` | M02–M04 | 代数/符号算子审计，机制验证系列 |
| `checkerboard_paper_benchmark` | M01–M04/F01–F06/P02 | 原稿数值复现；工作树输出不一定在原稿标签 |
| `checkerboard_theory_audit` | M02–M04 | 符号理论审计，不另算加载案例 |
| `checkerboard_thin_limit_compatibility` | M05/X01 | 近等距、薄极限、兼容投影；含候选实现 |
| `checkerboard_unified_membrane` | X01 | 历史统一膜形式尝试 |
| `checkerboard_xw_sensitivity` | M04/M06 | XW 响应/系数敏感性 |
| `intrinsic_m_beta` | X01/M06 | 历史 intrinsic 候选的 beta 扫描，不是生产 beta_M 的独立调参证明 |
| `intrinsic_m_completion` | X01/M03 | 历史 intrinsic 试验补全 |
| `intrinsic_m_general_locking` | X01/M05 | 历史替代膜形式锁定实验 |
| `m3_ab_xa_locking_audit` | X01/M05 | 历史候选能量的 M3 锁定审计 |
| `phase_kernel_isometry_theorem` | M03/M05 | 相位核/等距机制验证，非新壳加载案例 |

**边界：** 旧仓库 30 个实验子目录均列于上表。X01 汇总的是**一条替代能量/无锁定探索研究线**，其中每个候选公式的数值扫描是实验系列，不自动增加当前正式 CB-Full 论文的案例数。`archive` 内容异质，逐文件见历史 manifest；其内部个别项目不能仅靠目录名证明通过力学验证。旧仓库没有找到已完成的额外非工程逆向设计入口；本项目的 I01–I03 均在当前 Round 5。
