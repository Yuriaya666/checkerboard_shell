# 最小 forward 依赖链

本次抽取针对“给定 checkerboard geometry + material + load，计算 deformation”这一条主线。旧论文仓库中的 benchmark、输出和实验脚本没有进入新工程。

## 来源与对应模块

| forward 阶段 | 新工程模块 | 作用 |
| --- | --- | --- |
| geometry construction / mesh generation | `geometry/mesh.py`（结构化矩形网格构造） | 建立参数化 `(ny+1)×(nx+1)` raw positions、quad 拓扑和四条边界 |
| B/W midpoint construction | `geometry/core.py`（B 对角切向、W 平均切向、edge midpoint） | 计算 checkerboard B 面与 W 点的局部 frame |
| metric / curvature | `geometry/core.py`、`geometry/boundary.py`（第一/第二基本形式和边界闭合） | 得到 `a_B`、`a_W`、对称 `b_W` 和 `chi` |
| `X_M` | `geometry/intrinsic.py`（B-frame connection 与 M-character moment） | 论文膜能量的内禀补全通道 |
| `X_W^rel` | `geometry/boundary.py`（参考 frame 的 polar transport） | 论文弯曲兼容通道 |
| material / energy | `material/isotropic.py`、`mechanics/energy.py`（StVK/Koiter 分块） | 组合 `membrane_b + membrane_xm + bending + bending_compatibility` |
| residual / tangent | `mechanics/model.py`（`residual()` / `tangent()`） | 由总势能自动微分，保持能量-残量一致 |
| Newton | `solver/newton.py`（论文 Newton 线搜索实现的最小复制） | Hessian 特征分解下降方向、Armijo/残量线搜索 |
| visualization | `visualization/plotting.py`（可选绘图） | 参考/变形 quad 网格三维图 |

## 公开入口

`CheckerboardShell(mesh, material, boundary_conditions, load).solve()`（参数化正向静力求解入口）返回 `SimulationResult`（变形、位移、能量和收敛记录）。几何、材料、约束、载荷彼此独立，可在未来由 design parameters 外层驱动。

