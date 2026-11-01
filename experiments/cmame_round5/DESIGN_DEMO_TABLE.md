# Design-demo summary (Round 5)

| # | demo | design variables | free DOFs | key result | iterations | figure |
|---:|---|---:|---:|---|---:|---|
| 1 | adjoint gradient check | 9 | 297 | max relative error vs FD `1.5e-12` | — | — |
| 2 | actuation shape control | 9 nodal forces | 297 | target displacement error `6.5e-5`; force error `0.174` (non-unique) | 13 | `Figure_R5_actuation_design.png` |
| 3 | reference-shape optimization | 12 height controls | 297 | target displacement error `1.6e-4`; height field non-unique (diagnostic) | 20 | `Figure_R5_shape_design.png` |
| 4 | large-rotation shape design | 2 initial heights | 46 | tip error `1.9e-6`; `h_opt ≈ h_true` to `1e-6` | 8 | `Figure_R5_large_rotation_design.png` |

Demos 1--3 are linear reference-tangent problems; demo 4 is nonlinear
equilibrium at `M/Mmax=0.25`.  The 3V adjoint/design variables are the point:
no rotational design variables are used.  Actuation forces and reference
heights are not unique from a displacement target; the primary metric is target
displacement (demo 3) or target curvature/tip (demo 4).

## Same-target design-loop cost (9 actuators)

| method | design vars | free DOF | iterations | work | total time | target error |
|---|---:|---:|---:|---:|---:|---:|
| Checkerboard + adjoint | 9 | 297 | 13 | 79 linear solves | **1.05 s** | **0.0065%** |
| Checkerboard + finite difference | 9 | 297 | 12 | 438 linear solves | 2.08 s | 0.0066% |
| MITC4 + finite difference | 9 | 594 | 17 | 533 forward analyses | 11.85 s | 0.0026% |

Three layers: (1) same-method adjoint vs FD: 79 vs 438 solves, 5.5x fewer,
5.5x faster optimization; (2) state size: 297 vs 594 free DOF, about 1/2;
(3) black-box MITC4 FD: 11.85 s vs 1.05 s, 11.3x total.  The missing
same-gradient-method row is Checkerboard FD vs MITC4 FD: 2.08 vs 11.85 s,
5.7x, in FD_VS_FD.md.  The MITC4 comparison is
implementation-level: it rebuilds the model per analysis and uses FD; a MITC4
adjoint would reduce its work units.
