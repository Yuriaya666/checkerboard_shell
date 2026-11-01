# Same-method: adjoint vs finite difference

Same Checkerboard model, same target, same L-BFGS settings, same linear
solver.  Model/tangent assembly is timed once and included in `total`;
`optimization` is the design loop only.

| method | iterations | linear solves | optimization (s) | total (s) | target error |
|---|---:|---:|---:|---:|---:|
| adjoint | 11 | 68 | 0.18 | 1.08 | 0.0067% |
| finite_difference | 10 | 391 | 1.01 | 1.91 | 0.0067% |

Model/tangent assembly: 0.89 s.
Linear-solve ratio: 5.8x; optimization-time ratio: 5.6x.
