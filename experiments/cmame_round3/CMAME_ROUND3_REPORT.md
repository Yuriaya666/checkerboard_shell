# CMAME Round 3 report

## Q1. Does the evidence show that beta_M=8mu is not benchmark-specific tuning?

Yes, for the tested manuscript responses, with a deliberately limited scope.
`beta_M=8mu` is first fixed by the trace-free-shear normalization of the
chosen `X_M`, independently of a load or benchmark.  The recovered-driver
sweep then shows that `c_beta=1` is not a sharp displacement optimum:

| recovered case | maximum response change vs `c_beta=1`, `c_beta` in 0.25--4 | maximum over 0--8 | `E_XM/E_int` at 1 |
|---|---:|---:|---:|
| plate n=16 | 0 | 0 | 0 |
| Scordelis--Lo n=16 | 1.63e-5 | 1.36e-3 | 3.15e-7 |
| M3 n=12 | 4.21e-10 | 5.04e-10 | 4.20e-24 |

The plate result is exact because the recovered manuscript driver is the
scalar transverse branch, in which the membrane-completion channel is not
active.  For roof and M3, completion-added `DX_M` stiffness increases exactly
linearly with `c_beta`, while normal load-bearing displacement remains on a
broad plateau.  Thus the data support “not benchmark-fitted” and “low normal-
response sensitivity on these cases,” not “beta_M is unimportant.”

Removing the completion does not create a new full-tangent nullity for the
curved roof/M3 problems under their benchmark BC.  In the same plate n=16
flat in-plane membrane operator, however, `c_beta=0` increases nullity from 5
to 6; the isolated extra direction has base stiffness `9.62e-20` and receives
`8.27e-5` completion stiffness at `c_beta=1`.  This is the correct restricted
evidence for recovery of the one-dimensional flat quotient defect.

## Q2. What are now the strongest comparative advantages?

1. **Position-only state economy with competitive accuracy per state.** Across
   the 12 recovered plate/roof/M3 mesh pairs, Checkerboard uses 50.5--55.1% of
   the MidedgeTan raw scalar states while producing convergent continuum-
   reference errors of comparable order.  This is a state-count result, not a
   CPU or memory superiority claim.
2. **Exact physical quotient and finite-configuration gauge invariance.** Raw
   representatives can differ strongly while their midpoint surfaces overlap;
   after the production optimization the maximum relative tangent change is
   `1.25e-15` for `alpha/L<=1`.
3. **Exact planar local geometry.** Every nondegenerate Varignon B face has an
   unambiguous plane, frame, and normal for arbitrary raw positions.  This is
   geometric well-posedness, not a promise of distortion-insensitive error.
4. **Sparse algebraic structure without hidden dense X_M coupling.** The exact
   nonlinear tangent retains `2c X_M D2X_M` locally and represents the sole
   global term as one rank-one border.  Measured augmented nnz grows nearly
   linearly through n=128.

The current data do **not** establish cross-backend runtime or memory
superiority, general distorted-mesh robustness, MITC4 superiority, arbitrary
variable-spacing validity, or locking-free convergence.  In particular, the
M3 surface-preserving redistribution sensitivity remains material: at
delta=.30 the Checkerboard errors are 20.12% (smooth) and 28.89% (random).

## Q3. What is the strongest CMAME claim supported by the current data?

> Checkerboard midpoint geometry yields a structured, position-only shell
> discretization whose physical local patches are exactly planar and whose raw
> checkerboard gauge is preserved to numerical roundoff at finite
> configurations.  On the recovered plate, Scordelis--Lo, and M3 benchmark
> sequences, it achieves competitive continuum-reference accuracy with about
> half the raw scalar states of the verified MidedgeTan baseline.  The unique
> one-dimensional flat membrane quotient defect is removed by a
> quotient-minimal rank-one observation—minimal only in the linear
> observation-rank sense—and the exact nonlinear production tangent is a
> bounded-stencil sparse local operator plus one global rank-one update.  These
> results establish state-count economy and sparse algebraic structure, while
> not implying CPU superiority, distortion-insensitive approximation, or
> locking-free behavior.

This wording preserves the strongest validated contributions and explicitly
separates geometric well-posedness from approximation robustness, and state-
count/algebraic economy from measured runtime.
