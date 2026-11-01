# Final exact-tangent and nonlinear-design reliability audit

Date: 2026-10-30. This is the closing audit for the two existing engineering
specimens. The original strict comparison remains in
`efficiency_report.md`; results below supersede its automotive endpoint and
branch conclusions. All timed jobs used one OpenBLAS/OpenMP/MKL thread on the
same Intel Xeon Gold 6148 host. Aircraft values are three fresh-process
medians. Automotive complete optimizations are one full run per method;
their workflow ratios are descriptive, without a repeat-run uncertainty
estimate.

## A. Why the exact tangent is expensive

The accepted CB-Full energy has a local shell Hessian plus the exact
`2c DX_M DX_Mᵀ` rank-one correction. The production local Hessian uses
147 graph colours at radius 3. The original batch size 7 performs 21
batched reverse-mode Hessian-vector calls. A 20×40 cProfile run measured
3.059 s for tangent construction, of which 2.906 s was the coloured
Hessian-vector phase; graph creation took 0.017 s, residual evaluation
0.048 s, and rank vector evaluation 0.015 s
(`results/hessian_profile_20x40_legacy.json` and its `.txt` profile).
At the old 20×40 aircraft screen, sparse factorization took about
0.17 s and all candidate right-hand-side solves about 0.08 s.
Thus the dominant cost is repeated PyTorch reverse-mode execution, not
the 3V sparse linear algebra. A `torch.func` forward-over-reverse trial
was slower for a representative 21-direction batch, so it was not adopted.

## B–C. Exact implementation change and numerical equivalence

`mechanics/low_rank.py` now caches only grid topology: colour seeds and
COO row/column maps. The default uses 21 directions per batched reverse
pass, with the original 7-direction assembly still callable as
`assembly_strategy="legacy"`. Configuration-dependent geometry, strains,
energy, residual and Hessian values are always recomputed. The sparse
local-plus-rank-one definition is unchanged. The load-factor residual is
evaluated in `model.low_rank_tangent`; an algebraically identical correction
formerly applied by its Newton caller was removed. A non-unit-load test
checks the reduced sparse residual against independent AD residuals.

The old and cached tangents agree at reference, sinusoidal deformations
up to 20 mm in the unit test, and a 90° rigidly rotated/deformed configuration.
Energy objectivity, force, local matrix, rank vector, Hessian-vector action
and symmetry are covered by the regression suite. In the 31-layout aircraft
audit, **every layout and both targets have exactly equal saved strokes,
reactions, errors and objective values** for old versus new assembly
(`results/aircraft_efficiency_equivalence.json`). At one saved automotive
input, old/new cold solves have identical displacement to stored precision
and the same six Newton iterations
(`results/automotive_tangent_before_after.json`).
The complete old/new automotive consistent-gradient runs have identical
35-iteration paths and final controls.

## D. Aircraft: exact shared 31-layout screen

| Mesh | CB original exact tangent | CB cached exact tangent | MITC4 shared | CB original/cached tangent assembly |
|---|---:|---:|---:|---:|
| 12×24 | 1.599 s | 1.688 s | **0.497 s** | 1.390 / 1.481 s |
| 16×32 | 2.384 s | 2.310 s | **0.873 s** | 1.989 / 1.935 s |
| 20×40 | 3.876 s | 3.890 s | **2.033 s** | 3.133 / 3.127 s |

All entries include case construction, exact tangent assembly, shared
candidate responses and all 31 layout optimizations for two target
conditions. These are three-run medians; raw repeats, RSS, state DOFs and
factor/response times are in `results/tangent_aircraft_final.csv`.
The optimized CB implementation is 5.6% slower, 3.1% faster and 0.4%
slower than the original at the three meshes. This is **not a stable
aircraft speed improvement**. MITC4 remains 3.4×, 2.6× and 1.9× faster
than the cached CB implementation. At 20×40, CB has 2583 state DOFs
versus MITC4's 5166, yet roughly 3.13 s of its 3.89 s screen is spent
building the tangent. The 3V state and sparse-plus-rank-one solver reduce
matrix work but do not outweigh this element/AD cost. Different
reinforcement rotation kinematics and OpenSees/PyTorch implementations
also prevent a pure shell-element speed attribution. The previously
validated `[0,2,4]` nonlinear design is unchanged; the timing result
applies to **linear layout screening only**.

## E. Automotive gauge, gradients and KKT reliability

The old automotive panel fixed two datum nodes on the same checkerboard
sublattice. An exact relative-sublattice x-translation remained free.
At the ASD-designed controls, two old cold/ASD-initial states differed
by a 0.956 mm x-shift on precisely the opposite sublattice, while their
y/z fields agreed below `1e-15 m` and total internal-plus-bracket
potentials differed by only `1.1e-15 J`. This was a **zero-energy gauge**,
not evidence of two distinct stable equilibrium branches. A single
additional x datum on the other sublattice removes it. Its reaction is
about `8e-9 N` at the final CB design. The same locator is now imposed
on the ASDShellQ4 optimization and common reference: at one saved input,
it carries about 82 N and changes displacement by up to 0.10 mm, so a
comparison omitting it would not have matched boundary conditions.
From both tested initial
states, the gauged model reaches the same ASD-control edge response
(`0.57369846625` versus `0.57369846623 mm`); the minimum reduced
tangent eigenvalue is about 90 N/m at that probe
(`results/automotive_efficiency_branch_energy.json`). This is a
necessary positioning constraint for the position-only model, not a
new bracket or a change to the shell energy. It shifts the nominal
edge response and therefore requires the full reoptimization below.
The locator is a shared physical boundary condition in the final
comparison, even though its CB reaction at the CB optimum is negligible.

For the new CB consistent endpoint, independent cold-state central
differences at 0.05 and 0.01 mm give maximum objective-gradient errors
`1.53e-5` and `6.11e-7` objective units/mm. The maximum shape-Jacobian
errors decline from `4.41e-5` to `1.77e-6` mm/mm; force-Jacobian errors
from `3.28e-4` to `1.30e-5` N/mm. This supports the implicit
derivative, including the bracket reaction path
(`results/automotive_reliability_cb_consistent_zero.json`).
The same 0.01 mm check at the FD endpoint gives `6.12e-7`
objective-gradient discrepancy.

Both final CB designs activate the third variable's lower -1.5 mm stroke
bound and seventh variable's upper +1.5 mm bound. All nonlinear
shape/reaction inequalities have at least 0.438 mm margin in their
native mixed constraint vector (force margins are larger). At the
consistent endpoint, the maximum free/projection KKT gradient is
`4.73e-7` objective units/mm; the active-bound gradient signs are
correct, the reduced equilibrium tangent minimum eigenvalue is
`92.47 N/m`, and the equilibrium residual is below the solver tolerance.
At the FD endpoint, the projected gradient is `1.48e-5`, consistent
with its 0.05 mm FD derivative accuracy. These are verified local
stationarity results for the deterministic cold path, not proofs of
global optimality or of agreement between different shell formulations.

## F. Automotive complete optimization and common reference

| Method, 12×8 design grid | Complete time | SLSQP iterations | Equilibria | Newton steps | Design objective | 24×20 ASD reference objective | Reference max edge / bracket |
|---|---:|---:|---:|---:|---:|---:|---:|
| CB consistent, cached tangent | **200.3 s** | 35 | 50 | 317 | 0.031260202 | 0.037786 | 0.2569 mm / 88.55 N |
| CB centered FD, cached tangent | 2509.4 s | 35 | 766 | 4420 | 0.031260204 | 0.037785 | 0.2568 mm / 88.55 N |
| ASDShellQ4 corotational centered FD | 3676.5 s | 41 | 708 | 93541* | **0.015068** | **0.022173** | 0.3527 mm / 67.05 N |

All three controls pass the **same 24×20 corotational ASDShellQ4**
research limits of 0.6 mm edge error, 300 N bracket resultant and
±1.5 mm stroke. Reference acceptance solves take 141.6, 143.3 and
184.8 s respectively, separate from optimization time. Including one
reference acceptance check gives 341.9, 2652.7 and 3861.3 s.
The CB-FD and
CB-consistent designs are nearly identical in both design and common
reference objectives; CB-consistent is **12.53× faster** as a complete
same-discretization optimization. CB-FD is **1.47× faster** than the ASD-FD
workflow to an accepted engineering design, but the ASD optimum has a
materially lower reference objective. CB-consistent versus ASD-FD is
18.35× by time to accepted design, mixing both solver and gradient
strategies. Neither cross-shell ratio is an equal-objective optimum
speedup. With the common locator, the ASD design gives 0.581 mm maximum
error in the separately refined CB model, inside but closer to the
0.6 mm limit than its own reference response. This remaining difference
does not establish global interchangeability between formulations.

The independent ASD cold-path centered-difference gradient check gives
projected infinity norms `2.10e-5` and `3.80e-5` objective units/mm
at 0.05 and 0.01 mm steps respectively. Its bounds and nonlinear
inequalities are inactive. The non-monotone difference indicates a
noise floor in these nonlinear FD probes; the SLSQP success flag and
engineering acceptance are verified, but ASD stationarity is **not
independently certified to the CB gradient tolerance**. This further
limits an equal-optimum timing interpretation
(`results/automotive_asd_stationarity.json`).

The one-to-one old/new CB consistent-gradient comparison is more causal:
with identical 35 iterations, 50 equilibria, 317 Newton steps and final
controls, the legacy assembly takes **266.3 s** versus **200.3 s**
cached (**1.33×** total speedup). Tangent assembly falls from
220.3 to 160.3 s; the single-input cold-forward three-run medians fall
from 4.549 to 3.446 s (**1.32×**) with identical six-step path,
and factorization remains about 0.0145 s. The 20×40 aircraft
single-tangent performance does not show this gain, so the improvement
is mainly valuable where small/medium tangents are repeatedly rebuilt.

At the **same final CB controls**, three cold forward runs have median
19.20 s for CB and 4.75 s for matched-locator ASDShellQ4
(`results/automotive_efficiency_forward.json`). CB takes 33 Newton
steps over direct/continuation attempts and 37 tangent assemblies,
versus 122 reported OpenSees convergence iterations over two load
increments. Its tangent assembly alone takes about 16.8 s. Thus
the smaller CB state is **not a faster nonlinear forward solve at this
particular design**. The 200 s versus 3676 s whole-optimization
difference is dominated by the 50 versus 708 equilibrium analyses
and by the different gradient strategies and nonlinear paths. The
4.549/3.446 s old/new matched-path tangent probe above uses a different
saved input; its causality check must not be substituted for this
final-design forward cost.

\*The installed OpenSees `testIter()` total sums its nonlinear convergence
iterations over load increments. OpenSees does not expose separately
comparable numeric assembly/factorization counts here. The ASD optimizer
**was rerun from zero** with the same extra x locator, identical design
objective/constraints, 2-step continuation and 0.05 mm centered FD.
Its previous ungauged 3559.9 s result is retained separately as
`results/automotive_efficiency_asd_fd_12x8_zero_pre_common_locator.json`
and is excluded from the final fair ratio.

## G–H. Attribution and manuscript conclusion

| Source | Measured contribution | Interpretation |
|---|---|---|
| Position-only 3V state, sparse plus rank-one tangent | Half the shell state DOFs; 20×40 CB factorization ~0.17 s, versus ~3.13 s tangent assembly | Compact algebra is real, but **no aircraft total-time advantage** follows |
| Cached colouring and larger AD batches | Automotive matched-path 266.3→200.3 s; fixed forward 4.549→3.446 s; aircraft no stable gain | Implementation improvement, not new mechanics or a universal speedup |
| Exact implicit gradient | Matched CB design objectives with 50 versus 766 equilibria and 200 versus 2509 s | Strong differentiable-framework benefit; implicit differentiation is general |
| Shared layout response | Both CB and MITC4 use it; MITC4 faster at all tested aircraft meshes | Generic reuse, not a Position-Only-specific benefit |
| Conventional nonlinear shell workflow | ASD-FD 3676.5 s with matched locator and lower reference objective; at final CB controls ASD forward 4.75 s versus CB 19.20 s | Different element, code and nonlinear iteration path; no isolated element speed ratio |

The evidence now supports a precise CMAME statement: the position-only
formulation has a compact exact tangent and provides a practical
differentiable basis for accepted inverse designs; the automotive
consistent-gradient framework is substantially faster than finite
differences for the **same validated local optimum**. The present
PyTorch Hessian assembly prevents a general element-level speed
claim, and the shared MITC4 aircraft screen is faster. The automotive
cross-shell endpoint objectives differ, so its 1.47×/18.35× ratios
describe time to the stated engineering tolerances only. No remaining
identified error invalidates the two final CB automotive controls under
the common reference. The exact stationarity of the ASD-FD endpoint
remains unresolved at the tighter CB tolerance. More realistic hardware, connector and
production CAD studies are optional engineering extensions, not required
for this method-level conclusion.

## Reproduction and files

Run from the repository root in an environment with the project,
OpenSeesPy and pytest installed. The commands below reproduce the
measurements without a machine-specific Python path:

```bash
python -m pip install -e . openseespy pytest
export PYTHONPATH=src:.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m experiments.engineering_inverse.profile_hessian --strategy legacy
python -m experiments.engineering_inverse.profile_hessian --strategy cached
python -m experiments.engineering_inverse.efficiency_aircraft --methods cb_shared_sparse_legacy cb_shared_sparse mitc4_shared --repeats 3
python -m experiments.engineering_inverse.validate_aircraft_efficiency
python -m experiments.engineering_inverse.efficiency_automotive cb_consistent --start zero
python -m experiments.engineering_inverse.efficiency_automotive cb_consistent --start zero --hessian-strategy legacy
python -m experiments.engineering_inverse.efficiency_automotive cb_fd --start zero
python -m experiments.engineering_inverse.efficiency_automotive asd_fd --start zero
python -m experiments.engineering_inverse.efficiency_automotive_reference --start zero
python -m experiments.engineering_inverse.efficiency_forward
python -m experiments.engineering_inverse.validate_automotive_reliability --method cb_consistent
python -m experiments.engineering_inverse.validate_automotive_reliability --method cb_fd --steps 0.01
python -m experiments.engineering_inverse.validate_automotive_asd_stationarity
python -m experiments.engineering_inverse.validate_automotive_locator_effect
python -m experiments.engineering_inverse.automotive_branch_probe
python -m experiments.engineering_inverse.automotive_branch_energy
python -m experiments.engineering_inverse.make_tangent_reliability_figures
python -m pytest -q
```

The full regression suite reports **33 passed**. Final raw JSON, detailed
profile text and flat CSV are in `results/`; eight final PDF/SVG figure
pairs are in `figures/`. The relevant code changes are
`src/checkerboard_shell/mechanics/low_rank.py` (exact cached colouring),
`mechanics/model.py` (strategy selection and algebraically equivalent
load-factor handling), `automotive_panel.py` and `mount_springs.py`
(opposite-sublattice datum and reaction accounting),
`automotive_spring_design.py` (sensitivity with the same datum),
`efficiency_automotive.py` and `efficiency_aircraft.py` (matched
old/new strategy runs), plus profile, reliability, plotting and
regression scripts. The shell energy, Blind-Mode Completion, material,
brackets and engineering objectives/limits were not changed.
