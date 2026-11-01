# Engineering inverse design: strict computational efficiency study

Date: 2026-10-30. Timings were obtained on one Intel Xeon Gold 6148 CPU with
`OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1`. Each aircraft repeat
runs in a fresh process; preparation, construction, response generation and
all 31 layout optimizations are included. Raw records are in `results/`.
This study changes no shell energy, Blind-Mode Completion, material, loads,
engineering objectives or design constraints. The two specimens remain the
representative research models documented in `report.md`.

## §6.3 Aircraft trailing-edge skin: 31-layout linear screening

Five candidate stations give 31 nonempty layouts and two independent target
conditions. The timed task is **linear screening**, not a complete nonlinear
layout optimization. The previously verified three-station `[0,2,4]` solution
has separate nonlinear Checkerboard and corotational shell/frame checks in
`results/aircraft_nonlinear_final_12x24.json` and
`results/aircraft_corotational_reference.json`; no minimum-cardinality claim
is made.

The original Checkerboard direct path factors the layout-specific reduced
matrix for each of 31 layouts. Its shared path factors the spar-clamped base
once and computes candidate compliance columns. The new ShellMITC4 shared
path builds one spar-clamped OpenSees model, prescribes pressure and unit
candidate loads in successive time-series cases, and requests
[`Linear -factorOnce`](https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/algorithm/LinearAlgorithm.html).
Selected compliance submatrices then enforce the same patch strokes and
return the same resultant reactions as direct prescribed-DOF elimination.
The installed OpenSeesPy `SparseGeneral` system reports `numFact()==0`, so
the exact numeric factorization count is **not independently observable**;
the option's documented behavior and the response timings support the
intended reuse. Sparse numeric `K` is also not exposed through this interface:
the MITC4 nonzero counts below are connectivity-graph counts, not measured
numeric nonzeros.

The fifth method retains the exact Checkerboard local sparse plus rank-one
tangent in a bordered SuperLU system. It removes the former dense outer
product in this linear screening code. It is a measured implementation
improvement, separate from the shared-compliance algorithm. All 31 layouts
and both targets were compared numerically between direct, shared and sparse
paths; the largest within-formulation discrepancy is under `5.6e-5 mm` in
maximum error and `5.2e-5 N` in patch reaction
(`results/aircraft_efficiency_equivalence.json`).

| Mesh | CB direct | CB shared, dense | CB shared, sparse + rank one | MITC4 direct | MITC4 shared |
|---|---:|---:|---:|---:|---:|
| 12×24 | 4.08 s | 1.67 s | 1.56 s | 9.63 s | **0.505 s** |
| 16×32 | 12.03 s | 2.59 s | 2.31 s | 19.23 s | **0.908 s** |
| 20×40 | 42.22 s | 4.88 s | 4.00 s | 32.94 s | **2.04 s** |

Values are three-run medians. The old comparison of CB shared against MITC4
direct suggested a 5.8–7.4× advantage in these reruns, but it mixed different reuse
strategies. Under shared operators, MITC4 is **3.09×, 2.54× and 1.96×
faster** than the improved CB path at the three meshes. The ratio declines
with mesh size over the tested range; extrapolation beyond these meshes is
unsupported. The generic reuse strategy accelerates CB direct by 2.6–10.6×
and MITC4 direct by 16.2–21.2×, depending on mesh. The two direct paths do
different numbers of model rebuilds (32 CB factorizations versus 111 legacy
OpenSees analyses); their ratio is a workflow comparison only.

| Mesh | CB state / free DOF | MITC4 state / free DOF | CB local numeric nnz | CB bordered nnz | MITC4 graph nnz | CB sparse / MITC4 shared peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| 12×24 | 975 / 825 | 1950 / 1800 | 69,663 | 58,036 | 86,724 | 563 / 532 MiB |
| 16×32 | 1683 / 1485 | 3366 / 3168 | 121,389 | 106,590 | 157,140 | 598 / 541 MiB |
| 20×40 | 2583 / 2337 | 5166 / 4920 | 187,967 | 170,656 | 248,292 | 621 / 557 MiB |

The fewer CB state variables and stored sparse entries do **not** translate
to lower whole-process peak RSS in this Python/PyTorch implementation.
PyTorch, OpenSees and runtime allocations are included; peak RSS is not a
standalone element storage metric. At 20×40, sparse/rank-one storage reduces
CB peak RSS from about 824 MiB (dense shared path) to 621 MiB and total time
from 4.88 to 4.00 s. The change passes the same response/reaction regression
test and does not alter the tangent definition.

At 20×40, the sparse CB base shell tangent assembly takes 3.19 s of its
4.00 s screen; its bordered factorization takes 0.169 s and all 84 candidate
RHS solves 0.077 s. The MITC4 shared model construction takes 0.033 s,
the first pressure analysis (assembly, factorization and solve combined)
0.296 s, and 84 candidate load cases 1.49 s. These stages show that CB has
a cheap sparse algebra phase once its tangent exists, while PyTorch automatic
differentiation and colored Hessian assembly dominate the complete screen.
The OpenSees first analysis does not expose separate element assembly and
numeric factorization times; the candidate cases include OpenSees load-step
and result-extraction overhead. Thus the difference in algebra phase cannot
be attributed solely to the 3V/6V DOF ratio. The coupled reinforcements also
use position-line curvature in CB and rotational frame DOFs in OpenSees;
their kinematics are close but not identical (`results/beam_kinematics_audit.json`).

### Engineering outcome held separate from screening time

The saved `[0,2,4]` nonlinear CB design meets both targets (12×24 maximum
shape errors 0.412 and 0.599 mm; tolerance 0.8 mm; RMS tolerance 0.25 mm).
The 20×40 independent corotational shell/frame check gives RMS 0.181 and
0.206 mm and maxima 0.412 and 0.644 mm. Its material, geometry, loads and
attachment locations match, while the beam/shell rotation coupling remains
an explicitly measured model difference. The present timing table neither
includes these nonlinear checks nor proves that every linear-screen layout
would pass a nonlinear reference.

## §6.4 Automotive panel: complete nonlinear optimization

The comparison uses the existing 0.80×0.45 m representative panel, 0.7 mm
steel, four finite-stiffness bracket regions, eight adjustment variables,
the same weighted edge objective and SLSQP with `ftol=1e-10`, up to 80
iterations, ±1.5 mm strokes, 300 N bracket resultants and ±0.6 mm edge
component constraints. The common start is zero adjustment. For both finite
difference runs, a centered 0.05 mm step and one shared perturbed state per
sign/direction supplies objective, shape and force derivatives. This avoids
penalizing either solver for duplicate optimizer callbacks. The independent
reference is ASDShellQ4 with corotational kinematics at 24×20, with 16 load
steps for the acceptance check. Diagnostic runs showed repeated failure of
one-step full-load trials at nonzero FD probes while two-step continuation
converged; the timed nonlinear reference optimization therefore starts with
2 steps and retries with 4/8/16 steps if needed. OpenSees' documented
`NormUnbalance` default is the global 2-norm, whereas CB checks the largest
free residual component. The timed ASDShellQ4 path explicitly selects the
[max-norm option](https://opensees.github.io/OpenSeesDocumentation/user/manual/analysis/test/NormUnbalance.html)
at `1e-6` with up to 80 iterations per increment; older reference scripts
retain their original default for reproducibility. At the saved 12×8 design,
changing the OpenSees norm from 2 to max changed the converged displacement
by at most `8.8e-11 m` and bracket reaction by `7.7e-8 N` while reducing two
step Newton iterations from 133 to 124
(`results/automotive_efficiency_norm_equivalence.json`). A
separate saved-design check compares the converged one- and 16-step responses.

The complete optimizer and common-reference records appear in
`results/automotive_efficiency_*.json` and
`results/automotive_efficiency_common_reference_zero.json`. Every new design
input starts from the reference configuration in all three methods. Earlier
warm-start CB runs were retained only as diagnostics because using the last
optimizer trial as the next equilibrium initial condition made the objective
path dependent in this nonlinear model (`results/automotive_efficiency_cb_*_warm_diagnostic.json`).

| Method, 12×8 design grid | Total time | SLSQP iterations | Equilibria | Newton steps | Gradient evaluations | Design-grid objective | Common 24×20 objective | Common max error / force |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| CB consistent derivative | **238 s** | 19 | 52 | 296 | 19 | 0.03254 | 0.03784 | 0.198 mm / 73.9 N |
| CB centered FD | 1,384 s | 20 | 415 | 1,773 | 20 | 0.03564 | 0.05007 | 0.174 mm / 19.3 N |
| ASDShellQ4 corotational centered FD | 3,560 s | 40 | 680 | 90,752* | 40 | 0.01439 | **0.02803** | 0.330 mm / 78.2 N |

All three optimizers report success and all three **independently pass** the
24×20 ASDShellQ4 shape limit 0.6 mm, bracket resultant limit 300 N and
stroke bounds ±1.5 mm. The objective values are weighted squared edge-error
values plus the identical stroke penalty, not physical strain energies.
The optimization times above exclude the separate high-resolution acceptance
solve, which takes 143/152/178 s for CB-consistent/CB-FD/ASD-FD controls.
Including exactly one common-reference acceptance solve gives **381/1,536/
3,738 s** for the three respective complete design-and-check workflows.
The separate high-resolution CB cross-check is diagnostic and is excluded
from those totals. Reference solve time varies with the design input and is
reported as a single run, whereas the aircraft timings are three-run medians.
The reference objectives differ materially: the ASD design is best by this
measure and is not an equal-quality endpoint to either CB design. The
comparison of **time to a reference-accepted engineering design** is
therefore 2.57× for CB-FD versus ASD-FD and 14.95× for CB-consistent versus
ASD-FD, but neither is a rigorously controlled time-to-equal-objective
speedup. The former uses the same gradient strategy but different nonlinear
discretizations and optimizer paths; the latter adds a gradient-method
benefit. Within CB, the consistent run uses 5.81× less total time than CB-FD
and reaches a better common-reference objective. This is an observed
workflow result for this start and tolerance, not a general algorithmic
guarantee.

\*`OpenSees testIter()` summed over all load increments reports 90,752. This
is far larger than CB's 1,773 FD Newton steps and is a major difference in
the nonlinear solver paths; it cannot be assigned to shell state dimension
alone. OpenSees does not expose its sparse tangent assembly and numeric
factorization counts through the installed interface, so those entries are
unavailable rather than estimated.

| Same 12×8 saved CB input, cold solve | CB | ASDShellQ4 |
|---|---:|---:|
| State DOF / free shell DOF | 351 / 345 | 702 / 696 |
| Numeric local shell nnz / structural shell graph nnz | 23,823 + rank-one vector of 346 entries | 33,030 graph entries; numeric nnz unavailable |
| Median complete forward time, three repeats | **3.93 s** | 4.93 s |
| Tangent assemblies / reported Newton iterations | 6 / 5 | unavailable / 124 over two load steps |
| CB tangent assembly / numeric factorization, median | 3.58 s / 0.012 s | unavailable separately |
| Peak process RSS in this forward experiment | 550 MiB | 571 MiB |

Both cold forward runs use the same saved physical input and pass their own
equilibrium tests; however, ASDShellQ4 uses two load increments while CB's
direct solve uses one full-load Newton path. The median forward ratio is only
1.25× in CB's favor. The CB factorization is cheap, but tangent AD assembly
occupies most of its forward time. OpenSees' model construction averages only
about 0.005 s after first use, whereas its nonlinear analysis accounts for
nearly all of its forward time. Peak RSS is whole-process, including Python,
PyTorch and OpenSees, and should not be interpreted as element storage.

### Gradient and solver cost attribution

The complete cold-start CB consistent run spends 226 s in forward states and
12.0 s in 19 implicit derivative computations (mean 0.632 s each). Its 52
equilibria include continuation retries; 348 tangent assemblies take 205 s,
while 296 numeric factorizations take 0.715 s. CB-FD spends 1,110 s inside
20 centered FD gradient computations (mean 55.5 s each) and performs 415
equilibria; 2,188 tangent assemblies take 1,257 s. ASD-FD spends 3,353 s in
40 centered gradients (mean 83.8 s each) and performs 680 equilibria. The
ASD model construction totals only 3.44 s over those analyses; its nonlinear
analysis time is 3,550 s. These totals are nested: FD gradient time is part
of forward time, and tangent/factorization times are part of nonlinear
analysis time. They must not be added as disjoint bars.

The same-input one-step versus 16-step ASDShellQ4 check differs by at most
`4.3e-11 m` in displacement and `6.3e-8 N` in bracket resultant
(`results/automotive_efficiency_step_equivalence.json`), although many
optimizer probes require two steps to converge. A 12×8 saved-input tolerance
check finds that loosening the reference max-norm threshold from `1e-6` to
`1e-5` reduces two-step analysis time from 4.90 to 4.27 s while changing
displacement by at most `7.5e-10 m`; the published optimization retains
`1e-6` (`results/automotive_efficiency_tolerance_sensitivity.json`).

### Important nonlinear qualification

Both cold-start CB runs finish with SLSQP's success flag, but a fresh
consistent-gradient evaluation at their endpoints gives infinity norms
`0.00339` and `0.0200` in the objective units per millimetre of stroke.
Their strokes remain at least 0.246 and 1.138 mm from the ±1.5 mm bounds,
and their shape and force constraints are inactive. These nonzero gradients
mean the endpoints are **not independently verified stationary points**.
Small design changes and different initial displacement fields can select
different converged equilibrium branches. In one probe at the ASD optimum,
12×8 CB converges from a cold state to 0.588 mm maximum edge error, and from
the ASD displacement as initial state to 0.715 mm, with residuals below
`2.2e-8 N` in both (`results/automotive_efficiency_branch_probe.json`). At
24×20, the ASD-designed controls give 0.330 mm in
the common ASDShellQ4 reference but 0.648 mm in the independently refined
CB solve, which exceeds the 0.6 mm shape limit. This is a material
model/discretization and branch-selection limitation; it prevents claiming
the two nonlinear shells give interchangeable design sensitivities over the
whole feasible region. The accepted common-reference designs remain valid
for the stated 24×20 ASDShellQ4 criterion.

## Efficiency attribution and publishable claims

| Factor | Aircraft evidence | Automotive evidence | Supported inference |
|---|---|---|---|
| Position-only 3V state | Half the total DOF; lower CB sparse stored entries, but slower matched shared screen | 351 versus 702 state DOF; same-input forward 3.93 versus 4.93 s | Compact state and cheap CB factorization; **no intrinsic whole-workflow speedup proved** |
| Local sparse plus rank-one tangent | Sparse CB screen 4.00 versus dense CB 4.88 s at 20×40; RSS 621 versus 824 MiB | CB cold forward factorization ~0.012 s versus tangent assembly ~3.58 s | Exact storage/solver improvement within CB; AD assembly dominates total |
| Shared candidate compliance | CB direct/shared-sparse 42.22/4.00 s; MITC4 direct/shared 32.94/2.04 s at 20×40 | Not the repeated-layout problem | Generic reuse benefits **both** discretizations; cannot be credited to CB alone |
| Consistent design derivative | Not used in linear screening | CB implicit 0.632 s/evaluation, centered FD 55.5 s/evaluation; complete runs 238/1,384 s with differing endpoints | Strong measured benefit of this CB *implementation*; implicit differentiation is not unique to CB |
| Whole-task engineering acceptance | `[0,2,4]` nonlinear design separately passes reference; timing measures only linear screen | All three designs pass one 24×20 reference; common objective differs | CB is faster for reference-accepted automotive workflows here; not an equal-objective or element-only speedup |

The defensible CMAME statement is that this Position-Only implementation
offers a compact state and an exact sparse-plus-low-rank tangent, and that
its differentiable inverse-design framework can produce reference-accepted
designs at lower measured workflow cost in the automotive example. The
aircraft shared-operator experiment is a **negative efficiency result** for
the current implementation: the optimized OpenSees MITC4 screen is faster
at all three tested meshes despite twice the state DOF. Thus no general
claim that Position-Only Checkerboard is faster than conventional shell
elements is supported. The automotive workflow advantage combines shell
representation, very different nonlinear solver iteration counts, Python/
PyTorch versus OpenSees implementation, and gradient strategy; these cannot
be disentangled into a unique percentage assigned to the shell formulation.

## Reproduction and outputs

The efficiency-specific code changes are `mitc4_shared.py` (shared OpenSees
responses), `aircraft_sparse.py` (exact sparse plus rank-one CB solve),
`efficiency_aircraft.py` and `validate_aircraft_efficiency.py` (timing and
all-layout equivalence), `efficiency_automotive.py` and its reference and
validation companions (cold nonlinear optimizations and common checks),
`make_efficiency_figures.py` (tables and vector figures), and
`run_efficiency.py` (one entry point). `mount_springs.py` records solver
stages and permits an explicit OpenSees residual norm in the new comparison.
The shell energy and Blind-Mode Completion definitions were untouched.

The standalone entry is `python -m experiments.engineering_inverse.run_efficiency`
from the repository root after installing the project with NumPy, SciPy,
PyTorch, Matplotlib, pytest and OpenSeesPy (for example,
`python -m pip install -e . openseespy pytest`). It sets one BLAS/OpenMP/MKL
thread in its worker processes; `results/efficiency_environment.json`
records package versions and hardware. Aircraft: five methods, three meshes,
three fresh-process repeats each; raw per-repeat JSON and
`results/efficiency_aircraft_summary.csv`. Automotive: three complete cold
optimizations; raw histories, counters and stage times in method JSON plus
`results/efficiency_automotive_summary.csv`. The common reference,
stationarity, step/tolerance diagnostics and pre-abort records are also
under `results/`. All final figures are SVG/PDF in `figures/` and are
reproducible with `python -m experiments.engineering_inverse.make_efficiency_figures`.
The regression suite reports **31 passed** after the code changes.

Outstanding limitations are (i) the distinct aircraft reinforcement
rotation kinematics, (ii) no OpenSees numeric sparse-K/factorization audit,
(iii) the nonlinear automotive branch/stationarity issue, (iv) different
optimizer endpoints and reference objectives, and (v) one complete
automotive timing run per method versus three repeats for aircraft. These
limit statistical and element-only timing claims; they do not erase the
recorded reference-accepted engineering designs or the within-method
reuse/gradient measurements.
