# CMAME Round 2: driver recovery, exact sparse/rank-one production path, and revalidated experiments

## Read this first

Round 2 resolves the two implementation blockers found in the first audit and
recovers the actual manuscript benchmark definitions.  It does **not** support
an unqualified robustness or performance-superiority claim:

1. The recovered Checkerboard drivers reproduce every requested manuscript
   value.  The Round-1 plate/roof/M3 displacement studies used different
   observables, boundary-value problems, or loads and must not be cited.
2. The production tangent now retains the exact finite-configuration
   `X_M D2X_M` term while storing the only global contribution as one vector
   and one scalar.  It does not form a dense outer product.
3. The four dense boundary reconstruction matrices have been removed from the
   production path and replaced by fixed-width local stencils with identical
   formulas.
4. n=64 and n=128 now run as measured jobs.  Storage is sparse and nearly
   linear, but the current Python/AD implementation is not faster or more
   memory-efficient than the separate C++ LibShell backend.  State-count
   economy is real; end-to-end computational economy is not yet demonstrated.
5. Surface-preserving redistribution is mild for plate and roof but exposes a
   material M3 weakness: at n=12 and delta=0.30 the Checkerboard error grows
   from 2.20% to 20.12% (smooth) or 28.89% (random).  MidedgeTan also degrades,
   to 30.22% and 43.74%, respectively.  This result must not be hidden.
6. No provenance-closed MITC4 implementation was found.  None was invented.

No core shell formula, beta value, quadrature, boundary closure coefficient,
or load physics was changed.  The implementation changes are algebraically
equivalent storage/differentiation changes.  Raw CSV/JSON, all three timing
repeats, logs, vector figures, preview figures, source hashes, and baseline
adapter provenance are retained below this directory.

## Q1. Was the original manuscript driver found?

Yes.  A clean historical repository was found at
`/home/song/projects/CB-shell`, tag `v1.0-arxiv`, commit
`093ad89aa0f56b513dd89e207cfb1d69ca4a21f0`.  Its tracked benchmark modules,
archived outputs, figure scripts, and LibShell wrappers—not the target values—
were used to recover the drivers.  Source and executable SHA-256 hashes are in
[`driver_recovery.json`](00_driver_recovery/raw_data/driver_recovery.json).

The sparse/rank-one rerun gives:

| benchmark | n | recovered Checkerboard response | error against manuscript reference |
|---|---:|---:|---:|
| simply supported plate | 8 / 16 / 32 / 64 | 4.68433868e-4 / 4.50210236e-4 / 4.45377092e-4 / 4.44069708e-4 | 5.5961% / 1.4881% / 0.3986% / 0.1039% |
| Scordelis--Lo | 12 / 16 / 24 / 32 | 0.33703868 / 0.32322295 / 0.31347934 / 0.31008252 | 11.4546% / 6.8859% / 3.6638% / 2.5405% |
| MacNeal--Harder M3 | 6 / 8 / 12 / 16 | 0.054027889 / 0.078700913 / 0.091462717 / 0.093418247 | 42.2295% / 15.8473% / 2.2015% / 0.11046% |

The plate reference is `4.4360891e-4`, the roof error uses the historical
`0.3024` while retaining modern `0.3005924566`, and M3 uses `0.09352155`.
Sparse versus historical-dense responses agree exactly for plate n=8/16 and
to `1.31e-9` and `7.64e-10` relative for M3 n=8/16; the latter is the linear
solve sensitivity of the constrained ill-conditioned system, far below
manuscript plot precision.

The official LibShell `MidedgeAngleTanFormulation` at commit
`6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf` was also recovered.  Its plate and
roof curves reproduce the requested manuscript baseline numbers.  The raw
rows are in
[`midedgetan_verified_rerun.csv`](00_driver_recovery/raw_data/midedgetan_verified_rerun.csv).

## Q2. What caused the Round-1 discrepancies?

They were experimental-harness/driver errors, not formulation errors.  The
complete itemized audit is
[`benchmark_driver_discrepancy.md`](reports/benchmark_driver_discrepancy.md).

- Plate Round 1 measured a five-node, gauge-annihilating smoothed center;
  the manuscript measures the raw transverse center vertex.  At n=8 these
  are `4.50514559e-4` (1.56%) and `4.68433868e-4` (5.596%), respectively.
- Roof Round 1 solved a full-length model with different diaphragm components,
  area lumping, and a smoothed midpoint.  The manuscript solves the half-length
  domain, imposes `u_x=0` at the symmetry plane, `u_y=u_z=0` at the physical
  diaphragm, uses tensor-product load integration, and reads the raw free-edge
  node on the symmetry line.
- M3 Round 1 used raw pins instead of physical-midpoint symmetry plus derivative
  constraints, partially selected the gauge through those pins, reversed the
  point-load signs, replaced third-order `(15/8,-5/4,3/8)` point operators by
  incident-edge averages, and reported `-u_Ax` while labeling it `-u_By`.

There was no hidden n versus n+1 mesh convention.  All manuscript `n` values
count logical cells per direction.

## Q3. What grid scope does the current formulation actually have?

The current derivation requires a structured `(i,j)` topology and one constant
pair `hx,hy` per discretization.  It permits a non-rectangular or curved
reference embedding sampled over that uniform logical grid.  It is not a
derived arbitrary variable-step method.

Formula-by-formula:

1. B diagonals use fixed logical directions `(hx,hy)` and `(-hx,hy)` and one
   `ell=sqrt(hx^2+hy^2)`.
2. W normal differences are fixed centered differences divided by that `ell`.
3. Boundary moments assume B centers at `1/2,3/2,5/2`, which uniquely gives
   value weights `(15/8,-5/4,3/8)` and derivative weights `(-2,3,-1)/h`.
4. `gamma_W` contains fixed `hy/hx` and `hx/hy` factors from the auxiliary-Q1
   patch integral.
5. `X_M` uses a global alternating character but its mixed difference is still
   the fixed logical operator; geometric area weighting does not create a
   variable-step formula.
6. `X_W^rel` uses the same fixed reconstruction on complete interior W patches.

Consequently the Round-1 spacing-ratio study is **OUTSIDE CURRENT FORMULATION
SCOPE** and must not be presented as a method failure or success.  A local-step
generalization would require new difference weights, boundary Vandermonde
weights, and a re-integrated `gamma_W`; that is a future formulation change,
not an implementation patch.  See
[`formulation_scope.md`](reports/formulation_scope.md).

## Q4--Q6. Exact derivative structure and nonlinear implementation

The source energy is

`E_XM = (t/4) (8 mu) X_M^2 = c X_M^2`, with `c = 2 t mu`.

For `g = D X_M` and `H_X = D2 X_M`, its exact contribution is

`R_XM = 2 c X_M g`,

`K_XM = 2 c [g g^T + X_M H_X]`.

### Locality audit

On a finitely deformed curved M3 n=6 state:

- `X_M` equals the sum of local plaquette contributions to `1.36e-20`;
- `g` has global support (147/147 nonzeros);
- `g g^T` is fully populated but has rank exactly one;
- `H_X` has 6,669 numerical nonzeros out of 21,609 on this small patch;
- every `H_X` coupling lies within two raw-node steps in each logical direction.

Thus `D2X_M` has the same type of bounded local connectivity as the ordinary
shell tangent.  Only `g g^T` is global.  Details are in
[`xm_derivative_structure.md`](reports/xm_derivative_structure.md).

### Production representation

The new tangent is stored as

`K = K_local + alpha g g^T`, with `alpha = 2c = 4 t mu`,

where `K_local` contains ordinary membrane and bending terms,
`X_W^rel`, and the full nonlinear `2c X_M H_X` term.  Local Hessian columns are
recovered by exact AD Hessian-vector products with a conservative bounded
stencil coloring.  The global vector is not inserted as a dense outer product.

The solve uses the exact indefinite bordered system

```
[ K_local      g    ] [dq] = [rhs]
[   g^T    -1/alpha ] [lambda] [ 0 ]
```

with the mechanical linear constraints added as further KKT rows.  Its Schur
complement is exactly `K_local + alpha g g^T`; no SPD assumption or
Sherman--Morrison approximation is made.

Finite-configuration checks in
[`exactness.json`](02_low_rank_xm/raw_data/exactness.json) give:

| check | relative discrepancy |
|---|---:|
| nonlinear cantilever full tangent | 1.47e-16 |
| nonlinear cantilever tangent action | 7.72e-16 |
| nonlinear cantilever Newton increment | 2.54e-12 |
| archived nonlinear hemisphere final-state tangent | 1.37e-16 |
| archived nonlinear hemisphere tangent action | 6.36e-16 |

The retained `2c X_M H_X` Frobenius norms are nonzero: `0.01179` on the
cantilever and `73.38` on the hemisphere.  The dense and sparse cantilever load
paths both converge in 10 Newton iterations to the same final displacement and
residual.  The hemisphere comparison uses the original archived final state;
it validates the tangent there but does not rerun its full continuation.

**Answer:** the exact nonlinear tangent is strictly local sparse plus one
global rank-one update, and the nonlinear consistent term is preserved.

## Q7. Dense boundary matrices

The production `fast_boundary_geometry` no longer owns four dense
`[V x number_of_B_faces]` matrices.  Each W site stores at most a 3-by-3 list
of B-face indices and four weight lists.  Edge reconstruction, tensor-product
corners, dual fractions, and orientation signs are unchanged.

On a finite deformed curved M3 n=4 case, old dense versus new local gives zero
energy difference, `3.27e-17` relative residual difference, and `6.80e-18`
relative tangent difference.  Storage results are:

| n | old dense storage | new measured local storage | status |
|---:|---:|---:|---|
| 16 | 2.26 MiB | 0.099 MiB | both measured |
| 32 | 34.0 MiB | 0.374 MiB | both measured |
| 64 | 528.1 MiB | 1.45 MiB | both measured |
| 128 | 8.13 GiB | 5.71 MiB | old projected from exact shape; new measured |

At n=64, median geometry evaluation changes from 0.0121 s to 0.00775 s.
The historical dense code remains only as an audit artifact; it is absent from
the production geometry path.  See
[`Figure S1`](01_sparse_boundary/figures/Figure_S1_boundary_storage.pdf) and
the raw boundary JSON files.

## Q8. Production scaling and the actual cost of X_M

Every plate and M3 row at n=16,32,64,128 completed.  Each timing row has one
excluded warm-up and three preserved formal repeats in a fresh process with
eight fixed compute threads.  MidedgeTan uses a separate C++/LibShell backend;
its CPU and RSS are auxiliary, not a controlled backend comparison.

At M3 n=128:

| version | raw / active states | local nnz | augmented nnz | assembly median | factor median | peak RSS |
|---|---:|---:|---:|---:|---:|---:|
| CB-0 | 49,923 / 49,151 | 3,721,041 | 3,861,337 | 7.12 s | 18.60 s | 2,078 MiB |
| CB-M | 49,923 / 49,151 | 3,721,041 | 3,961,184 | 7.55 s | 17.38 s | 2,081 MiB |
| CB-Full | 49,923 / 49,151 | 3,721,041 | 3,961,184 | 11.49 s | 16.39 s | 2,073 MiB |
| MidedgeTan | 99,331 / 98,300 | not exposed | not exposed | different backend | different backend | 1,183 MiB |

The CB-M border adds exactly `2N+1 = 99,847` nonzeros at n=128, only 2.59%
over the constrained CB-0 system.  Evaluating `g` takes 0.094 s, 1.25% of the
CB-M assembly median.  CB-M peak RSS is 0.17% above CB-0, within process noise.
Log--log augmented-nnz slopes versus scalar states are about 1.09 for the plate
branch and 1.01 for M3, consistent with near-linear storage.

The 3V state advantage is concrete: CB-Full has 50.3% of MidedgeTan's raw
states at n=128, and M3 has 50.0% of its active states.  It does **not** yet
translate into a production memory/CPU advantage.  For M3 n=128 the Python/AD
worker peaks at 2.07 GiB versus 1.18 GiB for the C++ baseline process and takes
about 28.0 s for assembly+factor+solve versus 14.5 s total baseline process
time.  This is not a fair kernel comparison, but it rules out claiming a
measured implementation advantage.  The current M3 bottlenecks are sparse
factorization and AD local-Hessian assembly, not the rank-one vector/update.

For plate n=128, CB-Full uses 49,923 raw states, 395,641 nnz, 569 MiB RSS, and
1.87 s assembly+factor+solve; MidedgeTan uses 99,331 states, 1,255,561 nnz,
457 MiB child peak RSS, and 2.22 s total process time.  Again, backend differences
preclude a fair CPU claim.

See [`Figure C1`](03_scaling/figures/Figure_C1_time_scaling.pdf),
[`Figure C2`](03_scaling/figures/Figure_C2_memory_scaling.pdf),
[`Figure C3`](03_scaling/figures/Figure_C3_nnz_scaling.pdf), and
[`scaling.json`](03_scaling/raw_data/scaling.json).  All memory/time values in
those files are explicitly measured; the only projected boundary number is
the separately labeled old dense n=128 allocation.

## Q9. Revalidated surface-preserving distortion

The recovered manuscript geometry, load, BC, and observables were used at
plate n=16, roof n=16, and M3 n=12.  Interior samples were moved by a smooth
field and a fixed-seed (`20260922`) low-pass random tangential field at
delta=0,.05,.10,.20,.30.  Boundaries remain fixed; all nodes remain exactly on
the same plane, cylinder, or sphere; all logical cells remain valid.  The
unmodified official MidedgeTan formulation receives the identical vertex and
continuum-load tables through a provenance-recorded adapter.  Its uniform gate
reproduces all historical responses; the largest relative difference is
`2.0e-10`.

| case | method | uniform error | smooth delta=.30 | random delta=.30 |
|---|---|---:|---:|---:|
| plate n=16 | CB-Full | 1.488% | 1.601% | 1.543% |
|  | MidedgeTan | 1.248% | 1.354% | 1.267% |
| roof n=16 | CB-Full | 6.886% | 7.502% | 6.635% |
|  | MidedgeTan | 5.504% | 5.852% | 5.956% |
| M3 n=12 | CB-Full | 2.201% | 20.124% | 28.895% |
|  | MidedgeTan | 3.741% | 30.219% | 43.742% |

At delta=.30 the minimum scaled Jacobians are 0.996/0.686 (plate
smooth/random), 0.995/0.844 (roof), and 0.982/0.889 (M3); no failure was
discarded.  Smooth M3 has an area ratio 1.55 and maximum physical aspect ratio
2.98, while its scaled Jacobian remains high.  This illustrates why delta,
area ratio, aspect ratio, and scaled Jacobian are all retained rather than
cherry-picking one quality measure.

**Conclusion:** plate and roof remain stable with modest response variation;
M3 is materially redistribution-sensitive for both formulations.  Checkerboard
has smaller high-distortion M3 errors in this sample, whereas MidedgeTan has
smaller plate/roof errors.  There is no consistent evidence that exact
Varignon planarity by itself creates a general robustness advantage.

See [`Figure D1`](04_distortion_revalidated/figures/Figure_D1_distortion_comparison.pdf),
[`Figure D2`](04_distortion_revalidated/figures/Figure_D2_distorted_meshes.pdf),
and the two method-specific raw JSON/CSV files.  The old nonuniform spacing
study remains a development diagnostic outside current scope.

## Q10. Baseline status

MidedgeTan was fully recovered, uniformly verified, included in scaling, and
run on the same 30 distorted meshes.  Its driver source, build log, binary hash,
LibShell commit, and uniform gates are in [`05_baselines`](05_baselines/).

MITC4 remains a documented gap.  No trusted implementation exists in the
repositories or dependencies, so no homemade approximation is reported.  A
future comparison should use a documented Simo-type MITC4 in a validated code
with identical meshes, loads, physical BC, and observables.  See
[`BASELINE_GAP.md`](05_baselines/BASELINE_GAP.md).

## Regression and consistency checks

After optimization, the selected curved/deformed checks give:

- reference energy `8.93e-29` and residual infinity norm `1.37e-13`;
- finite rigid-translation and rigid-rotation energy changes of `1.49e-13`
  and `1.24e-13` relative on the deformed curved audit state;
- exact tangent symmetry to stored precision;
- residual directional derivative/tangent-action error `3.91e-9` at finite
  difference step `1e-7`;
- dense/sparse energy, residual, and tangent-action equivalence as above;
- identical 10-iteration nonlinear cantilever load paths;
- gauge alpha/L from 0 to 1 on undeformed M3, linearly deformed M3, and the
  nonlinear cantilever.

The rerun gauge maxima are `8.27e-15` physical midpoint motion and
`1.25e-15` relative tangent change.  On deformed states the largest relative
energy change is `1.65e-12`; the reference-state residual change normalized by
`||K|| L` is `3.57e-17`.  Thus geometry and tangent invariance remain at
roundoff after both optimizations.  Raw rows are in
[`regression_results.json`](reports/regression_results.json).

## Q11. Which Round-1 results remain usable?

| Round-1 result | Round-2 disposition |
|---|---|
| gauge invariance | **Retain.** Independently rerun after both optimizations; same roundoff-level conclusion. |
| beta_M sensitivity | **Development/Appendix only, pending recovered-driver rerun.** The algebraic target-stiffness result is valid, but benchmark response rows used reconstructed drivers. Do not cite their roof/M3 displacement insensitivity as manuscript-driver evidence. |
| curved spectrum | **Retain as a limited Appendix diagnostic.** It does not depend on the three recovered benchmark drivers; no new full-tangent near-null mode was found under its stated clamp. It is not a theorem. |
| free-boundary spectrum | **Retain as a limited Appendix diagnostic.** The exact local boundary rewrite leaves the mathematics unchanged, but the stronger 2x2 corner quotient and missing baseline spectra remain limitations. |
| distortion | **Delete/supersede.** Use only the Round-2 recovered-driver, surface-preserving comparison. Mark spacing-ratio results outside scope. |
| scaling | **Delete/supersede.** Round-1 dense projections and reference-linear prototype are replaced by measured n=16--128 exact nonlinear-capable production data. |

## Q12. CMAME manuscript recommendation

### Recommended main-text figures

1. **Round-2 Figure D1, surface-preserving distortion comparison.** Claim:
   stable solves on the sampled valid grids, modest plate/roof sensitivity,
   and a clearly disclosed M3 limitation; no universal superiority claim.
2. **Round-2 Figure C3 plus a compact C2 panel.** Claim: global completion is
   stored as an exact rank-one border, augmented nnz remains near-linear, and
   n=128 is feasible.  State that measured Python/AD peak memory is not better
   than the separate C++ baseline.
3. **Round-1 G1/G2 gauge figure with Round-2 regression values.** Claim:
   visibly different raw gauges yield identical physical geometry and
   roundoff-level internal tangent invariance at finite configurations.
4. **Existing manuscript convergence plot or Round-2 D3 state/error plot.**
   Claim: recovered Checkerboard and MidedgeTan reproduce the published
   uniform curves; 3V is a state-count economy, not automatically CPU economy.

### Appendix/Supplement

- exact tangent derivation and cantilever/hemisphere equivalence table;
- boundary storage Figure S1 and old/new equivalence;
- all timing repeats and backend caveats;
- beta sweep only with an explicit reconstructed-driver warning, or rerun it
  later with recovered drivers before making a benchmark-response claim;
- curved and free-edge spectra with their BC/quotient limitations;
- complete distortion quality metrics and all failed-case policy;
- MITC4 baseline gap.

### Delete or replace

- delete the Round-1 reconstructed plate/roof/M3 displacement tables;
- delete Round-1 D1--D3 from any manuscript draft and replace with Round-2 D1;
- delete Round-1 dense scaling projections as evidence of production cost;
- do not present spacing-ratio 2/4/8 as within the current formulation;
- do not retain claims that Checkerboard is generally faster, lower-memory,
  or more distortion-robust than standard shells.

### Claim revisions

- **Abstract:** safe claims are position-only 3V state count, exact raw-gauge
  invariance, and an exact sparse-local plus rank-one implementation.  Do not
  claim demonstrated CPU/memory superiority or general distorted-grid
  robustness.
- **Introduction:** describe Varignon planarity as an exact geometric property
  and motivation, not as proof of better distorted-mesh accuracy.  State the
  structured constant-logical-step scope.
- **Discussion:** report the M3 20--29% Checkerboard errors at delta=.30, the
  comparable/worse MidedgeTan degradation, the n=128 measured costs, and the
  current AD/factorization bottlenecks.  Separate mathematical state/storage
  structure from backend performance.

### Data-based submission position

The method is now benchmark-reproducible and its implementation matches its
mathematical low-rank structure.  That is a substantial CMAME-strengthening
result.  The defensible novelty is geometric/formulation structure plus exact
gauge and sparse-rank-one realization—not across-the-board robustness or speed.
The most important disclosed weakness is M3 sensitivity to interior spherical
sample redistribution.  A later trusted MITC4 study would strengthen external
positioning, but it should not be simulated by an unvalidated substitute.

## New summary comparison table

The machine-readable n=128 table is
[`Table_C1_summary.csv`](03_scaling/raw_data/Table_C1_summary.csv).  It reports
raw/active states, local and augmented nnz where exposed, rank-one storage,
measured peak RSS, observable, and continuum error.  Missing MidedgeTan M3 nnz
is explicitly empty because the historical driver did not expose it; no value
was inferred.
