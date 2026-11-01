# CMAME Round 4: completion of the manuscript experiments

> Date: 2026-10-26  
> Scope: close the experimental gaps in `docs/CMAME_PROJECT_SYNTHESIS.md`
> without changing the Checkerboard formulation or any accepted coefficient.

## Read first

This round adds seven things that the arXiv manuscript and the Round-2/3
evidence package did not yet have in one place:

1. a **full 3V simply supported plate** solve on the production energy, with a
   fair active-DOF comparison to MidedgeTan;
2. a **single-adapter, measured** 12-row MidedgeTan comparison;
3. a **dual roof reference** (historical `0.3024` and modern `0.3005924566`);
4. **reruns of the three finite-deformation benchmarks** with the current
   `checkerboard_shell.mechanics.energy.paper_components`;
5. **exact production-tangent checks** at the nonlinear final states of all
   three benchmarks;
6. a **documented MITC4 standard-shell baseline** via OpenSeesPy
   `ShellMITC4` (6 DOF/node);
7. **M3 dual observables** — raw equator vertex (main, same as MidedgeTan)
   and physical-midpoint reconstruction (secondary).

The nonlinear reruns reproduce the archived arXiv values to solver precision.
No shell formula, coefficient, quadrature, boundary closure, or load physics
was changed.

The MITC4 comparison is **mixed**, and this is the main new caveat for the
paper: MITC4 is substantially more accurate on the plate (n=8--32) and roof at
the same mesh, but uses 6 DOF per node; Checkerboard is better on M3 at
n=12--16 and on the plate at n=64.  A blanket “competitive accuracy per state”
claim is no longer supported; the defensible claim is the 3V state, the
blind-mode completion, the sparse rank-one tangent, and competitive accuracy
on specific problems.

## 1. Full 3V simply supported plate

The recovered manuscript plate driver is the scalar transverse branch: it
solves only the vertical DOF at interior nodes (`n=16` active `225`).  The full
production model is solved here with `w=0` on all four edges, the in-plane
field free except for three rigid pins, and the same interior nodal load.

| n | full 3V response | full error | scalar response | scalar error | MidedgeTan response | MT error | full active | MT active |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 8 | 4.68218e-4 | 5.5475% | 4.68434e-4 | 5.5961% | 4.65579e-4 | 4.9526% | 208 | 257 |
| 16 | 4.50032e-4 | 1.4479% | 4.50210e-4 | 1.4881% | 4.49145e-4 | 1.2480% | 800 | 1025 |
| 32 | 4.45208e-4 | 0.3605% | 4.45377e-4 | 0.3986% | 4.44996e-4 | 0.3126% | 3136 | 4097 |
| 64 | 4.43903e-4 | **0.0664%** | 4.44070e-4 | 0.1039% | 4.43956e-4 | 0.0782% | 12416 | 16385 |

Observations:

- The full 3V and the scalar branch differ by a nearly constant `1.7e-7`
  (about `0.04%` of the response).  The difference is systematic, not solver
  noise; the full production model is slightly softer.  The scalar branch is
  an independent hand-derived operator and is retained only as a cross-check.
- The proper comparison is **raw states** (3V) or **active all-DOF states**.
  Checkerboard uses 50.5--53.9% of MidedgeTan raw states, and 75.8--80.9% of
  its active states on this plate.  The earlier `225 vs 1025` active-state
  comparison was not apples-to-apples.
- At `n=64`, the full 3V Checkerboard error (0.0664%) is slightly smaller than
  MidedgeTan (0.0782%); at coarse meshes MidedgeTan is slightly better.
- The interior nodal load sums to `(n-1)^2/n^2 * 1e-8`; the boundary strip load
  is carried by the supports.  This is the same load definition as the
  recovered scalar driver.  The full 3V in-plane residual is about `3e-7`
  (0.07% of the transverse response), consistent with the discrete
  membrane--bending coupling rather than a constraint error.

Figure: `figures/Figure_R4_full_plate.png`.

## 2. Unified MidedgeTan comparison

All MidedgeTan rows below were rerun through the rebuilt official LibShell
`MidedgeAngleTanFormulation` adapter on the same uniform vertex and load
tables.  Raw and active states are measured from the adapter output.

| benchmark | n | Checkerboard response (error) | MidedgeTan response (error) | CB raw/active | MT raw/active |
|---|---:|---:|---:|---:|---:|
| plate | 8 | 4.68218e-4 (5.5475%) | 4.65579e-4 (4.9526%) | 243/208 | 451/257 |
| plate | 16 | 4.50032e-4 (1.4479%) | 4.49145e-4 (1.2480%) | 867/800 | 1667/1025 |
| plate | 32 | 4.45208e-4 (0.3605%) | 4.44996e-4 (0.3126%) | 3267/3136 | 6403/4097 |
| plate | 64 | 4.43903e-4 (0.0664%) | 4.43956e-4 (0.0782%) | 12675/12416 | 25091/16385 |
| roof | 12 | 0.337039 (11.455% / 12.125%) | 0.274884 (9.099% / 8.553%) | 507/468 | 963/912 |
| roof | 16 | 0.323223 (6.886% / 7.529%) | 0.285756 (5.504% / 4.936%) | 867/816 | 1667/1600 |
| roof | 24 | 0.313479 (3.664% / 4.287%) | 0.295748 (2.200% / 1.612%) | 1875/1800 | 3651/3552 |
| roof | 32 | 0.310083 (2.541% / 3.157%) | 0.299811 (0.856% / 0.260%) | 3267/3168 | 6403/6272 |
| M3 | 6 | 0.0530257 (43.301% / 42.230%) | 0.0531934 (43.12%) | 147/107 | 267/212 |
| M3 | 8 | 0.0777953 (16.816% / 15.847%) | 0.0811630 (13.21%) | 243/191 | 451/380 |
| M3 | 12 | 0.0909082 (2.794% / 2.201%) | 0.0970198 (3.741%) | 507/431 | 963/860 |
| M3 | 16 | 0.0930731 (0.479% / 0.1105%) | 0.0982063 (5.009%) | 867/767 | 1667/1532 |

For roof, the two errors are relative to `0.3024` and `0.3005924566`,
respectively.  For M3, the two errors are the raw equator vertex (main, the
same discrete observable as MidedgeTan) and the physical-midpoint moment
reconstruction (second), respectively.

What the unified run changes:

- The plate rows now use the complete 3V production solve.
- Every MidedgeTan row is measured, not taken from a historical CSV.
- M3 raw states are measured from the adapter (`963` at `n=12`), not inferred
  from `6n^2+8n+3`.
- The MidedgeTan M3 sequence is confirmed non-monotonic with measured states:
  `3.741%` at `n=12` and `5.009%` at `n=16`.  The cause is not isolated in
  this round (mesh, discrete boundary/director treatment, or the baseline's
  own convergence behavior); it is a property of the baseline comparison and
  should be stated explicitly in the paper.
- M3 is now compared at the **raw equator vertex** for both methods.  The
  Checkerboard physical-midpoint reconstruction is reported in parentheses:
  at `n=16` the raw observable gives `0.4795%` and the midpoint reconstruction
  gives `0.1105%`, versus MidedgeTan `5.009%`.  The conclusion is unchanged,
  but the headline number is no longer observable-mismatched.
- The roof result is honest only with the modern reference: Checkerboard is
  `3.157%` high at `n=32` against `0.3005924566`, while MidedgeTan is
  `-0.260%`.  The historical `0.3024` flatters Checkerboard by about `0.6`
  percentage points.

Figure: `figures/Figure_R4_unified_comparison.png`.

## 3. Nonlinear benchmarks on the production energy

The archived drivers were ported to call
`checkerboard_shell.mechanics.energy.paper_components` and
`make_targets`.  At a representative state the ported energy agrees with the
archived implementation to `5.6e-17` absolute, and every reported metric below
reproduces the archived value to solver precision.

### 3.1 End-moment cantilever

| mesh | max turning angle | RMS path error | max path error | final tip error | status |
|---|---:|---:|---:|---:|---|
| 16x4 | 6.49756 rad (372.3°) | 9.886% | 13.255% | 13.255% | VALIDATED |
| 24x4 | 6.29911 rad (360.9°) | 3.894% | 5.135% | 5.135% | VALIDATED |
| 32x4 | 6.28042 rad (**359.841°**) | **2.442%** | **2.792%** | **2.792%** | VALIDATED |

The `32x4` values match the arXiv text (`359.841°`, RMS `2.442%`, max
`2.792%`).  The final membrane-energy fraction is `0.432%`; `X_M` and
`X_W^rel` are numerically inactive on this pure-bending branch.

### 3.2 Alternating-force hemisphere

| n | final A error | final B error | joint RMS | max | steps | status |
|---:|---:|---:|---:|---:|---:|---|
| 6 | 38.242% | 30.722% | 34.934% | 30.613% | 62 | VALIDATED_WITH_FINITE_MESH_ERROR |
| 8 | 19.844% | 9.712% | 12.993% | 13.371% | 80 | VALIDATED_WITH_FINITE_MESH_ERROR |
| 12 | 6.765% | 1.223% | 3.413% | 4.559% | 84 | VALIDATED_WITH_FINITE_MESH_ERROR |
| 16 | 2.759% | 0.149% | 1.736% | 1.859% | 85 | VALIDATED_WITH_FINITE_MESH_ERROR |

The `n=12` values match the archived `6.765% / 1.223% / 3.413% / 4.559%`;
the `n=16` values match the arXiv text (`2.759% / 0.149% / 1.736% / 1.859%`).

### 3.3 Open-ended cylinder pull-out

| mesh | final w_A error | final -u_B error | final -u_C error | combined RMS | max | steps | status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| 8x8 | 18.935% | 14.078% | 5.093% | 13.648% | 18.371% | 45 | VALIDATED_WITH_FINITE_MESH_ERROR |
| 12x12 | 11.509% | 8.050% | 3.039% | 7.340% | 9.781% | 45 | VALIDATED_WITH_FINITE_MESH_ERROR |
| 16x16 | 8.094% | 5.316% | 1.861% | 4.483% | 5.316% | 46 | VALIDATED_WITH_FINITE_MESH_ERROR |

The `8x8` values match the archived `18.935% / 14.078% / 5.093% / 13.648%`.
The `12x12` values match the archived `11.509% / 8.050% / 3.039% / 7.340%`; the `16x16` values match the archived `8.094% / 5.316% / 1.861% / 4.483%` (the paper text reports the same fine-mesh final errors and combined RMS).

### 3.4 Exact production tangent at nonlinear final states

The archived final configurations were reused only as states.  At each state,
the current production tangent `K_local + alpha g g^T` was reduced to the
benchmark's kinematic subspace and compared with a dense AD Hessian of the
current production internal energy.

| case | reduced dimension | relative tangent difference | `X_M` at final state | rank-one storage (full state) |
|---|---:|---:|---:|---:|
| cantilever 32x4 | 62 | **1.64e-16** | -6.2e-35 | 496 |
| hemisphere n=12 | 431 | **1.49e-15** | -9.66e-6 | 508 |
| open cylinder 8x8 | 168 | **4.95e-15** | 1.78e-4 | 244 |

This closes the “production rank-one tangent only checked on a small
cantilever” gap: the exact nonlinear tangent is verified at a nearly
360-degree cantilever state, a strongly pinched curved hemisphere, and a
localized open-cylinder pull-out state.

## 4. MITC4 standard-shell baseline

### 4.1 Implementation and provenance

The standard-shell baseline uses OpenSeesPy **3.7.1**, element
`ShellMITC4` with `ElasticMembranePlateSection`, 6 DOF per node (three
translations and three rotations).  This is a documented MITC4
implementation, not a home-made element.  The same structured quad meshes,
materials, loads, and physical boundary conditions as the recovered drivers
are used.

The element is a **Mindlin--Reissner** shell.  Its plate response therefore
contains a transverse-shear contribution that the Kirchhoff--Love
Checkerboard and MidedgeTan models do not have.  This is reported explicitly;
it is not a locking error of MITC4.

### 4.2 Plate

| n | MITC4 response | MITC4 error | raw states | CB / MT error | CB / MT raw states |
|---:|---:|---:|---:|---:|---:|
| 8 | 4.41884e-4 | 0.389% | 486 | 5.547% / 4.953% | 243 / 451 |
| 16 | 4.43931e-4 | 0.073% | 1734 | 1.448% / 1.248% | 867 / 1667 |
| 32 | 4.44953e-4 | 0.303% | 6534 | 0.361% / 0.313% | 3267 / 6403 |
| 64 | 4.45982e-4 | 0.535% | 25350 | 0.066% / 0.078% | 12675 / 25091 |

The MITC4 error decreases to 0.073% at n=16 and then rises to 0.535% at n=64.
This is the Mindlin shear floor, not a convergence failure.  A thin-limit
check confirms it:

| thickness | n=32 | n=64 |
|---|---:|---:|
| t=0.01 | 0.303% | 0.535% |
| t=0.001 | 0.028% | 0.001% |
| t=0.0001 | 0.031% | 0.008% |

At `t=0.001` the MITC4 plate error is 0.001% at n=64; the 0.5% floor at
`t=0.01` is the shear deformation omitted by the thin-plate reference.

### 4.3 Roof

| n | MITC4 historical error | MITC4 modern error | CB historical / modern | MT historical / modern | raw states (CB/MT/MITC4) |
|---:|---:|---:|---:|---:|---:|
| 12 | 0.294% | 0.306% | 11.455% / 12.125% | 9.099% / 8.553% | 507 / 963 / 1014 |
| 16 | 0.497% | 0.101% | 6.886% / 7.529% | 5.504% / 4.936% | 867 / 1667 / 1734 |
| 24 | 0.516% | 0.083% | 3.664% / 4.287% | 2.200% / 1.612% | 1875 / 3651 / 3750 |
| 32 | 0.481% | 0.118% | 2.541% / 3.157% | 0.856% / 0.260% | 3267 / 6403 / 6534 |

MITC4 converges to `0.30084--0.30095`, i.e. `0.08--0.12%` from the modern
reference and `0.48--0.52%` from the historical one.  This is an independent
confirmation that the modern reference is the correct converged value and
that the Checkerboard roof response is genuinely high by about 3% at n=32.
MITC4 is more accurate than both Checkerboard and MidedgeTan on this problem,
at twice the raw state count of Checkerboard.

### 4.4 M3

| n | MITC4 raw error | MITC4 smoothed error | CB raw / midpoint | MT raw | raw states (CB/MT/MITC4) |
|---:|---:|---:|---:|---:|---:|
| 6 | 59.379% | 59.055% | 43.301% / 42.229% | 43.122% | 147 / 267 / 294 |
| 8 | 32.787% | 32.737% | 16.816% / 15.847% | 13.215% | 243 / 451 / 486 |
| 12 | 9.478% | 9.628% | 2.794% / 2.201% | 3.741% | 507 / 963 / 1014 |
| 16 | 3.615% | 3.769% | 0.479% / 0.110% | 5.009% | 867 / 1667 / 1734 |
| 24 | 1.038% | 1.144% | — | — | 3750 |
| 32 | 0.472% | 0.545% | — | — | 6534 |
| 48 | 0.155% | — | — | — | 14406 |
| 64 | 0.053% | — | — | — | 25350 |

MITC4 is **worse** than Checkerboard and MidedgeTan on the coarse M3 meshes
used in the manuscript.  It converges slowly because the raw corner point
load at B is singular: 0.47% at n=32, 0.155% at n=48, and 0.053% at n=64.
The Checkerboard raw-vertex value at n=16 is already 0.48%.  A smoothed
physical-point load does not change this conclusion.  MITC4 therefore does
not provide a coarse-mesh advantage on M3; the Checkerboard/MidedgeTan point
operators do.

### 4.5 What the MITC4 comparison changes

- **Do not claim uniform accuracy superiority.**  MITC4 is substantially more
  accurate on the plate at n=8--32 and on the roof at all manuscript meshes;
  Checkerboard is better on M3 at n=12--16 and on the plate at n=64.
- **Do claim a state-count difference.**  MITC4 uses 6 DOF per node; the
  Checkerboard uses 3, i.e. half.  Accuracy per state is problem-dependent,
  not uniformly better.
- **Do use MITC4 to settle the roof reference.**  Its converged modern value
  independently supports `0.30059` and confirms the Checkerboard roof
  limitation.
- **Do state the theory difference.**  MITC4 is Mindlin--Reissner; the
  Checkerboard and MidedgeTan are Kirchhoff--Love/Koiter.  The MITC4 plate
  shear floor is not a locking failure.

Combined table: `COMBINED_COMPARISON_TABLE.md`.  
Figure: `figures/Figure_R4_combined_three_methods.png`.

## 5. Figures and tables

- `figures/Figure_R4_full_plate.png` — full 3V vs scalar vs MidedgeTan.
- `figures/Figure_R4_unified_comparison.png` — 12-row error vs raw states.
- `figures/Figure_R4_nonlinear_paths.png` — cantilever, hemisphere, open
  cylinder paths.
- `UNIFIED_COMPARISON_TABLE.md`, `MITC4_TABLE.md`,
  `COMBINED_COMPARISON_TABLE.md`, `NONLINEAR_SUMMARY.md` — machine-readable
  summaries.
- `figures/Figure_R4_combined_three_methods.png` — Checkerboard /
  MidedgeTan / MITC4 error vs raw states.
- `raw_data/tangent_final_states.json` — exact tangent checks at nonlinear
  final states.

## 6. What changes in the manuscript

### 6.1 Recommended replacements

- **Plate table/figure:** use the full 3V rows as the main result; keep the
  scalar branch only as a cross-check in the supplement.  Replace the “225
  active states” comparison with the full 3V active count (`208/800/3136/12416`).
- **Roof table:** show both references, or use the modern `0.3005924566` and
  state why.  Do not report only the `0.3024` error.
- **M3 table:** main comparison uses the raw equator vertex for all methods;
  the Checkerboard physical-midpoint value is reported in parentheses.  Add a
  footnote that the MidedgeTan M3 sequence is non-monotonic.
- **Nonlinear section:** cite the Round-4 reruns as “current production energy”
  and keep the archived reference curves.
- **Add a MITC4 baseline table.**  Report it as a 6-DOF Mindlin--Reissner
  standard element; state that it is more accurate on the plate/roof but uses
  twice the raw states, and that it converges slowly on M3 because of the raw
  corner point load.

### 6.2 Claims that are now better supported

- The full 3V position-only state is competitive with MidedgeTan, and uses
  `75.8--80.9%` of its active states on the plate.
- The nonlinear benchmarks are reproducible with the current package, not only
  with the historical repository.
- The MidedgeTan comparison is measured, not inferred, for all 12 rows.
- MITC4 independently confirms the modern roof reference and provides a
  standard-shell accuracy/state-count reference.  The comparison is mixed:
  MITC4 wins plate (n=8--32) and roof; Checkerboard wins M3 (n=12--16) and
  plate n=64.

### 6.3 Claims that remain unsupported

- General robustness under surface-preserving redistribution: M3 still shows
  `20.12%` (smooth) and `28.89%` (random) at `delta=0.30`.
- CPU or peak-memory superiority: the current Python/AD implementation is
  slower and uses more RSS than the C++ LibShell baseline.
- Uniform accuracy superiority over MITC4: the MITC4 comparison is mixed and
  theory-dependent (Mindlin vs Kirchhoff--Love).
- Arbitrary variable spacing / unstructured topology: outside the constant
  `hx,hy` formulation.
- General curved-shell locking-free convergence: only the aligned
  generalized-cylinder `O(h^2)` result is proven.

## 7. Remaining engineering

1. **MITC4 baseline — done.**  OpenSeesPy `ShellMITC4` with
   `ElasticMembranePlateSection` (6 DOF/node) is used as a documented
   standard-shell baseline.  The results and provenance are in
   `raw_data/mitc4.{csv,json}` and `MITC4_TABLE.md`.
2. **M3 observable unification — done.**  The main table now uses the raw
   equator vertex for both methods; the Checkerboard physical-midpoint value
   is retained as a second observable.
3. **Plate load normalization.**  The recovered scalar driver applies
   `1e-8/n^2` at interior nodes, i.e. the interior part of the pressure load;
   boundary load goes directly into the supports.  The full 3V run uses the
   same load.  For the paper, state the nodal load definition explicitly.
4. **General gauge quotient API.**  The package still requires benchmark-side
   gauge rows.  Expose the rigid+gauge basis / orthogonal complement in the
   public API, or document the construction.
5. **Clean release.**  Commit `experiments/cmame_round4`, tag the workspace,
   update `README.md`, and attach `artifact_manifest.json` to the submission.

## 8. Provenance

- LibShell: `evouga/libshell`, commit
  `6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf`, rebuilt from the pinned
  upstream tarball.
- Eigen: 3.4.0 headers.
- MITC4: OpenSeesPy 3.7.1 (Zhu, McKenna & Scott, *SoftwareX* 19 (2022)
  101195), element `ShellMITC4`, section `ElasticMembranePlateSection`,
  6 DOF per node.  Installed in `/tmp/openseespy_pkg`; recorded in
  `logs/environment.json`.  The same structured meshes, materials, loads, and
  physical BCs are used.  `ShellMITC4` is the Dvorkin--Bathe MITC4 element.
- The adapter reproduces the Round-2 uniform gates (plate n=16, roof n=16,
  M3 n=12) to better than `1e-9` relative.
- Environment, git state, binary hash, and every artifact hash are in
  `logs/environment.json` and `artifact_manifest.json`.
