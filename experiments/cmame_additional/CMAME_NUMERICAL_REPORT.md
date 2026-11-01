# CMAME additional numerical experiments — audit report

## Read this first: material limitations and failures

The six experimental groups have been executed and their inputs, numerical
outputs, vector/preview figures, eigenshapes, and logs are saved in this
directory. **This is not yet a reviewer-ready three-method comparison.** This
checkout contains neither the paper's original three benchmark drivers nor a
verified MidedgeTan, MITC4, IGA or CAS implementation. Its declared
dependencies contain no shell-element backend. Therefore the tables explicitly
mark both requested external baselines **UNAVAILABLE**. No conclusion that
Checkerboard is better, faster, or more robust than those methods is supported.
Do not submit Figures D1–D3 or Table C1 as a three-method comparison.

Other limitations to flag before manuscript use:

1. The three reference-displacement comparisons use newly reconstructed
   geometry, quadrilateral area-lumped dead loads and gauge-invariant physical
   observables; original paper input files were not available. All displacement
   studies except the separate cantilever are **geometrically linear solves at
   the reference tangent**, not nonlinear Newton studies. The requested
   nonlinear convergence counts therefore do not exist for those rows.
2. Smooth and random distortion at nominal n=8, plus one refinement series,
   were computed. Fixed-seed random delta=0.40 inverts a corner of **all
   three** grids (minimum scaled Jacobian −0.00635, −0.00313, −0.17286 for
   plate, roof, M3). Those three cases remain `INVALID_MESH`, not silently
   omitted. Last valid sampled random level is delta=0.30.
3. The production `fast_boundary_geometry` stores four dense
   `[V × number_of_B_faces]` reconstruction matrices and the Newton solver
   explicitly diagonalizes a dense Hessian. n=64/128 scaling jobs were stopped
   by an *a priori* memory guard; their byte projections are not timings or
   failed Newton solves. A rank-one-aware **reference-linear** prototype was
   validated separately, but it still assembles the local tangent using dense
   AD; no end-to-end large-scale improvement is claimed.
4. The free-edge spectrum needs a 2×2 raw corner block to remove mixed rigid
   and raw-gauge null modes. It is stronger than pure minimal physical BC,
   though its physical size decreases with refinement. Without baseline spectra,
   a uniquely Checkerboard-specific edge mechanism cannot be ruled out
   categorically on other boundary closures or patches.
5. An initial mesh-quality diagnostic reversed the signed corner Jacobian,
   marking all legal cells invalid; corrected only the experimental diagnostic,
   then regenerated every result. Initial random perturbations had twice the
   smooth-field RMS; matched componentwise RMS and regenerated the grid sweep.
   The first sparse spectrum attempt stalled for 6085 ARPACK iterations on a
   large zero eigenspace; subsequently used a dense symmetric LAPACK subset
   solve. An initial free-strip BC left mixed gauge/rigid zero modes, then was
   replaced by the documented 2×2 corner block. These were **experiment-harness
   issues, not changes to shell formulation**. Specific locations,
   before/after numbers, and corrected reruns are in the
   [bug audit](harness_bug_audit.json), `primary_run.log`, `spectral_run.log`,
   and the scripts; final CSVs contain the corrected numerical batches.

### Reproducibility and exact problem statements

`environment.json` records git commit, Python/PyTorch/NumPy versions, threads,
precision and random seed 20260921. `consistency_checks.json` records the
automatic checks; 20 project and harness tests passed. The documented
`README.md` at the root of this experiment directory contains exact commands.
Each subdirectory has `config/experiment.json`, `raw_data/`, `figures/`,
`logs/run.log`, `README.md`, and `artifact_manifest.json`. Figure files are
PNG, PDF and SVG; every curve has CSV or JSON backing values. The stored NPZ
files provide raw reference/deformed geometry and the first six modes for each
spectral operator. No existing paper results were overwritten.

The actual reconstructed benchmarks were:

| Case | Geometry/material/load and constraints | Measured observable and continuum target |
|---|---|---|
| Plate | Unit square, E=1, ν=.3, t=.01; area-lumped constant global −z pressure 10⁻⁸; boundary `uz=0`, minimum in-plane rigid pins. | Minus the symmetric midpoint reconstruction around the center; target 4.4360891×10⁻⁴. |
| Roof | Entire 80° cylindrical arc, R=25, axial length 50, E=4.32×10⁸, ν=0, t=.25; global −z pressure 90; end diaphragm `ux=uz=0`, one axial translation pin; two straight edges otherwise free. | Minus global `uz` at the reconstructed free-edge axial midpoint; classical 0.3024 and modern approximately 0.30059 (both retained). |
| M3 | Quarter of R=10 hemisphere with 18° polar cutout, E=6.825×10⁷, ν=.3, t=.04; opposite equator corner physical-midpoint point forces of magnitude 1, symmetry-plane `uy=0` at φ=0 and `ux=0` at φ=90°, one z pin. | Minus physical midpoint-combination x displacement at inward point B; user-provided 0.09352155. |

The roof's geometry, diaphragm constraints and −z surface load agree with the
[original shell obstacle-course description](https://ntrs.nasa.gov/api/citations/19950006750/downloads/19950006750.pdf)
and [COMSOL's verified shell model](https://doc.comsol.com/6.4/doc/com.comsol.help.models.sme.scordelis_lo_roof/scordelis_lo_roof.html).
The roof's 0.3024 historical value and differing thin-shell convergence value
are discussed in the [shell benchmark reassessment](https://link.springer.com/article/10.1007/s11831-022-09798-5).
The M3 quarter symmetry, polar opening, and ±1 quarter loads follow the
[Sandia benchmark statement](https://www.osti.gov/servlets/purl/1029819);
the exact target 0.09352155 was supplied in this task, and must be checked
against the actual paper's M3 conventions before merging numerical errors.
The finite corner observable is a gauge-invariant interpolation; it need not
match a raw nodal observable used by another code exactly.

One scalar state is each Cartesian raw vertex coordinate; total count = 3V.
"Active" means only physical essential constraints removed (including the
available gauge reduction where reported). The dense reference tangent is
solved by a positive-eigenvalue pseudoinverse to avoid forcing a solution along
an exact raw gauge. All loads and observable weights annihilate the alternating
raw gauge; this is asserted by tests, not assumed. MidedgeTan and MITC4 have
an explicit `baseline_adapter.py` contract but no invented implementation.

## Q1–Q2: distortion, nonuniform spacing, and Varignon planarity

[Raw grid and displacement data](01_mesh_distortion/raw_data/all_cases.csv)
contain **75 Checkerboard conditions** (72 converged, three inverted). At
n=8, the full method's relative continuum-target errors (%) are:

| Benchmark | uniform | smooth δ=.40 | random δ=.30 (last valid) | spacing max/min 2 / 4 / 8 | directional nx/ny=8 |
|---|---:|---:|---:|---:|---:|
| Plate | 1.56 | 5.51 | 4.40 | 6.52 / 26.37 / 49.25 | 30.99 |
| Roof | 18.76 | 22.02 | 8.34 | 6.06 / 8.63 / 25.20 | 38.77 |
| M3 | 4.04 | 5.45 | 22.01 | 8.41 / 12.26 / 13.73 | 4.43 |

These differences are **signed effects on absolute continuum error**, not
comparisons to uniform-grid Checkerboard as an artificial reference. The roof
apparently improving under random distortion is not evidence of robustness:
its undistorted n=8 error already has a pronounced directional bias, and the
distortion can accidentally compensate it. Strong smooth distortion stays
solvable with positive B/W frame determinants and no reference-normal flips.
Nonuniform spacing is a clear weakness on this mesh and closure: plate error
grows from 1.56% to 49.25% at spacing ratio 8 even though all corner scaled
Jacobians remain positive. For the M3 random level .30, the minimum scaled
Jacobian is only 0.241 and error reaches 22.0%. [D1](01_mesh_distortion/figures/Figure_D1.pdf),
[D2](01_mesh_distortion/figures/Figure_D2.pdf),
[D3](01_mesh_distortion/figures/Figure_D3.pdf) and
[D4](01_mesh_distortion/figures/Figure_D4.pdf) show these trends and representative
meshes. Check the NPZ coordinates for possible artifacts beyond a scalar error
before making any visual smoothness claim.

**Q1:** stability is demonstrated for the sampled *valid* distortions; relative
accuracy is mixed and nonuniform-grid plate accuracy degrades substantially.
Whether it is better, worse or comparable to MidedgeTan/MITC4 is **unknown**.
**Q2:** exact planarity of the edge-midpoint Varignon B face is an algebraic
property retained on these grids, but no controlled cross-method evidence shows
a *perceptible robustness advantage*. Do not claim one from these plots.

Uniform-grid convergence offers a necessary check: at n=12 the plate error
is 0.91%, roof 8.01%, and M3 0.110%. The roof's remaining 8% discrepancy and
the missing paper benchmark files counsel additional validation before using
roof distortion curves as a central CMAME result.

## Q3: exact raw-gauge invariance in finite configurations

[G1](02_gauge_invariance/figures/Figure_G1.pdf) overlays the edge midpoints of
visibly different raw M3 grids; [G2](02_gauge_invariance/figures/Figure_G2.pdf)
and [all numerical errors](02_gauge_invariance/raw_data/invariance_errors.csv)
include the undeformed curved reference, linearly deformed M3 and a **converged
nonlinear** cantilever final state (10 Newton iterations, maximum displacement
0.273). For (A=(1,0.7,-0.4)/\|(1,0.7,-0.4)\|), alpha/L=0 to 1:

| Quantity | Largest change across all states/amplitudes |
|---|---:|
| Edge midpoint physical coordinates, absolute | 7.73×10⁻¹⁵ |
| B first form, relative | 7.29×10⁻¹⁶ |
| W second form, relative | 1.16×10⁻¹⁴ |
| X_M, absolute | 2.69×10⁻¹⁶ |
| X_W^rel, absolute norm | 1.96×10⁻¹⁴ |
| Internal tangent, relative Frobenius norm | 1.23×10⁻¹⁵ |
| Internal energy, relative on deformed states | 7.34×10⁻¹³ |
| Internal residual, relative on deformed states | 1.03×10⁻¹¹ |
| Internal residual, absolute on reference M3 | 2.44×10⁻⁸ |
| Internal residual over `||K|| × characteristic length`, maximum | 3.59×10⁻¹⁷ |

**Q3:** geometric data and tangent meet the double-precision roundoff scale;
finite residual and energy are small but **not all literal 10⁻¹⁶ relative**.
Relative energy/residual differences at the stress-free reference are
mathematically undefined, not zero: their absolute/scaled numbers are given.
The direct raw dead-load potential was not compared because arbitrary raw
nodal forces need not share this gauge, although the three constructed
benchmark loads do annihilate it.

## Q4–Q5: rank-one X_M and actual computational economy

`X_M` is evaluated as one global alternating weighted scalar; PyTorch AD
differentiates its quadratic energy in both residual and explicit Hessian.
At a stress-free reference, (X_M=0), and the incremental tangent is precisely
(2c(DX_M)(DX_M)^T). The [rank-one audit](03_scaling/raw_data/rank_one_audit.json)
finds rank 1 with relative reconstruction errors 5.35×10⁻¹⁵ (plate) and
2.17×10⁻¹⁴ (M3). It fills **100% of M3 n=16 explicit tangent entries**
(CB-0 density 0.0839; CB-M and CB-Full density 1.000 at the stated
1e-12-relative nnz threshold). On the plate, physical pins leave a smaller
global support, and density grows from 0.0166 (CB-0) to 0.0448 (CB-Full).

The exact sparse augmented reference-linear representation uses the **same
energy and tangent**, placing only one gradient vector in an augmented
saddle-point system. [Equivalence data](03_scaling/raw_data/low_rank_equivalence.csv):

| n=16 | Dense tangent | Sparse augmented storage | Dense nnz → augmented nnz | Relative displacement disagreement | Naive solve median → sparse factor+update median |
|---|---:|---:|---:|---:|---:|
| Plate | 5.09 MB | 0.154 MB | 28,364 → 12,581 | 6.15×10⁻¹⁵ | 0.0213 s → 0.00534 s |
| M3 | 5.52 MB | 0.719 MB | 690,561 → 59,672 | 3.66×10⁻¹¹ | 0.0254 s → 0.0153 s |

For M3 the sparse solve's relative linear residual is 2.86×10⁻⁸ (ill
conditioning/finite precision), compared with near-machine displacement
agreement on plate; M3 energy differs by about 1.43×10⁻¹⁰ in its units.
This is a demonstrably equivalent **reference-tangent** representation, not a
production nonlinear solver. Away from reference, differentiation also yields
the term (2c X_M D^2X_M); an exact large-deformation implementation must
keep it, not simply reuse reference (gg^T). The dense boundary operators and
the dense-AD assembly remain unresolved.

Measured n=16 post-warmup explicit tangent assembly medians (plate CB-0/M/Full)
are approximately 0.417/1.126/2.770 s; matrix-free tangent actions are
0.0279/0.0339/0.0458 s. At n=32 CB-Full tangent actions are 0.112 s (plate)
and 0.0881 s (M3). These are *component timings* in one backend, not fair
external CPU benchmarks. Four boundary reconstruction matrices alone require
about 554 MB at n=64 and 8.72 GB at n=128; dense reduced tangents would be
roughly 1.23 GB and 19.5 GB for the plate respectively, excluding AD graphs
and factorizations. [C1](03_scaling/figures/Figure_C1.pdf),
[C2](03_scaling/figures/Figure_C2.pdf),
[C3](03_scaling/figures/Figure_C3.pdf) and
[Table C1](03_scaling/raw_data/Table_C1.csv) distinguish measured and projected
data.

**Q4:** yes, the current explicit Hessian is dense from global X_M on a curved
patch; an exactly equivalent low-rank reference implementation avoids
*representing the rank-one term densely*, but does **not yet** remove other
production dense costs. **Q5:** 3V is confirmed as a *state count*, but no
MidedgeTan/MITC4 DOF or CPU baselines exist and no end-to-end sparse solver was
measured. Therefore computational economy relative to standard shells is
**unproven**, especially for n≥64.

## Q6: β_M sensitivity, no-benchmark-tuning check

The [raw sweep](04_beta_sensitivity/raw_data/beta_sensitivity.csv) and
[B1](04_beta_sensitivity/figures/Figure_B1.pdf),
[B2](04_beta_sensitivity/figures/Figure_B2.pdf),
[B3](04_beta_sensitivity/figures/Figure_B3.pdf) include c_beta=0,.25,.5,1,2,4,8.
The largest *relative change of the load-bearing displacement versus c=1*
among c=.25–4 is:

| Case | n=8 | n=12 |
|---|---:|---:|
| M3 | 2.62×10⁻¹⁰ | 6.16×10⁻¹⁰ |
| Roof | 1.38×10⁻⁴ | 1.09×10⁻⁵ |

The **completion-added** z_M quadratic stiffness is exactly linear in c_beta
by construction; the total directional stiffness also contains other terms.
In these constrained finite meshes, setting c=0 lowers a positive eigenvalue
(M3 n=8: 1.166 to 0.110; n=12: 0.626 to 0.0573) but does **not** generate
an additional exact nullity under the present physical BC. Thus the heuristic
claim that an exact membrane blind mode “reappears” at c=0 is **not verified
under these BC**. The broad, flat displacement response does support that
c=1 is not a narrow benchmark-tuned optimum for these two *linear* observables;
it does not establish parameter insensitivity for nonlinear or different loads.

## Q7: positive/negative curvature blind-mode spectra

The [positive-K sphere and negative-K saddle spectra](05_curved_spectrum/raw_data/spectra.csv)
contain 20 eigenvalues per operator per n=8/16/32 and saved mode vectors.
With the documented two-layer perimeter essential clamp, the full-tangent
smallest values are positive:

| Geometry | n=8 | n=16 | n=32 |
|---|---:|---:|---:|
| Sphere | 7.64×10⁻⁴ | 2.05×10⁻⁴ | 5.01×10⁻⁵ |
| Saddle | 1.52×10⁻³ | 3.64×10⁻⁴ | 8.74×10⁻⁵ |

Their ratios to the largest tangent eigenvalue fall approximately by four
per refinement (e.g. sphere 3.96×10⁻³ → 9.49×10⁻⁴ → 2.24×10⁻⁴), consistent
with physical refinement scaling; no full-tangent mode lies below the declared
relative near-null cutoff 10⁻⁹. The alternating score alone grows in some
positive modes (up to .67 sphere and .71 saddle at n=32); **do not call these
spurious merely from their pattern**, because their stiffness is positive and
the strong clamp changes the operator's nullspace. [K1](05_curved_spectrum/figures/Figure_K1.pdf),
[K2](05_curved_spectrum/figures/Figure_K2.pdf),
[stored soft-mode shapes](05_curved_spectrum/figures/Figure_K3.pdf) and
[Table K1](05_curved_spectrum/raw_data/Table_K1.csv) support the finite-patch
statement only. Some isolated observation operators have zero/near-zero
channels, as expected; the *combined tangent* is what answers the physical
mechanism question.

**Q7:** no new clearly physical, persistent checkerboard-like **zero** mode
was found in the tested clamped patches; this is not a general curved-surface
theorem or evidence for unconstrained curved shells.

## Q8: free-edge localized modes

[Free-strip spectra and mode metrics](06_boundary_modes/raw_data/spectra.csv)
show positive full-tangent minima under the final quotient-space corner clamp:

| Strip | n=8 | n=16 | n=32 | alternating score of n=32 softest |
|---|---:|---:|---:|---:|
| Flat | 1.44×10⁻⁷ | 2.67×10⁻⁸ | 5.12×10⁻⁹ | 5.08×10⁻⁵ |
| Cylinder | 1.49×10⁻⁷ | 2.41×10⁻⁸ | 4.32×10⁻⁹ | 6.16×10⁻⁴ |

For the softest full modes, the two-layer boundary displacement fraction
**decreases** from .898 to .535 (flat) and .876 to .517 (cylinder) as n
increases. Boundary moment reconstruction and X_W^rel completed interior
patches do not create an observed persistent exact checkerboard/free-edge
zero mode in these two sampled strips. Pure bending observations have null
spaces and should not be mistaken for the full shell tangent. [E1](06_boundary_modes/figures/Figure_E1.pdf),
[E2](06_boundary_modes/figures/Figure_E2.pdf), and
[E3](06_boundary_modes/figures/Figure_E3.pdf) contain the refinement and saved
representative shapes. **Q8:** the specific hypothesized edge mode was not
observed, but the stronger-than-minimal corner clamp and missing baseline
spectra prevent a universal negative conclusion.

## Cross-cutting checks and recorded implementation issues

At the curved M3 reference, internal energy is 1.27×10⁻²⁸ and residual norm
5.82×10⁻¹³; rigid translation energy 8.07×10⁻²⁴ and randomly oriented finite
rotation energy 2.04×10⁻²³; tangent relative asymmetry 1.61×10⁻¹⁶; selected
finite-difference directional tangent/residual error 2.13×10⁻⁹. Reference
B-frame metric determinant stays positive (minimum .152); reference normal
flips were not found. A selected alpha=1 raw-gauge check, normalized to
appropriate geometric/tangent scales, has largest finite error 4.13×10⁻¹⁵.
Additional per-distortion B/W determinant and normal-flip values are in the
case CSV. These checks passed; they are **selected cases**, not exhaustive
automatic checks at every configuration and load step.

No core implementation bug was silently fixed. The three issues listed in
the opening section belong to the experiment harness and eigensolver/constraint
selection. On the accepted energy side, the **dense boundary operator storage
and dense explicit Hessian are concrete implementation scalability deficits**;
they have not been changed in `src/`. A full production low-rank extension
must preserve the nonlinear (X_M D^2X_M) term and the paper's reference
relative curvature; the prototype does not yet do that.

## Q9–Q10: manuscript selection and submission advice

Recommended **conditional** main-text figures (3–5) and the claims each could
support:

1. [G1](02_gauge_invariance/figures/Figure_G1.pdf) +
   [G2](02_gauge_invariance/figures/Figure_G2.pdf) as one compact paired
   figure: visually different raw configurations yield the same physical
   midpoint geometry; internal tangent invariance is close to roundoff while
   residual error is reported honestly.
2. [D1](01_mesh_distortion/figures/Figure_D1.pdf) together with
   [D2](01_mesh_distortion/figures/Figure_D2.pdf): stability on sampled legal
   distorted grids **and** loss of accuracy on severe nonuniform spacing.
   Main text only after original benchmark settings and at least MidedgeTan
   (preferably MITC4) are added at identical meshes.
3. [C3](03_scaling/figures/Figure_C3.pdf) with the low-rank storage rows from
   [Table C1](03_scaling/raw_data/Table_C1.csv): global X_M makes the current
   curved explicit tangent dense, and the algebraically exact rank-one
   reference representation removes that representation cost. **Do not** claim
   production CPU/peak-memory advantage until the dense boundary and nonlinear
   assembly are optimized and benchmarked.
4. [B1](04_beta_sensitivity/figures/Figure_B1.pdf) with
   [B2](04_beta_sensitivity/figures/Figure_B2.pdf): normal *linear* response is
   broad in β while the rank-one target stiffness rises linearly. Appropriate
   for Discussion/Appendix if the paper already contains equivalent parameter
   material.
5. [K1](05_curved_spectrum/figures/Figure_K1.pdf) and
   [K2](05_curved_spectrum/figures/Figure_K2.pdf), or a two-panel subset:
   finite clamped positive/negative-K patches have no additional full-tangent
   near-null mode under the stated threshold; keep full spectra and E-series
   boundary modes in Supplement.

The new [summary comparison table](summary_comparison_table.csv) is a
**transparent placeholder**: n=12 Checkerboard displacement/state/error and
explicit UNAVAILABLE baseline rows. It must be populated from validated
MidedgeTan/MITC4 adapters before it can serve as a CMAME quantitative
comparison. Detailed gauge, β, K and E tables belong in Appendix/Supplement;
failed/inverted cases and complete raw data must remain available.

Recommended claim revisions:

- **Abstract:** retain “position-only 3V states” and observed finite raw-gauge
  invariance; remove “more robust on distorted/nonuniform meshes,” “faster
  than standard shells,” and any universal curved/free-edge stability claim.
- **Introduction:** say exact Varignon planarity is a geometric guarantee and
  motivates testing robustness, not that it automatically improves error.
- **Discussion:** explicitly state the 49% nonuniform plate error, the random
  M3 error increase to 22%, dense tangent/boundary operators, missing external
  baselines, and the limited BC/refinement range of spectral conclusions.
  Separate mathematical 3V state economy from measured implementation cost.

**Submission decision (Q10):** the independent gauge demonstration and honest
β/spectral appendices are valuable, but a CMAME submission claiming comparative
robustness or computational economy is **premature** with this checkout.
Recover the validated paper benchmark drivers, integrate a trustworthy
MidedgeTan implementation and (if available) MITC4 on the *same* distorted
grids, then implement sparse boundary reconstruction and an exact nonlinear
rank-one-aware tangent before drawing timing conclusions. The most conspicuous
present weakness is severe nonuniform-grid plate degradation **despite**
positive Jacobians; it deserves analysis rather than tuning β, quadrature,
stencil or load definitions to conceal it.
