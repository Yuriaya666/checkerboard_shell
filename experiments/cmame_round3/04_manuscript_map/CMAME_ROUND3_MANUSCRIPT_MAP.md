# CMAME Round 3 manuscript revision map

This map is keyed to the current English manuscript
`/home/song/projects/CB-shell/arxiv_submission/Checkerboard_Shells.tex`.  It
does not alter the formulation.  It closes the manuscript around recovered
benchmark drivers, exact sparse/rank-one implementation evidence, and the
observed M3 redistribution limitation.

## Abstract

### Retain

- The exact planar Varignon midpoint geometry.
- The physical midpoint quotient and position-only unknowns.
- The flat blind-mode classification and the definitions of `X_M` and
  `X_W^rel`.
- The limited aligned-generalized-cylinder `O(h^2)` approximation statement.

### Revise

- Replace “`X_M` removes the target membrane defect while leaving standard
  load-bearing responses essentially unchanged” by a quantified statement:
  “On the recovered plate n=16, roof n=16, and M3 n=12 drivers, varying
  `beta_M/(8mu)` from 0.25 to 4 changes the reported displacement by at most
  `1.63e-5` relative to the nominal response; a flat membrane-operator test
  recovers one additional null direction when the completion is removed.”
- If “compact” or “efficient” is used, attach it only to state count or sparse
  algebraic structure, not runtime.

### Delete

- Any claim that the method is generally robust under distorted meshes.
- Any implication that exact B-face planarity makes approximation error
  distortion-insensitive.
- “locking-free” without the fixed theorem qualifier.

### Add

- One sentence: “The exact tangent is implemented as a bounded-stencil sparse
  local part plus one global rank-one update, without forming a dense outer
  product.”
- If space permits, one limitation phrase: “Surface-preserving spherical
  redistribution can materially degrade the M3 response.”

## Introduction

### Retain

- Motivation from Varignon planarity, common midpoint carrier, and the absence
  of independent director/rotation states.
- Comparison of variable types with mixed, IGA, and midedge approaches.
- The contribution “competitive accuracy per state,” now backed by the final
  recovered-driver table.

### Revise

- At the current sentence “If both geometric quantities can be evaluated
  robustly from current nodal positions,” replace “robustly” by “unambiguously
  and objectively.”  The proven property is geometric well-posedness, not
  approximation robustness.
- State the domain explicitly: structured logical topology with constant
  `hx,hy` per discretization.  Arbitrary variable parameter spacing is not
  derived.
- Split “computational advantage” into:
  1. `3V` state-count economy;
  2. sparse local plus one rank-one algebraic structure;
  3. measured runtime, for which no superiority is claimed.

### Delete

- Generic “low-cost” or “more efficient” wording lacking one of the three
  qualifiers above.
- Any prospective MITC4 comparison stated as completed evidence.

### Add

- Cite the final comparison table: Checkerboard uses 50.5--55.1% of the raw
  MidedgeTan scalar states over the reported mesh sequences while giving
  benchmark errors of comparable order.

## Section 2 — Checkerboard representation and discrete shell geometry

### Retain

- The physical edge-midpoint surface and exact raw gauge.
- Exact planarity of every Varignon B face for arbitrary raw configurations.
- W-centered curvature from neighboring B normals.
- The accepted moment boundary reconstruction and complete-interior policy for
  `X_W^rel`.

### Revise

- After the B/W definitions, add: “Exact planarity guarantees a unique local
  plane and normal; it does not imply distortion-independent approximation
  accuracy.”
- State that the displayed difference formulas, `gamma_W`, and the boundary
  moments use constant `hx,hy`; variable-step weights require a new derivation.

### Delete

- None of the mathematical definitions.

### Add

- A short implementation note or Appendix pointer: boundary closure is stored
  as at most a 3-by-3 indexed local stencil per W site.  The old dense operator
  was an implementation artifact, not part of the formulation.

## Section 3 — physical blind modes and variational observability

### Retain

- The quotient-space distinction between raw gauge and physical blind modes.
- The flat B-metric one-dimensional quotient defect.
- The two flat symmetric-W curvature defects.

### Revise

- Every use of “minimal completion” should read, or be immediately defined as:
  **“quotient-minimal only in the linear observation-rank sense required to
  remove the one-dimensional flat membrane quotient defect.”**
- Do not generalize the flat nullity theorem to arbitrary curved references.

### Delete

- Any unqualified suggestion that `X_M` is a uniquely determined continuum
  coefficient or a global minimality result across nonlinear formulations.

### Add

- The Round-3 flat n=16 in-plane operator diagnostic: nullity changes from 6
  to 5 when `c_beta` changes from 0 to 1; the one removed direction has base
  stiffness `9.62e-20` and completion-added stiffness `8.27e-5` in the stated
  normalization.  The remaining shared nullity contains rigid and raw-gauge
  directions.

## Section 4 — completion and finite-configuration energy

### Retain

- `beta_M=8mu` as the a priori trace-free-shear normalization paired with the
  chosen `X_M` normalization.
- `beta_W^0=mu gamma_W^0`, surface-polar transport, and the accepted energy.
- The explicit statement already present in the Appendix that `beta_M` is not
  uniquely dictated by the continuum Gram law.

### Revise

- “Neither coefficient depends on benchmark tuning” should be supported by
  two distinct facts: the analytical normalization and the Round-3 response
  plateau.  Do not replace this by “the coefficient is unimportant.”
- Describe the tangent exactly:
  `K_XM=2c[g g^T + X_M D2X_M]`, `c=2tmu` at `c_beta=1`.
- State that `D2X_M` is local and the only global term is `g g^T`.

### Delete

- Any reference-state-only approximation that omits `X_M D2X_M` when
  describing the finite-configuration implementation.

### Add

- Exact bordered solve and nonlinear verification: cantilever and hemisphere
  tangent relative differences `1.47e-16` and `1.37e-16`; the nonlinear local
  term is nonzero in both cases.
- Round-3 beta result: maximum relative displacement change over 0.25--4 is
  zero for the scalar plate branch, `1.63e-5` for roof, and `4.21e-10` for M3.
  Curved full-tangent nullity does not increase at `c_beta=0` under the tested
  BC; the exact additional defect is demonstrated only in the flat membrane
  operator.

## Numerical Results

### Retain

- Geometric convergence, recovered uniform plate/roof/M3 convergence, and the
  three nonlinear benchmark paths, with their existing limitations.
- MidedgeTan as the sole validated external baseline.
- Gauge-invariance figure and post-optimization values.

### Revise

- Replace every Round-1 reconstructed-driver displacement plot/table by the
  Round-3 final comparative evidence table and Figure B.
- Add the exact sparse/rank-one scaling figure.  Report n=128 feasibility and
  near-linear nnz growth, but explicitly note that the current Python/AD M3
  worker is slower and uses more RSS than the separate C++ baseline.
- Replace the Round-1 distortion plots by Round-2/3 Figure D.
- Replace the Round-1 beta plot by Round-3 Figure E.
- Use “converges on the tested sequences” rather than “converges robustly.”

### Delete

- Round-1 distortion, scaling, and beta figures based on reconstructed
  drivers.
- The nonuniform-spacing ratios as evidence about the current method; mark
  them outside the constant-spacing formulation scope.
- CPU or memory comparisons in the main evidence table.
- Any MITC4 row or placeholder interpreted as data.

### Add

- Final 12-row Checkerboard/MidedgeTan response, error, raw-state, and active-
  state table.
- Surface-preserving distortion limitation: at M3 n=12 and delta=.30,
  Checkerboard errors are 20.12% (smooth) and 28.89% (random), versus 30.22%
  and 43.74% for MidedgeTan.  Plate/roof variations remain modest.
- Beta sweep raw data and the flat quotient diagnostic.
- State that the plate manuscript driver is a scalar transverse branch, so its
  normal response is exactly independent of the membrane completion.

## Discussion

### Retain

- Position-only kinematics, common midpoint carrier, finite-rotation update,
  exact gauge, and competitive accuracy per state.
- The explicitly qualified fixed-`h` residual-locking discussion.

### Revise

- Replace a single “robustness” concept by:
  - **geometric well-posedness:** exact planar B patches and unique local frames
    on every nondegenerate sampled configuration;
  - **approximation robustness:** empirical and problem-dependent; mild for the
    tested plate/roof redistribution, poor for M3 at high redistribution.
- Replace a single “efficiency” concept by:
  - **state-count economy:** established, roughly half MidedgeTan raw states;
  - **algebraic structure:** established, local sparse plus one rank-one update;
  - **measured runtime:** no superiority established across different backends.
- Change “direct computational consequence” after state economy to “direct
  state-space consequence.”

### Delete

- “Computationally more efficient,” “low-cost,” or similar unqualified claims.
- “Locking-free.”  The available result is an aligned-cylinder approximation
  statement plus observed fixed-`h` residual locking that moves under
  refinement.
- Any statement that M3 compactness proves general robustness under mesh
  redistribution.

### Add

- The implementation bottleneck: sparse factorization and local AD Hessian,
  not the rank-one vector/update.
- The M3 redistribution weakness and constant-spacing scope as explicit
  limitations.
- MITC4 as future external validation, not missing data to be inferred.

## Conclusion

### Retain

- Exact midpoint geometry, physical quotient, blind-mode identification,
  position-only energy, and finite-configuration objectivity.

### Revise

- The strongest numerical conclusion should be “competitive continuum-
  reference accuracy per scalar state on the recovered benchmark sequences,”
  not “more efficient.”
- Describe completion as restoring the identified observations, not as proving
  universal stability or locking freedom.

### Delete

- Universal distorted-mesh robustness, CPU superiority, arbitrary-topology
  readiness, and locking-free conclusions.

### Add

- One closing limitation sentence: the current derivation is structured and
  constant-step, and M3 shows material sensitivity to spherical sample
  redistribution.
- One implementation sentence: an exact sparse-local plus rank-one nonlinear
  realization avoids hidden dense `X_M` coupling.

## Claim-by-claim disposition

| Claim | Disposition | Manuscript-safe wording |
|---|---|---|
| “robust under distorted meshes” | **Delete/rewrite** | “The local B geometry remains well-defined on all sampled nondegenerate meshes; displacement accuracy is problem-dependent and M3 is redistribution-sensitive.” |
| “computationally more efficient” | **Delete/rewrite** | “The method uses about half the raw scalar states of MidedgeTan on the reported sequences and has near-linear sparse/rank-one storage; no cross-backend runtime superiority is claimed.” |
| “low-cost” | **Use only with an object** | “The rank-one correction has low incremental storage/evaluation cost”; do not apply to the full solver. |
| “locking-free” | **Delete** | “The aligned generalized-cylinder branch has an `O(h^2)` near-isometry approximation, while fixed-`h` ultra-thin residual locking remains observable.” |
| “minimal completion” | **Strictly qualify** | “Quotient-minimal only in the linear observation-rank sense required to remove the one-dimensional flat membrane quotient defect.” |
| “position-only” | **Retain** | “The independent shell state consists of `3V` raw nodal position scalars; physical geometry is evaluated on the midpoint quotient.” |
| “competitive accuracy per state” | **Retain** | “On the recovered plate, roof, and M3 sequences, Checkerboard gives errors of comparable order using 50.5--55.1% of MidedgeTan raw scalar states.” |
| “exact planar local geometry” | **Retain with boundary** | “Every Varignon B face is exactly planar for arbitrary raw positions; this guarantees an unambiguous local frame, not distortion-insensitive approximation.” |
| “sparse/rank-one implementation” | **Retain/strengthen** | “The exact nonlinear tangent is a bounded-stencil sparse local part, including `2c X_M D2X_M`, plus one global rank-one update stored and solved through a border.” |

## Required manuscript figure/table routing

- Main text: Round-3 Figures A1/B/C/D and the final comparative evidence table.
- Optional main or Appendix: Figure E beta sensitivity.
- Appendix/Supplement: A2 quantitative gauge errors, boundary-storage figure,
  full scaling timing/RSS, curved/free-edge spectra, and all mesh-quality rows.
- Deprecated and excluded: the reconstructed-driver figure list in
  `03_final_figures/DEPRECATED_FIGURES.md`.
