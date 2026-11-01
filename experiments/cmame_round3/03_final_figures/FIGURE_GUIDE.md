# Final manuscript figure guide

## Figure A1 — gauge and physical quotient

**Caption.** Raw Checkerboard configurations at `alpha=0` and `alpha=L` under
the alternating gauge transformation.  Although the raw vertices differ
strongly, the physical edge-midpoint representation overlaps to numerical
precision.

**Strictly supports.** The raw vertex representation contains a gauge orbit
whose physical midpoint geometry is unchanged.

**Does not support.** Accuracy, stability, or convergence under physical mesh
distortion; the large raw displacement is a gauge change, not a changed mesh.

## Figure A2 — finite-configuration gauge-invariance errors

**Caption.** Changes of physical geometry, internal observations, energy,
residual, and tangent under `q -> q + alpha (-1)^(i+j) A`, evaluated through
`alpha/L=1`.  Round-2 post-optimization regression gives maximum physical
midpoint motion `8.27e-15` and maximum relative tangent change `1.25e-15`.

**Strictly supports.** Finite-configuration internal gauge invariance is
realized near double-precision roundoff in the tested states.

**Does not support.** Gauge invariance of an arbitrary external raw-nodal load
potential, or exact zero relative error for a stress-free zero denominator.

## Figure B — recovered benchmark accuracy per raw state

**Caption.** Continuum-reference displacement error versus raw scalar state
count for Checkerboard and the recovered official LibShell MidedgeTan
implementation on the manuscript plate, Scordelis--Lo roof, and M3 sequences.

**Strictly supports.** Both methods provide convergent benchmark responses;
Checkerboard is competitive per raw state and uses approximately 50--55% of
the MidedgeTan raw scalar states on these grids.

**Does not support.** CPU-time or memory superiority, a MITC4 comparison, or a
universal accuracy ranking across other loads and meshes.

## Figure C — exact sparse/rank-one algebraic scaling

**Caption.** Augmented-system nonzeros versus raw scalar states for the exact
local sparse plus rank-one Checkerboard tangent.  The M3 rank-one correction
is represented by one border rather than a dense outer product.

**Strictly supports.** Near-linear sparse storage growth and the absence of an
explicit dense `g g^T` production matrix through n=128.

**Does not support.** A runtime or peak-memory advantage over the separate C++
LibShell backend; current Python/AD factorization and assembly remain costly.

## Figure D — surface-preserving distortion

**Caption.** Displacement error for Checkerboard and MidedgeTan under smooth
and fixed-seed random redistribution while every node remains on the same
plane, cylinder, or sphere and boundaries remain fixed.  All sampled cells
are valid.

**Strictly supports.** Stable sampled solves and modest response changes for
plate/roof, together with a significant M3 redistribution sensitivity for
both methods.  At delta=.30, Checkerboard M3 errors are 20.12% (smooth) and
28.89% (random).

**Does not support.** General distortion-insensitive accuracy, a universal
Varignon robustness advantage, or validity on arbitrary variable parameter
spacing outside the constant-`hx,hy` derivation.

## Figure E — revalidated beta_M sensitivity

**Caption.** Response change relative to `c_beta=1`, completion-energy fraction,
and completion-added `DX_M` stiffness for `beta_M=c_beta(8mu)`.  The plate
response is the manuscript scalar transverse branch, where the membrane
completion channel is inactive.

**Strictly supports.** `c_beta=1` lies on a broad normal-response plateau for
the three tested manuscript drivers; completion-added target-channel stiffness
is linear in `c_beta`.  A separate flat plate membrane-operator diagnostic
shows exactly one additional null direction at `c_beta=0`.

**Does not support.** Irrelevance of beta_M, parameter independence under
arbitrary nonlinear/in-plane loading, or a curved-reference theorem that
removing completion must create a new exact null mode.

## Main-text versus supplement recommendation

- Main text: A1 (optionally A2 as a paired subfigure), B, C, and D.
- Add E only if the completion parameter is questioned prominently; otherwise
  place it in the Appendix with the flat-operator nullity table.
- Supplement: full A2 metrics, timing/RSS plots, boundary-storage plot,
  curved/free-edge spectra, and complete distortion quality data.
