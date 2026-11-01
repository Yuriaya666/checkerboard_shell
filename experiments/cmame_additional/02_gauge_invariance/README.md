# Experiment 2: finite-configuration raw gauge

`run_primary.py --only gauge` tests curved reference M3, linearly deformed M3
and a converged nonlinear cantilever for seven values of alpha/L.  `config/`
records the non-axial direction and amplitudes.  `raw_data/state_*.npz` saves
all raw reference/current positions, and `invariance_errors.csv` saves the
unrounded changes of edge midpoints, B metric, W curvature, completions,
internal energy, internal residual and tangent.

Relative internal energy/residual changes are undefined (`NaN`) when their
unshifted value is effectively zero.  Absolute differences and a tangent ×
characteristic-length normalized residual are retained instead.  External
dead-load potential is intentionally excluded; an arbitrary raw-node load is
not gauge invariant.  The selected benchmark loads are tested separately for
gauge annihilation.  Figures G1 and G2 have PNG/PDF/SVG versions.

