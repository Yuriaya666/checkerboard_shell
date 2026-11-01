# Experiment 5: positive/negative Gaussian curvature

`run_spectra.py --only curved` analyzes a spherical patch and a saddle at
n=8/16/32, with a two-raw-layer essential perimeter clamp.  The sampled
operators are J_B, J_B plus D X_M, J_bW, J_bW plus D X_W^rel, and the full
reference tangent.  `spectra.csv` stores the smallest 20 Hessian eigenvalues;
the observation Hessians are J-transpose-J Gram operators so their reported
singular values are the square roots of nonnegative eigenvalues.  `Table_K1.csv`
states **observed** near-nullity among the first 20, not exact total nullity.
Full mode shapes and geometry are in `raw_data/mode_shapes/*.npz`.

`Figure_K3` depicts stored softest full-tangent modes (independently scaled
for visibility).  The classification threshold and essential BC are recorded
in `config/experiment.json`.  A flat-state theorem is not inferred from these
two finite patches.

