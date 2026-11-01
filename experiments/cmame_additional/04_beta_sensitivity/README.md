# Experiment 4: membrane completion multiplier

`run_primary.py --only beta` sweeps c_beta = 0, .25, .5, 1, 2, 4, 8 for
M3 and roof at n=8/12.  Every case uses the same reference-tangent load,
geometry and physical constraints.  `beta_sensitivity.csv` and its JSON twin
record the observable, continuum-reference error, energy fraction, positive
spectrum, condition estimate, nullity, and quadratic stiffness in the
rank-one X_M direction.  `z_M_completion_added_stiffness` isolates the exactly
linear beta contribution; the *total* directional stiffness also includes
other observations and essential constraints.  Figures B1–B3 are backed by
these data and have PNG/PDF/SVG outputs.

