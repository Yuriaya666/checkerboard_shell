# Experiment 1: mesh distortion and spacing

`run_primary.py --only distortion` evaluates full Checkerboard on plate,
Scordelis–Lo roof and M3 using the geometrically linear reference tangent and
a spectral pseudoinverse.  `config/experiment.json` records the complete grid;
`raw_data/all_cases.csv` includes every success, inverted mesh and unavailable
baseline, mesh quality, observables, active states, energy ratios and lowest
positive eigenvalues.  Random offsets use seed 20260921 and match the smooth
field's interior RMS component by component.  `raw_data/representative_*.npz`
contains actual reference/deformed nodal coordinates and connectivity.

Figures D1–D4 are backed by the CSV/NPZ data.  D3 contains only Checkerboard
because neither verified external baseline is installed.  For this study the
roof and M3 reference comparison assumes the geometry/physical BC/load
reconstruction documented in the root report; these configurations were **not**
present in the original project.  No nonlinear Newton convergence is claimed.

