# Experiment 6: free longitudinal strip edges

`run_spectra.py --only boundary` uses a flat and cylindrical 2:1 physical
strip, n=8/16/32 cells in each parameter direction.  Both longitudinal
boundaries are free except the 2×2 raw corner block required to remove all
mixed rigid-motion/raw-gauge combinations.  The first attempt with only
isolated point constraints produced exact zero modes which were traced to
that constraint deficiency; the failed attempt remains in the root spectral
log.  This corner constraint's physical support shrinks under refinement,
but it is stronger than the ideal minimum physical BC.

`spectra.csv` contains the first 20 eigenvalues of the bending observation,
its completed variant, and the full tangent, together with physical midpoint
norm, alternating score, and two-layer boundary-localization ratio.  Full
mode shapes are in `raw_data/mode_shapes/*.npz`; Figures E1–E3 have
PNG/PDF/SVG versions.  MidedgeTan/MITC4 spectra are explicitly marked
`UNAVAILABLE`; free-edge softness therefore cannot be attributed uniquely to
Checkerboard from this experiment.

