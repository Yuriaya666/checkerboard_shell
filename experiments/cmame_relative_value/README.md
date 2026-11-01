# Relative numerical value study

This directory adds no inverse-design problem and does not alter the accepted CB-Full energy. It reuses the recovered plate, roof, M3, official LibShell MidedgeTan, and OpenSees ShellMITC4 drivers. `evidence_inventory.md` separates prior analytical results from new comparisons.

## Environment

Use a Python with NumPy, SciPy, PyTorch, Matplotlib, OpenSeesPy and pytest. Set `PYTHONPATH=src:.` from the repository root. If OpenSeesPy is installed in a separate directory, append that directory to `PYTHONPATH`. Set `OPENBLAS_NUM_THREADS=OMP_NUM_THREADS=MKL_NUM_THREADS=1` for comparable CPU observations.

Build the official LibShell adapter using the pinned LibShell/Eigen source trees recorded in `experiments/cmame_round2/05_baselines/midedgetan_distortion_provenance.json`:

```sh
LIBSHELL_SOURCE_ROOT=/path/to/libshell EIGEN_SOURCE_ROOT=/path/to/eigen \
  bash experiments/cmame_round2/05_baselines/build_midedgetan_distortion.sh /path/to/midedgetan_thickness
export MIDEDGETAN_THICKNESS_BIN=/path/to/midedgetan_thickness
```

The optional fifth CLI argument adds thickness selection to the existing adapter; omitting it retains every original benchmark thickness and output. The adapter has no benchmark-specific retuning.

Run from the repository root:

```sh
python -m experiments.cmame_relative_value.run_scan
python -m experiments.cmame_relative_value.run_scan --references --methods MidedgeTan MITC4
python -m experiments.cmame_relative_value.extend_thin_reference
python -m experiments.cmame_relative_value.run_completion
python -m experiments.cmame_relative_value.analyze
python -m experiments.cmame_relative_value.record_environment
python -m pytest -q tests
```

Each scan writes CSV/JSON after each case. The main scan uses four meshes per curved benchmark and five thicknesses. Reference checks use 48, 64 and 96 cells per direction for both external methods; the thinnest M3 case extends to 128. All M3 methods receive the same raw-node point loads and use the same raw-node output. Existing CB M3 constraints still act on physical midpoints while the external codes use their own symmetry-rotation constraints; this kinematic difference limits attribution.

The roof CB and MidedgeTan drivers integrate pressure over the parameterized cylinder area, while the existing MITC4 driver uses the chordal element area. `analyze.py` applies the exact uniform-cylinder chord-to-arc factor to MITC4 linear responses before cross-method comparisons. The untouched solver outputs remain in the raw CSV.

`thickness_mesh.svg` reports *scalar-response discrepancy to the finest MITC4 surrogate*, not certified continuum error. `reference_quality.csv` gives 64→96 refinement change, final refinement change and cross-method disagreement. Where reference quality is provisional, conclusions are restricted to fixed-grid response trends. The energy fractions in `completion_ablation.csv` and `thickness_scan.csv` sample the reference-tangent quadratic form along each linear response; `cylinder_isometry.csv` is a prescribed deformation rather than a load-equilibrium test.
