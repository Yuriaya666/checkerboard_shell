# Reproducing the positioning decision experiment

From the repository root, use a Python environment containing NumPy, SciPy,
PyTorch, Matplotlib, and OpenSeesPy. Install the package in editable mode or
set `PYTHONPATH=src:.` (also append any separate OpenSeesPy package directory).
Use the same interpreter for every command.

```bash
export PYTHONPATH=src:.
export OPENBLAS_NUM_THREADS=1 OMP_NUM_THREADS=1 MKL_NUM_THREADS=1
python -m experiments.cmame_positioning.diagnostic
python -m experiments.cmame_positioning.gauge_check
python -m experiments.cmame_positioning.run_reference
MPLBACKEND=Agg python -m experiments.cmame_positioning.plot_results
```

`diagnostic.py` writes the small n=8/12 fixed-input geometry and shell
comparisons. `run_reference.py` adds the specified mesh refinement and central
design-parameter difference quotients, saving results incrementally. It runs
many OpenSees solves and may take several minutes. `gauge_check.py` verifies
that the alternating raw-vertex gauge leaves the physical midpoint surface
and the accepted CB-Full energy unchanged. The plots read `summary.csv`.

All quantities are SI unless the plot axis explicitly says mm. The geometry
is a synthetic nonuniform-curvature diagnostic, **not** a reproduced object
from the Jiang et al. checkerboard design papers. ShellMITC4 is a linear shell
in these scripts; ASDShellQ4 is used in its default linear mode. This study
tests the claimed geometry-conversion necessity, not a nonlinear shell-element
accuracy ranking.

The full research decision and comparability limits are in [report.md](report.md).
