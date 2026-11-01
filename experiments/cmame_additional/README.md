# CMAME additional numerical experiments

This directory contains new, independent experiment drivers and **all recorded
outputs**, separate from the paper's original results.  The repository's
`src/checkerboard_shell` mathematics and default coefficients were not edited.

## Reproduction

Use an environment containing PyTorch, NumPy, SciPy, Matplotlib, pytest and
psutil.  The recorded executable was
`/home/song/miniconda3/envs/allatom/bin/python` (environment versions and git
commit in `environment.json`).  From the project root:

```bash
PYTHONPATH=src /home/song/miniconda3/envs/allatom/bin/python -m pytest -q tests experiments/cmame_additional/test_experiments.py
PYTHONPATH=src /home/song/miniconda3/envs/allatom/bin/python experiments/cmame_additional/run_primary.py --threads 4
PYTHONPATH=src /home/song/miniconda3/envs/allatom/bin/python experiments/cmame_additional/run_spectra.py --threads 4
PYTHONPATH=src /home/song/miniconda3/envs/allatom/bin/python experiments/cmame_additional/plot_modes.py
PYTHONPATH=src /home/song/miniconda3/envs/allatom/bin/python experiments/cmame_additional/finalize_artifacts.py
```

The spectral script can take 20 minutes or more.  `--only distortion|gauge|scaling|beta|consistency`
and `--only curved|boundary` rerun subsets, overwriting only the relevant new
experiment outputs.  Each experiment has `config/`, `raw_data/`, `figures/`,
`logs/`, `README.md`, and an `artifact_manifest.json`.  Plot files include
PNG, PDF and SVG.  Raw eigenmodes are compressed NumPy archives; curves and
tables have CSV/JSON backing data.  Root logs retain failed first attempts.

The validated MidedgeTan and MITC4 implementations presumed by the request
are absent here.  `baseline_adapter.py` is an explicit interface to add them;
no substitute, external published number, or incomparable CPU result is
represented as a computed baseline.  See [the numerical report](CMAME_NUMERICAL_REPORT.md)
before using any figure in a manuscript.

