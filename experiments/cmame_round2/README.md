# CMAME additional experiments — Round 2

This directory contains the benchmark-driver recovery, sparse boundary
closure, exact nonlinear low-rank tangent, and the revalidated scaling and
surface-preserving distortion studies.  It does not overwrite or reinterpret
the first-round data in `experiments/cmame_additional/`.

The manuscript drivers were recovered from the clean historical repository
`/home/song/projects/CB-shell`, tag `v1.0-arxiv`, commit
`093ad89`.  Full source hashes and numerical rerun records are stored under
`00_driver_recovery/`.

Primary entry point: [`CMAME_ROUND2_REPORT.md`](CMAME_ROUND2_REPORT.md).

Reproduction commands (from the repository root, using the recorded Python
environment) are:

```bash
python -m experiments.cmame_round2.run_driver_recovery
python -m experiments.cmame_round2.run_structure_audits
python -m experiments.cmame_round2.run_boundary_equivalence
python -m experiments.cmame_round2.run_exactness
python -m experiments.cmame_round2.run_scaling
python -m experiments.cmame_round2.run_distortion
python -m experiments.cmame_round2.run_midedgetan_distortion
python -m experiments.cmame_round2.run_regressions
python -m experiments.cmame_round2.make_figures
python -m experiments.cmame_round2.finalize_artifacts
```

Scaling uses one excluded warm-up and three stored formal repeats per row.
The MidedgeTan distortion build expects the verified LibShell and Eigen source
roots recorded by `05_baselines/build_midedgetan_distortion.sh`; both can be
overridden with environment variables documented in that script.
