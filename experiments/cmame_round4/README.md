# CMAME completion round (Round 4)

This round closes the experimental gaps identified in
`docs/CMAME_PROJECT_SYNTHESIS.md`.  It does not change the Checkerboard
formulation, the production energy, or any accepted coefficient.

## What is new

1. **Full 3V simply supported plate.** The recovered manuscript plate driver
   is the scalar transverse branch (only one vertical DOF per interior node).
   `common.full_plate()` solves the complete 3V production model
   (`paper_components`, `sparse_low_rank`) with `w=0` on all edges and three
   in-plane rigid pins.  It reports the raw 3V state count and the true active
   DOF count, so the comparison with MidedgeTan is on the same footing.
2. **One MidedgeTan adapter for every comparison row.** The official LibShell
   `MidedgeAngleTanFormulation` adapter is rebuilt from the pinned upstream
   commit and rerun on the same uniform logical grids for plate, roof, and M3.
   Raw and active states are measured from the adapter output, not inferred
   from a formula.
3. **Dual roof reference.** Every roof row reports errors against the
   historical `0.3024` and the modern `0.3005924566` reference.
4. **Nonlinear benchmarks on the current production energy.** The
   end-moment cantilever, alternating-force hemisphere, and open-ended
   cylinder pull-out drivers are ported to call
   `checkerboard_shell.mechanics.energy.paper_components`.  The archived
   kinematic reductions, loads, and reference curves are unchanged.
5. **Exact production-tangent checks** at the nonlinear final states of all
   three benchmarks;
6. **A documented MITC4 standard-shell baseline** via OpenSeesPy
   `ShellMITC4` (6 DOF/node), with the same meshes, loads, and physical BCs;
7. **M3 dual observables** — raw equator vertex (main, same as MidedgeTan)
   and physical-midpoint reconstruction (secondary).

## Reproduction

```bash
cd /home/song/projects/checkerboard_shell
PY=/home/song/miniconda3/envs/allatom/bin/python

# Rebuild the MidedgeTan adapter if needed (pinned LibShell commit)
LIBSHELL_SOURCE_ROOT=/tmp/libshell_dl/libshell_src \
EIGEN_SOURCE_ROOT=/tmp/eigen_dl/eigen-3.4.0 \
bash experiments/cmame_round2/05_baselines/build_midedgetan_distortion.sh \
     /tmp/checkerboard_round2_midedgetan_distortion

# 12-row unified comparison (plate n=8..64, roof n=12..32, M3 n=6..16)
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.unified_comparison

# Nonlinear benchmarks
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.cantilever --meshes 16 24 32
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.hemisphere --meshes 6 8 12 16
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.open_cylinder --meshes 8x8 12x12 16x16

# MITC4 baseline and combined three-method table
PYTHONPATH=/tmp/openseespy_pkg:src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.mitc4
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.combined_comparison

# Summaries and exact-tangent checks
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.plate_full
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.nonlinear_summary
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.nonlinear_figures
PYTHONPATH=src:. MPLBACKEND=Agg $PY -m experiments.cmame_round4.tangent_final_states
PYTHONPATH=src:. $PY -m experiments.cmame_round4.finalize
```

## Provenance

- LibShell: upstream `evouga/libshell`, commit
  `6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf`.
- Eigen: 3.4.0 headers.
- The rebuilt adapter's SHA-256 is recorded in
  `logs/environment.json`; the uniform values reproduce the Round-2
  historical gates to better than `1e-9` relative.
- Checkerboard: current workspace commit and uncommitted state are recorded in
  `logs/environment.json`.

## Outputs

- `UNIFIED_COMPARISON_TABLE.md` — 12-row Checkerboard/MidedgeTan table.
- `MITC4_TABLE.md` — OpenSeesPy `ShellMITC4` baseline.
- `COMBINED_COMPARISON_TABLE.md` — three-method table.
- `NONLINEAR_SUMMARY.md` — cantilever/hemisphere/open-cylinder metrics.
- `raw_data/unified_comparison.{csv,json}` — full long-form rows.
- `raw_data/tangent_final_states.json` — exact production-tangent checks.
- `figures/Figure_R4_full_plate.{png,pdf}` — full 3V vs scalar vs MidedgeTan.
- `figures/Figure_R4_unified_comparison.{png,pdf}` — error vs raw states.
- `figures/Figure_R4_nonlinear_paths.{png,pdf}` — three nonlinear paths.
- `figures/Figure_R4_combined_three_methods.{png,pdf}` — error vs raw states
  for all three methods.
- `raw_data/cantilever/`, `raw_data/hemisphere/`,
  `raw_data/open_cylinder/open_cylinder_pullout/` — per-mesh summaries and raw
  load paths.

## Caveats

- The M3 comparison still uses different discrete observables: Checkerboard
  reports a physical-midpoint moment reconstruction at the equator, while the
  LibShell driver reports the raw equator vertex.  Both target the same
  continuum value, but this should be stated explicitly in the paper.
- MITC4 is a Mindlin--Reissner element; its plate response at `t/a=0.01`
  contains a ~0.5% transverse-shear floor relative to the thin-plate
  reference.  A thin-limit check (`t=0.001`) reduces the error to 0.001%.
- The MITC4 comparison is mixed: MITC4 is more accurate on the plate
  (n=8--32) and roof at the same mesh, but uses 6 DOF/node; Checkerboard is
  better on M3 (n=12--16) and on the plate at n=64.  Do not claim uniform
  accuracy superiority over MITC4.
- MITC4 converges slowly on M3 with a raw corner point load; at `n=32` it
  reaches 0.47%, comparable to the Checkerboard raw-vertex value at `n=16`.
- MidedgeTan remains the only verified *midedge* baseline; MITC4 is the
  standard-shell baseline.
- The nonlinear path drivers use the archived dense reduced Newton solver, not
  the production sparse rank-one tangent.  The production tangent is checked
  independently at the archived final states of the cantilever (32x4),
  hemisphere (n=12), and open cylinder (8x8): the reduced
  `K_local + alpha g g^T` agrees with a dense AD Hessian to `1.6e-16`,
  `1.5e-15`, and `4.9e-15`, respectively (see
  `raw_data/tangent_final_states.json`).
