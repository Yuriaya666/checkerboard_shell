# CMAME manuscript-closing Round 3

Round 3 contains no new benchmark and no shell-formulation change.  It uses
the recovered manuscript drivers to revalidate `beta_M`, builds the final
Checkerboard/MidedgeTan evidence table, collects the manuscript figure set,
and provides a section-by-section revision map.

Primary outputs:

- `CMAME_ROUND3_REPORT.md`
- `01_beta_revalidated/raw_data/beta_revalidated.{csv,json}`
- `02_final_comparison/MANUSCRIPT_COMPARISON_TABLE.{md,tex}`
- `03_final_figures/FIGURE_GUIDE.md`
- `04_manuscript_map/CMAME_ROUND3_MANUSCRIPT_MAP.md`

Reproduction:

```bash
python -m experiments.cmame_round3.run_beta_revalidated
python -m experiments.cmame_round3.build_final_comparison
python -m experiments.cmame_round3.collect_final_figures
python -m experiments.cmame_round3.finalize_round3
```
