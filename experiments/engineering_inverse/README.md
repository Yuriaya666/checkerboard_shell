# Engineering inverse-design studies

The final exact-tangent and nonlinear-reliability audit is documented in
[tangent_reliability_report.md](tangent_reliability_report.md). It supersedes
the automotive endpoint and branch conclusions in the earlier
[strict efficiency report](efficiency_report.md). The final comparison uses
the same additional x locator in Checkerboard and ASDShellQ4, and all three
automotive designs pass the stated 24×20 reference limits. Earlier studies
remain in [report.md](report.md) and
[RESEARCH_REPORT.md](RESEARCH_REPORT.md); their older model conclusions should
not be cited as the final efficiency result.

```bash
python -m pip install -e . openseespy pytest
python -m experiments.engineering_inverse.run_efficiency
```

The command above runs the standard three-method workflow. For the
final matched-locator and old/new Hessian comparison, use the commands in
`tangent_reliability_report.md`; the complete automotive runs take hours on
the reported CPU. Parameters are in `config.json`, raw records and CSV tables
in `results/`, and vector PDF/SVG figures in `figures/`.
The final figures and tables alone can be rebuilt with
`python -m experiments.engineering_inverse.make_tangent_reliability_figures`.
