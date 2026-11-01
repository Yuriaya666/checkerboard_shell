# Actual-units aircraft skin case

Run:

```bash
cd /home/song/projects/checkerboard_shell
PYTHONPATH=/tmp/openseespy_pkg:src:. MPLBACKEND=Agg \
  /home/song/miniconda3/envs/allatom/bin/python \
  -m experiments.cmame_round7_aircraft.aircraft_actual
PYTHONPATH=src:. MPLBACKEND=Agg \
  /home/song/miniconda3/envs/allatom/bin/python \
  -m experiments.cmame_round7_aircraft.make_actual_figures
```

Outputs:

- `raw_data/aircraft_actual.json`
- `figures/Figure_R7_aircraft_actual_deformation.png`
- `figures/Figure_R7_aircraft_actual_profile.png`
- `AIRCRAFT_ACTUAL_REPORT.md`
