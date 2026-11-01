"""Repeat the complete 31-layout linear screen after factorization reuse."""
from __future__ import annotations

import json
import os
import platform
import time
from pathlib import Path
import numpy as np
import torch

from experiments.engineering_inverse.aircraft_skin import OUT,run


def measure(repeats=3,meshes=((12,24),(16,32),(20,40))):
    torch.set_num_threads(1)
    report=dict(environment=dict(machine=platform.machine(),
        processor=platform.processor(),python=platform.python_version(),
        openblas_threads=os.getenv("OPENBLAS_NUM_THREADS"),torch_threads=torch.get_num_threads()),
        timing_scope="model assembly + one spar-clamped factorization and compliance basis + complete 31-layout linear search; MITC4 complete 31-layout search with model rebuild; nonlinear checks excluded",
        before_single_runs_seconds={"12x24":3.62,"16x32":11.87,"20x40":39.92},
        before_source="prior RESEARCH_REPORT.md: 12x24 median of three, other grids single run",
        meshes={})
    for nx,ny in meshes:
        rows=[]
        for repeat in range(repeats):
            started=time.perf_counter()
            result=run(nx,ny,with_mitc4=True)
            cb=result["checkerboard"];mt=result["mitc4"]
            rows.append(dict(repeat=repeat,wall_seconds=time.perf_counter()-started,
                checkerboard_seconds=cb["total_seconds"],
                checkerboard_shell_assembly_seconds=cb["shell_assembly_seconds"],
                checkerboard_compliance_seconds=cb["compliance_preparation_seconds"],
                checkerboard_layout_seconds=cb["layout_search_seconds"],
                mitc4_seconds=mt["search_seconds"],
                mitc4_forward_analyses=mt["forward_analyses"],
                checkerboard_winner=(None if cb["winner"] is None else cb["winner"]["stations"]),
                mitc4_winner=(None if mt["winner"] is None else mt["winner"]["stations"])))
            print(nx,ny,repeat,rows[-1]["checkerboard_seconds"],
                  rows[-1]["mitc4_seconds"],flush=True)
            report["meshes"][f"{nx}x{ny}"]=dict(rows=rows,
                checkerboard_median_seconds=float(np.median([r["checkerboard_seconds"] for r in rows])),
                mitc4_median_seconds=float(np.median([r["mitc4_seconds"] for r in rows])))
            (OUT/"performance_round2.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":measure()
