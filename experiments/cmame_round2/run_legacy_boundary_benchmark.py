"""Measure the historical dense boundary path in isolated subprocesses."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
HISTORICAL = Path("/home/song/projects/CB-shell")


PROGRAM = r'''
import json, resource, time
import numpy as np, torch
from experiments.checkerboard_benchmark.mechanics.checkerboard_koiter import koiter_targets, _fast_boundary_geometry
n=int(__import__('sys').argv[1]); theta0=np.deg2rad(18.)
theta=torch.linspace(theta0,np.pi/2,n+1,dtype=torch.float64)
phi=torch.linspace(0,np.pi/2,n+1,dtype=torch.float64)
tt,pp=torch.meshgrid(theta,phi,indexing='ij')
q=10*torch.stack((torch.sin(tt)*torch.cos(pp),torch.sin(tt)*torch.sin(pp),torch.cos(tt)),-1)
hx=(np.pi/2)/n; hy=(np.pi/2-theta0)/n
started=time.perf_counter(); targets=koiter_targets(q,hx,hy); target_seconds=time.perf_counter()-started
times=[]
for _ in range(4):
 started=time.perf_counter(); g=_fast_boundary_geometry(q,hx,hy,targets); float(g.second_form.sum()); times.append(time.perf_counter()-started)
matrices=(targets.value_operator,targets.derivative_x_operator,targets.derivative_y_operator,targets.mixed_operator)
print(json.dumps({'n':n,'old_dense_bytes_measured':sum(x.numel()*x.element_size() for x in matrices),
'old_target_creation_seconds':target_seconds,'old_geometry_time_median_seconds':float(np.median(times[1:])),
'old_peak_rss_mb':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.}))
'''


def run():
    rows=[]
    environment=os.environ.copy();environment["PYTHONPATH"]=str(HISTORICAL)
    for n in (8,16,32,64):
        completed=subprocess.run([sys.executable,"-c",PROGRAM,str(n)],cwd=HISTORICAL,
            env=environment,text=True,capture_output=True)
        if completed.returncode==0:
            rows.append(json.loads(completed.stdout.strip().splitlines()[-1]))
        else:
            rows.append({"n":n,"status":"FAILED","reason":completed.stderr[-2000:]})
    rows.append({"n":128,"status":"NOT_ALLOCATED","old_dense_bytes_projected":8724676608,
                 "reason":"projection retained; allocating 8.72 GB solely to reconfirm the known shape was not justified"})
    out=ROOT/"01_sparse_boundary"/"raw_data";out.mkdir(parents=True,exist_ok=True)
    (out/"legacy_dense_boundary.json").write_text(json.dumps(rows,indent=2)+"\n")
    return rows


if __name__=="__main__":print(json.dumps(run(),indent=2))
