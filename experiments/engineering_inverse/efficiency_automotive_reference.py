"""Common refined nonlinear acceptance for all complete optimizer outputs."""
from __future__ import annotations

import json
import time
import argparse

import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case,metrics,PARAMS
from experiments.engineering_inverse.efficiency_automotive import OUT,METHODS
from experiments.engineering_inverse.mount_springs import checkerboard_forward,mitc4_forward


def run(nx=24,ny=20,reference_steps=16,start="zero",methods=METHODS):
    case=build_case(nx,ny)
    path=OUT/f"automotive_efficiency_common_reference_{start}.json"
    report=(json.loads(path.read_text()) if path.exists() else
        dict(mesh=[nx,ny],reference_steps=reference_steps,
        reference="ASDShellQ4 corotational; same reference for all designs",
        start=start,
        methods={}))
    for method in methods:
        source=json.loads((OUT/f"automotive_efficiency_{method}_12x8_{start}.json").read_text())
        x=np.asarray(source["adjustment_mm"])
        control=x*.001
        t=time.perf_counter()
        reference=mitc4_forward(case,control,corotational=True,
                                  analysis_steps=reference_steps,norm_type=0,
                                  gauge_x_fix=True)
        reference_seconds=time.perf_counter()-t
        if reference["status"]!=0:raise RuntimeError(f"reference failed for {method}")
        q=reference["displacement"].ravel()
        e=case["geometry_error"]+case["O"]@q
        quality=metrics(case,q,e)
        error_mm=e*1000
        weights=np.tile([8.,1.],nx+1)
        reference_objective=float((weights*error_mm)@error_mm/len(error_mm)+.0005*x@x)
        force=float(np.max(np.linalg.norm(reference["mount_reaction"],axis=1)))
        accepted=bool(quality["max_abs_error_mm"]<=
            PARAMS["max_edge_tolerance_mm"]+1e-8 and force<=
            PARAMS["reaction_limit_N"]+1e-8 and np.max(abs(x))<=1.5+1e-8)
        t=time.perf_counter()
        cb=checkerboard_forward(case,control)
        cb_seconds=time.perf_counter()-t
        if cb["status"]!="CONVERGED":raise RuntimeError(f"refined CB failed for {method}")
        qcb=cb["displacement"].ravel()
        ecb=case["geometry_error"]+case["O"]@qcb
        report["methods"][method]=dict(adjustment_mm=x.tolist(),
            reference_status=int(reference["status"]),
            reference_seconds=reference_seconds,
            reference_system_size=reference["system_size"],
            reference_accepted=accepted,
            reference_objective=reference_objective,
            reference_quality=quality,
            reference_max_reaction_N=force,
            reference_force_balance_N=reference["force_balance"].tolist(),
            refined_cb_seconds=cb_seconds,refined_cb_quality=metrics(case,qcb,ecb),
            refined_cb_max_reaction_N=float(np.max(np.linalg.norm(cb["mount_reaction"],axis=1))),
            refined_cb_residual_inf=float(cb["residual_inf"]),
            max_displacement_difference_mm=float(np.max(np.linalg.norm(
                reference["displacement"]-cb["displacement"],axis=1))*1000))
        path.write_text(json.dumps(report,indent=2)+"\n")
        print(method,accepted,quality["max_abs_error_mm"],force,flush=True)
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--start",choices=("zero","cb_feasible"),default="zero")
    parser.add_argument("--methods",nargs="+",choices=METHODS,default=list(METHODS))
    args=parser.parse_args()
    run(start=args.start,methods=args.methods)
