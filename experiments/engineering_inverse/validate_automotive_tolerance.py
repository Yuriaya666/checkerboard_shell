"""Sensitivity of ASD forward timing and outputs to residual tolerance.

This is diagnostic only. The published complete optimization retains the
strict 1e-6 maximum-component unbalance criterion.
"""
from __future__ import annotations

import json
import time
import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.automotive_spring_design import OUT
from experiments.engineering_inverse.mount_springs import mitc4_forward


def run():
    saved=json.loads((OUT/"automotive_spring_design_12x8.json").read_text())
    control=np.asarray(saved["adjustment_mm"])*.001
    case=build_case(12,8)
    output={}
    for tolerance in (1e-6,1e-5,1e-4):
        t=time.perf_counter()
        r=mitc4_forward(case,control,corotational=True,analysis_steps=2,
                        norm_type=0,residual_tolerance=tolerance)
        output[str(tolerance)]=dict(status=r["status"],seconds=time.perf_counter()-t,
            newton_iterations=r["newton_iterations_total"],
            displacement_m=r["displacement"].tolist(),
            reaction_N=r["mount_reaction"].tolist())
        if r["status"]!=0:raise RuntimeError((tolerance,r["status"]))
    base=output[str(1e-6)]
    for tolerance in (1e-5,1e-4):
        item=output[str(tolerance)]
        item["max_displacement_difference_m"]=float(np.max(np.abs(
            np.asarray(item["displacement_m"])-base["displacement_m"])))
        item["max_reaction_difference_N"]=float(np.max(np.abs(
            np.asarray(item["reaction_N"])-base["reaction_N"])))
    (OUT/"automotive_efficiency_tolerance_sensitivity.json").write_text(
        json.dumps(output,indent=2)+"\n")
    print({key:{k:v for k,v in row.items() if k not in
        ("displacement_m","reaction_N")} for key,row in output.items()})
    return output


if __name__=="__main__":run()
