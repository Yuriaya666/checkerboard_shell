"""Independent consistent-gradient check at saved CB optimization endpoints."""
from __future__ import annotations

import json
import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.automotive_spring_design import response_sensitivity,OUT
from experiments.engineering_inverse.mount_springs import checkerboard_forward,bracket_matrices


def run():
    case=build_case(12,8);K,D,W=bracket_matrices(case)
    weight=np.tile([8.,1.],13)
    output={}
    for start in ("zero","cb_feasible"):
        for method in ("cb_consistent","cb_fd"):
            path=OUT/f"automotive_efficiency_{method}_12x8_{start}.json"
            if not path.exists():continue
            r=json.loads(path.read_text())
            x=np.asarray(r["adjustment_mm"])
            state=checkerboard_forward(case,x*.001)
            if state["status"]!="CONVERGED":raise RuntimeError((start,method,state["status"]))
            J,_=response_sensitivity(case,state,K,D,W)
            e=(case["geometry_error"]+case["O"]@state["displacement"].ravel())*1000
            gradient=2*J.T@(weight*e)/len(e)+.001*x
            output[f"{start}_{method}"]=dict(gradient=gradient.tolist(),
                infinity_norm=float(np.max(abs(gradient))),
                min_distance_to_stroke_bound_mm=float(np.min(1.5-abs(x))),
                max_edge_error_mm=float(np.max(abs(e))),
                max_mount_reaction_N=float(np.max(np.linalg.norm(
                    state["mount_reaction"],axis=1))))
    (OUT/"automotive_efficiency_stationarity.json").write_text(
        json.dumps(output,indent=2)+"\n")
    print(json.dumps(output,indent=2))
    return output


if __name__=="__main__":run()
