"""Cold-path finite-difference KKT check for the matched ASDShellQ4 endpoint."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case, PARAMS
from experiments.engineering_inverse.mount_springs import mitc4_forward

OUT=Path(__file__).resolve().parent/"results"


def run(steps=(.05,.01)):
    case=build_case(12,8)
    source=json.loads((OUT/"automotive_efficiency_asd_fd_12x8_zero.json").read_text())
    x=np.asarray(source["adjustment_mm"],dtype=float)
    weights=np.tile([8.,1.],13)
    def state(controls):
        result=None
        for increments in (2,4,8,16):
            result=mitc4_forward(case,controls*.001,corotational=True,
                                 analysis_steps=increments,norm_type=0,
                                 gauge_x_fix=True)
            if result["status"]==0:break
        if result["status"]!=0:
            raise RuntimeError((controls.tolist(),result["status"]))
        error=(case["geometry_error"]+case["O"]@result["displacement"].ravel())*1000
        force=np.linalg.norm(result["mount_reaction"],axis=1)
        objective=float((weights*error)@error/len(error)+.0005*controls@controls)
        return objective,float(np.max(abs(error))),float(np.max(force))
    objective,error,force=state(x)
    report=dict(controls_mm=x.tolist(),objective=objective,
                max_edge_error_mm=error,max_mount_reaction_N=force,
                min_bound_distance_mm=float(np.min(1.5-abs(x))),
                min_shape_margin_mm=float(PARAMS["max_edge_tolerance_mm"]-error),
                min_force_margin_N=float(PARAMS["reaction_limit_N"]-force),
                gradients={})
    for h in steps:
        gradient=np.zeros(8)
        for i in range(8):
            xp=x.copy();xm=x.copy();xp[i]+=h;xm[i]-=h
            gradient[i]=(state(xp)[0]-state(xm)[0])/(2*h)
        projected=gradient.copy()
        projected[x<=-1.5+1e-7]=np.minimum(gradient[x<=-1.5+1e-7],0.)
        projected[x>=1.5-1e-7]=np.maximum(gradient[x>=1.5-1e-7],0.)
        report["gradients"][str(h)]=dict(full=gradient.tolist(),
            projected_infinity_norm=float(np.max(abs(projected))))
        (OUT/"automotive_asd_stationarity.json").write_text(
            json.dumps(report,indent=2)+"\n")
        print(h,report["gradients"][str(h)],flush=True)
    return report


if __name__=="__main__":
    run()
