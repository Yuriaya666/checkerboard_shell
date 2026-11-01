"""Document multiple converged CB responses at the ASD-optimized controls."""
from __future__ import annotations

import json
import numpy as np
import torch

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.automotive_spring_design import OUT
from experiments.engineering_inverse.mount_springs import checkerboard_forward,mitc4_forward


def run():
    case=build_case(12,8)
    source=json.loads((OUT/"automotive_efficiency_asd_fd_12x8_zero_pre_common_locator.json").read_text())
    controls=np.asarray(source["adjustment_mm"])*.001
    asd=mitc4_forward(case,controls,corotational=True,
                      analysis_steps=2,norm_type=0)
    if asd["status"]!=0:raise RuntimeError(asd["status"])
    initial=torch.as_tensor(asd["displacement"]).reshape(case["mesh"].shape)
    result={"input_adjustment_mm":source["adjustment_mm"],"states":{}}
    for name,guess in (("reference_cold",None),("asd_displacement",initial)):
        cb=checkerboard_forward(case,controls,initial=guess)
        if cb["status"]!="CONVERGED":raise RuntimeError((name,cb["status"]))
        q=cb["displacement"].ravel()
        error=(case["geometry_error"]+case["O"]@q)*1000
        result["states"][name]=dict(status=cb["status"],
            residual_inf_N=float(cb["residual_inf"]),
            maximum_edge_error_mm=float(np.max(abs(error))),
            maximum_mount_reaction_N=float(np.max(np.linalg.norm(
                cb["mount_reaction"],axis=1))),
            displacement_m=cb["displacement"].tolist())
    (OUT/"automotive_efficiency_branch_probe.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print({k:{f:v for f,v in state.items() if f!="displacement_m"}
           for k,state in result["states"].items()})
    return result


if __name__=="__main__":run()
