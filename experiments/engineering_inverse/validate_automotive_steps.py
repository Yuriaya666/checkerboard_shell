"""Verify one-step and staged corotational equilibria at saved designs."""
from __future__ import annotations

import json
import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.mount_springs import mitc4_forward
from experiments.engineering_inverse.automotive_spring_design import OUT


def run():
    case=build_case(12,8)
    saved=json.loads((OUT/"automotive_spring_design_12x8.json").read_text())
    samples={"zero":np.zeros(8),
             "saved_cb_design":np.asarray(saved["adjustment_mm"])*.001}
    result={}
    for name,control in samples.items():
        a=mitc4_forward(case,control,corotational=True,analysis_steps=1,norm_type=0)
        b=mitc4_forward(case,control,corotational=True,analysis_steps=16,norm_type=0)
        if a["status"]!=0 or b["status"]!=0:
            raise RuntimeError(f"OpenSees step comparison did not converge: {name}")
        result[name]=dict(max_displacement_difference_m=float(np.max(np.abs(
            a["displacement"]-b["displacement"]))),
            max_mount_reaction_difference_N=float(np.max(np.abs(
                a["mount_reaction"]-b["mount_reaction"]))),
            one_step_seconds=a["model_construction_seconds"]+a["analysis_seconds"]+
                a["extraction_seconds"],
            sixteen_step_seconds=b["model_construction_seconds"]+b["analysis_seconds"]+
                b["extraction_seconds"],
            one_step_last_norm=a["last_convergence_norm"],
            sixteen_step_last_norm=b["last_convergence_norm"])
        if result[name]["max_displacement_difference_m"]>1e-7 or (
                result[name]["max_mount_reaction_difference_N"]>1e-2):
            raise AssertionError((name,result[name]))
    (OUT/"automotive_efficiency_step_equivalence.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    return result


if __name__=="__main__":run()
