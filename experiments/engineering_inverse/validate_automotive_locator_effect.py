"""Show why the additional x datum must also be applied to ASDShellQ4."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.mount_springs import mitc4_forward

OUT=Path(__file__).resolve().parent/"results"


def run():
    case=build_case(12,8)
    saved=json.loads((OUT/"automotive_efficiency_cb_consistent_12x8_zero.json").read_text())
    controls=np.asarray(saved["adjustment_mm"])*.001
    old=mitc4_forward(case,controls,corotational=True,analysis_steps=2,norm_type=0)
    matched=mitc4_forward(case,controls,corotational=True,analysis_steps=2,
                         norm_type=0,gauge_x_fix=True)
    if old["status"] or matched["status"]:
        raise RuntimeError((old["status"],matched["status"]))
    result=dict(input_adjustment_mm=saved["adjustment_mm"],
                max_displacement_change_m=float(np.max(abs(
                    old["displacement"]-matched["displacement"]))),
                max_bracket_reaction_change_N=float(np.max(abs(
                    old["mount_reaction"]-matched["mount_reaction"]))),
                matched_locator_reaction_N=float(matched["gauge_reaction_N"]),
                old_newton_iterations=int(old["newton_iterations_total"]),
                matched_newton_iterations=int(matched["newton_iterations_total"]))
    (OUT/"automotive_locator_effect_asd.json").write_text(
        json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
    return result


if __name__=="__main__":
    run()
