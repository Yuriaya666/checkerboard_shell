"""Mesh sensitivity and common MITC4 assessment of mounting adjustments."""

import json
import time

import numpy as np

from experiments.engineering_inverse.automotive_panel import (
    OUT,build_case,metrics,mount_reactions)
from experiments.engineering_inverse.common import mitc4_prescribed_response


def run():
    source=json.loads((OUT/"automotive_12x8.json").read_text())
    cb_candidate=source.get("checkerboard_nonlinear_feasible",
                            source["checkerboard"])
    inputs={"checkerboard":np.asarray(cb_candidate["adjustment_mm"])*1e-3,
            "mitc4":np.asarray(source["mitc4"]["adjustment_mm"])*1e-3}
    result=dict(reference="OpenSees ShellMITC4 linear shell at three meshes",
                warning="Not a geometrically nonlinear independent reference",
                meshes={})
    for nx,ny in ((12,8),(16,12),(20,16)):
        case=build_case(nx,ny)
        entries={}
        for name,d in inputs.items():
            start=time.perf_counter()
            response,q,reactions=mitc4_prescribed_response(case["mesh"],
                case["material"],case["fixed"],case["design"],d,case["O"])
            error=case["geometry_error"]+response
            lookup={int(d):k for k,d in enumerate(case["fixed"])}
            mount=mount_reactions(case,reactions,lookup)
            entries[name]=dict(**metrics(case,q,error),
                maximum_mount_reaction_N=float(np.max(np.linalg.norm(mount,axis=1))),
                maximum_constraint_reaction_component_N=float(np.max(np.abs(reactions))),
                analysis_seconds=time.perf_counter()-start,
                free_dofs=6*case["mesh"].num_vertices-len(case["fixed"]))
        result["meshes"][f"{nx}x{ny}"]=entries
    (OUT/"automotive_reference.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__=="__main__":
    print(json.dumps(run()["meshes"]["20x16"],indent=2))
