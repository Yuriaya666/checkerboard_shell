"""Common independent MITC4 mesh check of aircraft stroke designs."""

import json
import time

import numpy as np

from experiments.engineering_inverse.aircraft_skin import OUT, build_case
from experiments.engineering_inverse.common import mitc4_prescribed_response


def evaluate(case,stations,strokes,mode):
    groups=[[3*n+2 for n in case["patches"][s]] for s in stations]
    fixed=sorted(set(case["spar"]+[d for g in groups for d in g]))
    start=time.perf_counter()
    response,q,reactions=mitc4_prescribed_response(case["mesh"],case["material"],
        fixed,groups,np.asarray(strokes),case["O"],case["load"],
        beam_lines=case["beam_lines"])
    elapsed=time.perf_counter()-start
    error=response-case["targets"][mode]
    lookup={d:k for k,d in enumerate(fixed)}
    net=[float(reactions[[lookup[d] for d in g]].sum()) for g in groups]
    return dict(rms_error_mm=float(np.sqrt(np.mean(error**2))*1e3),
                max_error_mm=float(np.max(np.abs(error))*1e3),
                reaction_N=net,analysis_seconds=elapsed,
                free_dofs=6*case["mesh"].num_vertices-len(fixed))


def run():
    source=json.loads((OUT/"aircraft_12x24.json").read_text())
    designs={}
    for method in ("checkerboard","mitc4"):
        winner=source[method]["winner"]
        if winner is not None:
            designs[method+"_linear"]=(winner["stations"],
                [mode["strokes_m"] for mode in winner["cases"]])
    final=OUT/"aircraft_nonlinear_final_12x24.json"
    if final.exists():
        result=json.loads(final.read_text())
        designs["checkerboard_nonlinear"]=(result["stations"],
            [mode["strokes_m"] for mode in result["cases"]])
    data=dict(reference="OpenSees ShellMITC4 + identical beam sections; linear kinematics",
              warning="This is not a geometrically nonlinear common reference.",
              meshes={})
    for nx,ny in ((12,24),(16,32),(20,40)):
        case=build_case(nx,ny)
        rows={}
        for name,(stations,inputs) in designs.items():
            rows[name]=[evaluate(case,stations,inputs[mode],mode)
                        for mode in range(2)]
        data["meshes"][f"{nx}x{ny}"]=rows
    (OUT/"aircraft_reference.json").write_text(json.dumps(data,indent=2)+"\n")
    return data


if __name__=="__main__":
    result=run()
    print(json.dumps(result["meshes"]["20x40"],indent=2))
