"""Audit all lower-cardinality aircraft layouts against nonlinear equilibrium."""

import json
import time

import numpy as np

from experiments.engineering_inverse.aircraft_skin import (
    OUT, build_case, layout_operator, nonlinear_optimize_controls,
    prepare_layout_compliance)
from experiments.engineering_inverse.common import checkerboard_displacement_operator


def run(nx=12,ny=24,max_count=2):
    source=json.loads((OUT/f"aircraft_{nx}x{ny}.json").read_text())
    case=build_case(nx,ny)
    base=checkerboard_displacement_operator(case["mesh"],case["material"],
        case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"])
    prepared=prepare_layout_compliance(case,base)
    path=OUT/f"aircraft_nonlinear_audit_{nx}x{ny}.json"
    report=(json.loads(path.read_text()) if path.exists() else
        dict(mesh=[nx,ny],max_count=max_count,
             reference="same coupled Checkerboard shell + beam line energy",
             statuses=[]))
    report["status_semantics"]=(
        "not_accepted_local means the tested local nonlinear design failed the "
        "engineering tolerances; it is not a global infeasibility certificate")
    for row in report["statuses"]:
        if row["status"]=="infeasible":row["status"]="not_accepted_local"
    path.write_text(json.dumps(report,indent=2)+"\n")
    start=time.perf_counter()
    entries=[r for r in source["checkerboard"]["search"]
             if r["count"]<=max_count]
    for entry in entries:
        selected=entry["stations"]
        if any(row["stations"]==selected for row in report["statuses"]):
            continue
        op=layout_operator(case,selected,prepared)
        cases=[None,None]
        # The span-gradient case is more demanding. A failure rules out layout.
        try:
            k=1
            d=np.asarray(entry["cases"][k]["strokes_m"])
            initial=op["q0"]+op["Q"]@d
            cases[k]=nonlinear_optimize_controls(case,selected,d,
                case["targets"][k],max_iterations=4,
                initial_displacement=initial)
            if cases[k]["success"]:
                k=0
                d=np.asarray(entry["cases"][k]["strokes_m"])
                initial=op["q0"]+op["Q"]@d
                cases[k]=nonlinear_optimize_controls(case,selected,d,
                    case["targets"][k],max_iterations=4,
                    initial_displacement=initial)
            status=("feasible" if all(c is not None and c["success"] for c in cases)
                    else "not_accepted_local")
        except Exception as exc:
            status="solver_failure"
            cases[k]=dict(error=f"{type(exc).__name__}: {exc}")
        report["statuses"].append(dict(stations=selected,count=len(selected),
                                        status=status,cases=cases))
        report["elapsed_seconds"]=time.perf_counter()-start
        path.write_text(json.dumps(report,indent=2)+"\n")
        print(selected,status,report["elapsed_seconds"],flush=True)
    return report


if __name__=="__main__":
    run()
