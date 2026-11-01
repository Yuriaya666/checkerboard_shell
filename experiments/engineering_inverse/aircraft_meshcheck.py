"""Same-stroke nonlinear Checkerboard response across aircraft mesh sizes."""
from __future__ import annotations

import json
import numpy as np

from experiments.engineering_inverse.aircraft_skin import (
    OUT,build_case,nonlinear_optimize_controls)


def run(meshes=((12,24),(16,32),(20,40))):
    final=json.loads((OUT/"aircraft_nonlinear_final_12x24.json").read_text())
    selected=final["stations"]
    path=OUT/"aircraft_checkerboard_meshcheck.json"
    report=(json.loads(path.read_text()) if path.exists() else
            dict(stations=selected,meshes={}))
    for nx,ny in meshes:
        if f"{nx}x{ny}" in report["meshes"]:continue
        case=build_case(nx,ny)
        rows=[]
        for k,old in enumerate(final["cases"]):
            d=np.asarray(old["strokes_m"])
            result=nonlinear_optimize_controls(case,selected,d,
                case["targets"][k],max_iterations=0)
            result.pop("displacement_m",None)
            rows.append(result)
            print(nx,ny,k,result["rms_error_mm"],result["max_error_mm"],
                  result["reactions_N"],flush=True)
        report["meshes"][f"{nx}x{ny}"]=rows
        path.write_text(
            json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":run()
