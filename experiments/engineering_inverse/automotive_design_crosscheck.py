"""Re-evaluate both saved automotive designs on the same fine shell models."""
from __future__ import annotations

import json
import numpy as np

from experiments.engineering_inverse.automotive_panel import build_case,metrics
from experiments.engineering_inverse.mount_springs import checkerboard_forward,mitc4_forward
from experiments.engineering_inverse.aircraft_skin import OUT


def run(meshes=((12,8),(20,16),(24,20))):
    designs={
        "checkerboard_nonlinear":json.loads((OUT/"automotive_spring_design_12x8.json").read_text())["adjustment_mm"],
        "asd_corotational_refined":json.loads((OUT/"automotive_corotational_design.json").read_text())["adjustment_mm"]}
    report=dict(designs=designs,meshes={})
    for nx,ny in meshes:
        case=build_case(nx,ny)
        path=OUT/"automotive_design_crosscheck.json"
        mesh=report["meshes"].get(f"{nx}x{ny}",{})
        for design_name,controls_mm in designs.items():
            if design_name in mesh:continue
            controls=np.asarray(controls_mm)*1e-3
            methods={}
            for name,solver,kwargs in (
                ("checkerboard",checkerboard_forward,{}),
                ("shellmitc4_linear",mitc4_forward,{}),
                ("asd_corotational",mitc4_forward,{"corotational":True})):
                solution=solver(case,controls,**kwargs)
                if solution["status"] not in ("CONVERGED",0):
                    raise RuntimeError(f"{name} failed on {nx}x{ny}: {solution['status']}")
                q=solution["displacement"].ravel()
                error=case["geometry_error"]+case["O"]@q
                methods[name]=dict(status=solution["status"],**metrics(case,q,error),
                    max_mount_reaction_N=float(np.max(np.linalg.norm(
                        solution["mount_reaction"],axis=1))),
                    force_balance_N=solution["force_balance"].tolist(),
                    current_moment_balance_Nm=solution["current_moment_balance"].tolist())
            mesh[design_name]=methods
            print(nx,ny,design_name,
                  [(k,v["max_abs_error_mm"],v["max_mount_reaction_N"])
                   for k,v in methods.items()],flush=True)
        report["meshes"][f"{nx}x{ny}"]=mesh
        path.write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":run()
