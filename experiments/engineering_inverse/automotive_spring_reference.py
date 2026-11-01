"""Fixed-design, common-input convergence audit for elastic panel brackets."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np

from experiments.engineering_inverse.automotive_panel import (
    build_case, metrics, membrane_principal_strain, PARAMS)
from experiments.engineering_inverse.mount_springs import (
    checkerboard_forward, mitc4_forward, K_BRACKET,
    bracket_outer_fiber_stress)

OUT=Path(__file__).resolve().parent/"results"


def _row(case,solution):
    if solution["status"] not in ("CONVERGED",0):
        raise RuntimeError(f"reference equilibrium did not converge: {solution['status']}")
    q=solution["displacement"].reshape(-1)
    error=case["geometry_error"]+case["O"]@q
    mount=solution["mount_reaction"]
    return dict(status=solution["status"],
        **metrics(case,q,error),
        mount_reaction_N=mount.tolist(),
        max_mount_reaction_N=float(np.max(np.linalg.norm(mount,axis=1))),
        max_flexure_root_stress_MPa=float(np.max(
            bracket_outer_fiber_stress(mount))/1e6),
        force_balance_N=solution["force_balance"].tolist(),
        spring_energy_J=float(solution["spring_energy"]),
        max_membrane_principal_strain=float(np.max(membrane_principal_strain(case,q))),
        min_gap_mm=float(np.min(PARAMS["target_gap_m"]-error[0::2])*1e3),
        residual_inf_N=solution.get("residual_inf"),
        shell_energy_J=solution.get("shell_energy"))


def run(meshes=((12,8),(16,12),(20,16),(24,20)),
        stiffness_factors=(.5,1.,2.)):
    design=json.loads((OUT/"automotive_spring_design_12x8.json").read_text())
    controls=np.array(design["adjustment_mm"])*1e-3
    report=dict(design_adjustment_mm=design["adjustment_mm"],
        stiffness_N_per_m=K_BRACKET,bracket_stiffness_factors=list(stiffness_factors),
        meshes={})
    for nx,ny in meshes:
        case=build_case(nx,ny)
        cb=checkerboard_forward(case,controls)
        mt=mitc4_forward(case,controls)
        asd_linear=mitc4_forward(case,controls,asd_linear=True)
        corot=mitc4_forward(case,controls,corotational=True)
        row=dict(checkerboard=_row(case,cb),mitc4_linear=_row(case,mt),
                 asd_linear=_row(case,asd_linear),
                 asd_corotational=_row(case,corot))
        report["meshes"][f"{nx}x{ny}"]=row
        print(nx,ny,[(k,v["max_abs_error_mm"],v["max_mount_reaction_N"],v["status"])
                     for k,v in row.items()],flush=True)
        (OUT/"automotive_spring_reference.json").write_text(json.dumps(report,indent=2)+"\n")
    case=build_case(12,8)
    report["stiffness_sensitivity"]={}
    for factor in stiffness_factors:
        k=K_BRACKET*factor
        cb=checkerboard_forward(case,controls,stiffness=k)
        mt=mitc4_forward(case,controls,stiffness=k)
        report["stiffness_sensitivity"][str(factor)]=dict(
            checkerboard=_row(case,cb),mitc4_linear=_row(case,mt))
    (OUT/"automotive_spring_reference.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":run()
