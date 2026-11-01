"""Fixed-input forward audit of the existing automotive hard-patch model."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import torch

from checkerboard_shell import BoundaryConditions, CheckerboardShell, NewtonSettings
from experiments.engineering_inverse.automotive_panel import (
    build_case, mount_reactions, prescribed_field, metrics, PARAMS)
from experiments.engineering_inverse.common import (
    checkerboard_displacement_operator, mitc4_prescribed_response)

OUT = Path(__file__).resolve().parent / "results"


def _moment(points, forces):
    return np.cross(points, forces).sum(axis=0)


def run(meshes=((12,8),(16,12),(20,16)), amplitude=1e-4):
    report = {"connection":"all patch nodes fixed in x, same prescribed y/z translation",
              "amplitude_m":amplitude, "meshes":{}}
    for nx,ny in meshes:
        case=build_case(nx,ny)
        xyz=case["mesh"].reference_vertices.detach().numpy().reshape(-1,3)
        inputs={"zero":[0]*8,
                "one_y":[amplitude,0,0,0,0,0,0,0],
                "one_z":[0,amplitude,0,0,0,0,0,0],
                "opposite_z":[0,amplitude,0,-amplitude,0,0,0,0],
                "four_z":[0,amplitude,0,amplitude,0,amplitude,0,amplitude]}
        existing=OUT/"automotive_12x8.json"
        if existing.exists():
            payload=json.loads(existing.read_text())
            inputs["old_optimized"]=(np.array(payload["checkerboard"]["adjustment_mm"])*1e-3).tolist()
        model=CheckerboardShell(case["mesh"],case["material"],
             BoundaryConditions(torch.zeros(case["mesh"].shape,dtype=torch.bool)),
             solver_settings=NewtonSettings(residual_tolerance=1e-5,
                                            maximum_newton_iterations=30))
        row={"mount_node_counts":[len(p) for p in case["mount_patches"]],
             "fixed_dofs":len(case["fixed"]), "cases":{}}
        for name,values in inputs.items():
            mask=torch.zeros(case["mesh"].shape,dtype=torch.bool)
            mask.reshape(-1)[case["fixed"]]=True
            model.boundary_conditions=BoundaryConditions(mask,prescribed_field(case,values))
            solved=model.solve()
            q=(solved.deformed_vertices-case["mesh"].reference_vertices).detach().numpy().reshape(-1,3)
            full_reaction=model.full_residual_and_reactions(solved.coordinates)[1].numpy().reshape(-1,3)
            cb_mount=mount_reactions(case,full_reaction)
            cb_datum=full_reaction[case["datum"]]
            cb_balance=full_reaction.sum(axis=0)
            cb_moment=_moment(xyz,full_reaction)
            lookup={int(d):i for i,d in enumerate(case["fixed"])}
            mt_response,mt_q,mt_reduced=mitc4_prescribed_response(case["mesh"],
                case["material"],case["fixed"],case["design"],values,case["O"])
            mt_full=np.zeros((case["mesh"].num_vertices,3))
            for dof,value in zip(case["fixed"],mt_reduced):mt_full[dof//3,dof%3]=value
            mt_mount=mount_reactions(case,mt_reduced,lookup)
            mt_datum=mt_full[case["datum"]]
            mt_balance=mt_full.sum(axis=0)
            mt_moment=_moment(xyz,mt_full)
            d=np.asarray(values)
            cb_energy=.5*sum(d[2*k+c]*cb_mount[k,c+1] for k in range(4) for c in (0,1))
            mt_energy=.5*sum(d[2*k+c]*mt_mount[k,c+1] for k in range(4) for c in (0,1))
            error_cb=case["geometry_error"]+case["O"]@q.ravel()
            error_mt=case["geometry_error"]+mt_response
            row["cases"][name]={
                "input_m":values,
                "checkerboard":{"status":solved.status,"residual_inf":solved.residual_inf,
                    "mount_reactions_N":cb_mount.tolist(),"datum_reactions_N":cb_datum.tolist(),
                    "force_balance_N":cb_balance.tolist(),"moment_balance_Nm":cb_moment.tolist(),
                    "energy_J":cb_energy,"internal_energy_J":solved.energies["total"],
                    "max_edge_error_mm":metrics(case,q.ravel(),error_cb)["max_abs_error_mm"],
                    "max_displacement_mm":float(np.max(np.linalg.norm(q,axis=1))*1e3),
                    "displacement_m":q.tolist()},
                "mitc4":{"mount_reactions_N":mt_mount.tolist(),"datum_reactions_N":mt_datum.tolist(),
                    "force_balance_N":mt_balance.tolist(),"moment_balance_translations_Nm":mt_moment.tolist(),
                    "energy_from_work_J":mt_energy,
                    "max_edge_error_mm":metrics(case,mt_q,error_mt)["max_abs_error_mm"],
                    "max_displacement_mm":float(np.max(np.linalg.norm(mt_q.reshape(-1,3),axis=1))*1e3),
                    "displacement_m":mt_q.reshape(-1,3).tolist()}}
            print(f"{nx}x{ny} {name}: |F| CB {np.max(np.linalg.norm(cb_mount,axis=1)):.3g} N; MITC4 {np.max(np.linalg.norm(mt_mount,axis=1)):.3g} N",flush=True)
        report["meshes"][f"{nx}x{ny}"]=row
        (OUT/"automotive_fixed_input_diagnostics.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":run()
