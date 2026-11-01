"""Cross-check the reference-feasible mounting candidate in nonlinear CB."""

import json
import numpy as np
import torch

from checkerboard_shell import BoundaryConditions,CheckerboardShell,NewtonSettings
from experiments.engineering_inverse.automotive_panel import (
    OUT,build_case,membrane_principal_strain,metrics,
    mount_reactions,prescribed_field)


def run():
    source=json.loads((OUT/"automotive_12x8.json").read_text())
    reference=json.loads((OUT/"automotive_reference.json").read_text())
    case=build_case()
    d=np.asarray(source["mitc4"]["adjustment_mm"])*1e-3
    mask=torch.zeros(case["mesh"].shape,dtype=torch.bool)
    mask.reshape(-1)[case["fixed"]]=True
    prescribed=prescribed_field(case,d)
    model=CheckerboardShell(case["mesh"],case["material"],
        BoundaryConditions(mask,prescribed),solver_settings=NewtonSettings(
            residual_tolerance=1e-6,maximum_newton_iterations=30))
    solved=model.solve()
    q=(solved.deformed_vertices-case["mesh"].reference_vertices).numpy().reshape(-1)
    error=case["geometry_error"]+case["O"]@q
    reaction=model.full_residual_and_reactions(solved.coordinates)[1].numpy().reshape(-1)
    mount_reaction=mount_reactions(case,reaction)
    data=dict(candidate_source="MITC4 constrained linear optimization; SLSQP did not declare convergence",
        adjustment_mm=(d*1e3).tolist(),
        checkerboard_nonlinear=dict(status=solved.status,iterations=solved.iterations,
            residual_inf=solved.residual_inf,
            **metrics(case,q,error),
            maximum_mount_reaction_N=float(np.max(np.linalg.norm(mount_reaction,axis=1))),
            mount_reaction_N=mount_reaction.tolist(),
            maximum_constraint_reaction_component_N=float(np.max(np.abs(reaction))),
            max_abs_membrane_principal_strain=float(np.max(membrane_principal_strain(case,q))),
            strain_energy_J=solved.energies["total"]),
        reference_mitc4_20x16=reference["meshes"]["20x16"]["mitc4"])
    (OUT/"automotive_crosscheck.json").write_text(json.dumps(data,indent=2)+"\n")
    return data


if __name__=="__main__":
    print(json.dumps(run(),indent=2))
