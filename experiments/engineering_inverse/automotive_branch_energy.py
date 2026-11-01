"""Independent branch energy and free-tangent stability audit at fixed controls."""
from __future__ import annotations

import json
from pathlib import Path
import numpy as np
import torch

from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.mount_springs import (
    K_BRACKET, bracket_matrices, checkerboard_forward)

OUT = Path(__file__).resolve().parent / "results"


def run():
    case = build_case(12, 8)
    Kspring, _, _ = bracket_matrices(case)
    source = json.loads((OUT / "automotive_efficiency_branch_probe.json").read_text())
    controls = np.asarray(source["input_adjustment_mm"]) * 0.001
    fixed = np.array(sorted([3*n+c for n in case["datum"] for c in range(3)]
                            +[case["gauge_x_dof"]]))
    free = np.setdiff1d(np.arange(3*case["mesh"].num_vertices), fixed)
    results = {}
    for name, record in source["states"].items():
        initial = torch.as_tensor(record["displacement_m"], dtype=torch.float64
                                  ).reshape(case["mesh"].shape)
        state = checkerboard_forward(case, controls, initial=initial)
        tangent, _, _ = assemble_paper_low_rank_tangent(
            torch.as_tensor(state["displacement"], dtype=torch.float64),
            case["mesh"].reference_vertices, case["mesh"].hx, case["mesh"].hy,
            state["model"].targets, case["material"])
        matrix = (tangent.local+Kspring)[free][:, free].tocsr()
        vector = tangent.vector[free]
        stiffness = matrix.toarray() + tangent.coefficient * np.outer(vector, vector)
        eigmin = float(np.linalg.eigvalsh(stiffness)[0])
        q = state["displacement"].ravel()
        edge = (case["geometry_error"]+case["O"]@q)*1000
        results[name] = dict(
            status=state["status"], residual_inf_N=float(state["residual_inf"]),
            total_internal_plus_bracket_energy_J=float(
                state["shell_energy"]+state["spring_energy"]),
            minimum_free_tangent_eigenvalue=eigmin,
            max_edge_error_mm=float(np.max(abs(edge))),
            max_bracket_reaction_N=float(np.max(
                np.linalg.norm(state["mount_reaction"], axis=1))))
    old = json.loads((OUT/"automotive_efficiency_branch_probe_pre_gauge.json").read_text())
    old_states=list(old["states"].values())
    old_q=[np.asarray(row["displacement_m"]).reshape(case["mesh"].shape)
           for row in old_states]
    reference=case["mesh"].reference_vertices
    targets=make_targets(reference,case["mesh"].hx,case["mesh"].hy)
    old_potentials=[]
    for q in old_q:
        energy=float(paper_components(
            reference+torch.as_tensor(q,dtype=torch.float64),reference,
            case["mesh"].hx,case["mesh"].hy,targets,case["material"])["total"])
        bracket=.5*sum(K_BRACKET*np.dot(state["weights"][p],
            (q.reshape(-1,3)[:,c]-controls[2*p+c-1])**2)
            for p in range(4) for c in (1,2))
        old_potentials.append(energy+bracket)
    difference=old_q[1]-old_q[0]
    results["old_ungauged_difference"]=dict(
        total_potential_J=old_potentials,
        potential_difference_J=float(old_potentials[1]-old_potentials[0]),
        maximum_x_difference_m=float(np.max(abs(difference[...,0]))),
        maximum_yz_difference_m=float(np.max(abs(difference[...,1:]))),
        checkerboard_odd_x_shift_m=float(difference[0,1,0]),
        checkerboard_even_x_shift_m=float(difference[0,0,0]))
    path = OUT / "automotive_efficiency_branch_energy.json"
    path.write_text(json.dumps(results, indent=2)+"\n")
    print(json.dumps(results, indent=2))
    return results


if __name__ == "__main__":
    run()
