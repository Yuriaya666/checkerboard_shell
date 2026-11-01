"""Cold-path gradient, reaction derivative, and KKT audit for saved CB designs."""
from __future__ import annotations

import argparse
import json
from pathlib import Path
import numpy as np
import torch

from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.engineering_inverse.automotive_panel import build_case, PARAMS
from experiments.engineering_inverse.automotive_spring_design import response_sensitivity
from experiments.engineering_inverse.mount_springs import bracket_matrices, checkerboard_forward

OUT = Path(__file__).resolve().parent / "results"


def run(source_method="cb_consistent", steps=(.05, .01), endpoint="zero"):
    case = build_case(12, 8)
    Kspring, D, W = bracket_matrices(case)
    weights = np.tile([8., 1.], 13)
    source = json.loads((OUT/f"automotive_efficiency_{source_method}_12x8_{endpoint}.json").read_text())
    x = np.asarray(source["adjustment_mm"], dtype=float)

    def evaluate(controls):
        state = checkerboard_forward(case, controls*.001)
        if state["status"] != "CONVERGED":
            raise RuntimeError((state["status"], controls.tolist()))
        q = state["displacement"].ravel()
        error = (case["geometry_error"]+case["O"]@q)*1000
        forces = state["mount_reaction"]
        objective = float((weights*error)@error/len(error)+.0005*controls@controls)
        return state, error, forces, objective

    state, error, forces, objective = evaluate(x)
    J, dF = response_sensitivity(case, state, Kspring, D, W)
    tangent, _, _ = assemble_paper_low_rank_tangent(
        torch.as_tensor(state["displacement"], dtype=torch.float64),
        case["mesh"].reference_vertices, case["mesh"].hx, case["mesh"].hy,
        state["model"].targets, case["material"])
    fixed=np.array(sorted([3*n+c for n in case["datum"] for c in range(3)]
                          +[case["gauge_x_dof"]]))
    free=np.setdiff1d(np.arange(3*case["mesh"].num_vertices),fixed)
    local=(tangent.local+Kspring)[free][:,free].toarray()
    vector=tangent.vector[free]
    min_eigenvalue=float(np.linalg.eigvalsh(
        local+tangent.coefficient*np.outer(vector,vector))[0])
    dF *= .001
    gradient = 2*J.T@(weights*error)/len(error)+.001*x
    reaction_gradient = -np.einsum("pc,pcj->pj", forces, dF)/np.linalg.norm(forces,axis=1)[:,None]
    constraint = np.r_[PARAMS["reaction_limit_N"]-np.linalg.norm(forces,axis=1),
                       PARAMS["max_edge_tolerance_mm"]-error,
                       PARAMS["max_edge_tolerance_mm"]+error]
    bound_distance = float(np.min(1.5-np.abs(x)))
    projected = gradient.copy()
    projected[x <= -1.5+1e-7] = np.minimum(gradient[x <= -1.5+1e-7],0.)
    projected[x >= 1.5-1e-7] = np.maximum(gradient[x >= 1.5-1e-7],0.)
    if np.min(constraint) <= 1e-5:
        raise RuntimeError("active nonlinear inequality needs multiplier audit")
    result = dict(source=source_method,controls_mm=x.tolist(),objective=objective,
                  gradient_per_mm=gradient.tolist(),gradient_infinity_norm=float(np.max(abs(gradient))),
                  min_inequality_margin=float(np.min(constraint)),
                  min_bound_distance_mm=bound_distance,
                  kkt_projected_gradient_infinity_norm=float(np.max(abs(projected))),
                  kkt_interpretation="nonlinear inequalities inactive; stroke bounds projected",
                  gauge_reaction_N=state["gauge_reaction_N"],
                  minimum_free_tangent_eigenvalue=min_eigenvalue,
                  max_edge_error_mm=float(np.max(abs(error))),
                  max_reaction_N=float(np.max(np.linalg.norm(forces,axis=1))),
                  finite_differences={})
    path = OUT/f"automotive_reliability_{source_method}_{endpoint}.json"
    for h in steps:
        fd_gradient = np.zeros(8)
        fd_shape = np.zeros_like(J)
        fd_force = np.zeros_like(dF)
        max_state_change = 0.
        for i in range(8):
            plus=x.copy();minus=x.copy();plus[i]+=h;minus[i]-=h
            sp,ep,Fp,fp=evaluate(plus)
            sm,em,Fm,fm=evaluate(minus)
            fd_gradient[i]=(fp-fm)/(2*h)
            fd_shape[:,i]=(ep-em)/(2*h)
            fd_force[:,:,i]=(Fp-Fm)/(2*h)
            max_state_change=max(max_state_change,
                float(np.linalg.norm(sp["displacement"]-state["displacement"])),
                float(np.linalg.norm(sm["displacement"]-state["displacement"])))
        result["finite_differences"][str(h)] = dict(
            objective_gradient=fd_gradient.tolist(),
            max_objective_gradient_difference=float(np.max(abs(fd_gradient-gradient))),
            max_shape_jacobian_difference=float(np.max(abs(fd_shape-J))),
            max_force_jacobian_difference_N_per_mm=float(np.max(abs(fd_force-dF))),
            max_state_change_m=max_state_change)
        path.write_text(json.dumps(result,indent=2)+"\n")
        print(source_method,h,result["finite_differences"][str(h)],flush=True)
    return result


if __name__ == "__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--method",choices=("cb_consistent","cb_fd"),
                        default="cb_consistent")
    parser.add_argument("--start",default="zero")
    parser.add_argument("--steps",nargs="+",type=float,default=[.05,.01])
    args=parser.parse_args()
    run(args.method,args.steps,args.start)
