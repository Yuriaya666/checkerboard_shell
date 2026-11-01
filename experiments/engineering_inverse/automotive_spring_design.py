"""Direct nonlinear adjustment optimization with finite-stiffness brackets."""
from __future__ import annotations

import json
import time
from pathlib import Path
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch
from scipy.optimize import minimize

from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.engineering_inverse.automotive_panel import build_case, metrics, PARAMS
from experiments.engineering_inverse.mount_springs import (
    K_BRACKET, bracket_matrices, checkerboard_forward, mitc4_forward)

OUT=Path(__file__).resolve().parent/"results"


def response_sensitivity(case, solution, Kspring, D, W,
                         hessian_assembly_strategy="cached"):
    """Implicit derivatives of edge observations and bracket resultants."""
    n=3*case["mesh"].num_vertices
    fixed=np.array(sorted([3*node+c for node in case["datum"] for c in range(3)]
                          +[case["gauge_x_dof"]]))
    free=np.setdiff1d(np.arange(n),fixed)
    tangent,_,_=assemble_paper_low_rank_tangent(
        torch.as_tensor(solution["displacement"]),
        case["mesh"].reference_vertices,case["mesh"].hx,
        case["mesh"].hy,solution["model"].targets,
        case["material"],variant="CB-Full",
        colour_batch=7 if hessian_assembly_strategy=="legacy" else 21,
        assembly_strategy=hessian_assembly_strategy)
    reduced=tangent.__class__((tangent.local+Kspring)[free][:,free].tocsr(),
                              tangent.vector[free],tangent.coefficient)
    # The bordered solve is more reliable than a Sherman--Morrison update
    # when the local block is poorly conditioned: two large terms can cancel.
    bordered=reduced.augmented_matrix()
    rhs=np.vstack((D[free,:],np.zeros((1,D.shape[1]))))
    dq=np.zeros((n,8));dq[free]=spla.splu(bordered).solve(rhs)[:len(free)]
    dF=np.zeros((4,3,8))
    nodal=dq.reshape(-1,3,8)
    for p in range(4):
        for component in (1,2):
            index=2*p+component-1
            dF[p,component,index]=K_BRACKET
            dF[p,component]-=K_BRACKET*(W[p]@nodal[:,component,:])
    return case["O"]@dq,dF


def nonlinear_design(nx=12,ny=8, max_iterations=80, verbose=False):
    case=build_case(nx,ny)
    Kspring,D,W=bracket_matrices(case)
    observations=case["O"]
    weight=np.tile([8.0,1.0],case["mesh"].nx+1)
    cached={"key":None,"value":None}
    counts={"forward":0,"linearizations":0}
    history=[]

    def state(x):
        key=tuple(np.asarray(x,dtype=float))
        if cached["key"]==key:return cached["value"]
        d=np.asarray(x)*.001
        previous=(None if cached["value"] is None else
                  torch.as_tensor(cached["value"]["solution"]["displacement"])
                  .reshape(case["mesh"].shape))
        solution=checkerboard_forward(case,d,initial=previous)
        counts["forward"]+=1
        if solution["status"]!="CONVERGED":
            raise RuntimeError(f"CB equilibrium failed: {solution['status']} at {x}")
        q=solution["displacement"].reshape(-1)
        error=(case["geometry_error"]+observations@q)*1000
        J,dF=response_sensitivity(case,solution,Kspring,D,W)
        mount=solution["mount_reaction"]
        current=dict(solution=solution,error_mm=error,J=J,
                     reaction=mount,reaction_jac=dF,
                     objective=float((weight*error)@error/len(error)+.0005*np.dot(x,x)))
        counts["linearizations"]+=1
        cached.update(key=key,value=current)
        history.append(dict(adjustment_mm=np.asarray(x).tolist(),
            objective=current["objective"],
            max_edge_error_mm=float(np.max(np.abs(error))),
            max_mount_reaction_N=float(np.max(np.linalg.norm(mount,axis=1))),
            residual_inf=float(solution["residual_inf"]),
            continuation_steps=int(solution["continuation_steps"])))
        if verbose:
            print("state",counts["forward"],"objective",current["objective"],
                  "max error",np.max(abs(error)),"reaction",np.max(np.linalg.norm(mount,axis=1)),flush=True)
        return current

    def fun(x):return state(x)["objective"]
    def jac(x):
        data=state(x)
        return 2*data["J"].T@(weight*data["error_mm"])/len(weight)+.001*np.asarray(x)
    def reaction_constraint(x):
        return (PARAMS["reaction_limit_N"]-
                np.linalg.norm(state(x)["reaction"],axis=1))
    def reaction_jac(x):
        data=state(x)
        F=data["reaction"]
        return (-np.einsum("pc,pcj->pj",F,data["reaction_jac"])*.001
                /np.maximum(np.linalg.norm(F,axis=1)[:,None],1e-12))
    def shape_constraint(x):
        e=state(x)["error_mm"]
        limit=PARAMS["max_edge_tolerance_mm"]
        return np.r_[limit-e,limit+e]
    def shape_jac(x):
        J=state(x)["J"]
        return np.vstack((-J,J))
    start=time.perf_counter()
    constraints=[dict(type="ineq",fun=reaction_constraint,jac=reaction_jac),
                 dict(type="ineq",fun=shape_constraint,jac=shape_jac)]
    result=minimize(fun,np.zeros(8),jac=jac,method="SLSQP",
        bounds=[(-1.5,1.5)]*8,constraints=constraints,
        options=dict(ftol=1e-10,maxiter=max_iterations))
    try:final=state(result.x)
    except RuntimeError:final=None
    elapsed=time.perf_counter()-start
    payload=dict(mesh=[nx,ny],method="Checkerboard nonlinear + area-integrated bracket",
        bracket_stiffness_N_per_m=K_BRACKET,
        adjustment_mm=result.x.tolist(),optimizer_success=bool(result.success),
        optimizer_message=str(result.message),iterations=int(result.nit),
        counts=counts,optimization_seconds=elapsed,history=history,
        final=None if final is None else dict(**metrics(case,final["solution"]["displacement"].reshape(-1),
                                                final["error_mm"]*.001),
            residual_inf=final["solution"]["residual_inf"],
            mount_reaction_N=final["reaction"].tolist(),
            maximum_mount_reaction_N=float(np.max(np.linalg.norm(final["reaction"],axis=1))),
            force_balance_N=final["solution"]["force_balance"].tolist(),
            spring_energy_J=float(final["solution"]["spring_energy"]),
            shell_energy_J=float(final["solution"]["shell_energy"]),
            displacement_m=final["solution"]["displacement"].tolist(),
            objective=final["objective"]))
    (OUT/f"automotive_spring_design_{nx}x{ny}.json").write_text(json.dumps(payload,indent=2)+"\n")
    return payload


if __name__=="__main__":
    result=nonlinear_design(verbose=True)
    print(json.dumps({k:result[k] for k in ("optimizer_success","optimizer_message",
                                             "adjustment_mm","counts","final") if k!="final"},indent=2))
    if result["final"]:print("max error",result["final"]["max_abs_error_mm"],
                             "max reaction",result["final"]["maximum_mount_reaction_N"])
