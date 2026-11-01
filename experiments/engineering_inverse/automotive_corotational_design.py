"""Local finite-difference design refinement with independent nonlinear shell.

Starts from the separately obtained Checkerboard feasible point. Thus this is
a cross-method feasibility/refinement study; its runtime must not be compared
with a from-zero Checkerboard optimization as an equal-start speed contest.
"""
from __future__ import annotations

import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

from experiments.engineering_inverse.automotive_panel import build_case,metrics,PARAMS
from experiments.engineering_inverse.mount_springs import mitc4_forward

OUT=Path(__file__).resolve().parent/"results"


def run(trust_radius_mm=.2,step_m=1e-5):
    case=build_case(12,8)
    source=json.loads((OUT/"automotive_spring_design_12x8.json").read_text())
    x0=np.asarray(source["adjustment_mm"])
    started=time.perf_counter()
    forward=0
    def solve(x):
        nonlocal forward
        result=mitc4_forward(case,np.asarray(x)*1e-3,corotational=True)
        forward+=1
        if result["status"]!=0:raise RuntimeError(f"corotational reference failed at {x}")
        return result
    baseline=solve(x0)
    q0=baseline["displacement"].ravel()
    f0=baseline["mount_reaction"]
    J=np.zeros((case["O"].shape[0],8));F=np.zeros((4,3,8))
    for k in range(8):
        probe=x0.copy();probe[k]+=step_m*1e3
        r=solve(probe)
        J[:,k]=case["O"]@(r["displacement"].ravel()-q0)/step_m
        F[:,:,k]=(r["mount_reaction"]-f0)/step_m
        print("basis",k,flush=True)
    weights=np.tile([8.,1.],case["mesh"].nx+1)
    error0=(case["geometry_error"]+case["O"]@q0)*1e3
    def state(x):
        delta=np.asarray(x)-x0
        return error0+J@delta,f0+np.einsum("pcj,j->pc",F,delta*1e-3)
    def objective(x):
        e,_=state(x)
        return float((weights*e)@e/len(e)+.0005*np.dot(x,x))
    def gradient(x):
        e,_=state(x)
        return 2*J.T@(weights*e)/len(e)+.001*x
    def force_fun(x):return 300-np.linalg.norm(state(x)[1],axis=1)
    def force_jac(x):
        f=state(x)[1]
        return -1e-3*np.einsum("pc,pcj->pj",f,F)/np.maximum(
            np.linalg.norm(f,axis=1)[:,None],1e-12)
    def shape_fun(x):
        e,_=state(x);return np.r_[.6-e,.6+e]
    def shape_jac(x):return np.vstack((-J,J))
    bounds=[(max(-1.5,value-trust_radius_mm),
             min(1.5,value+trust_radius_mm)) for value in x0]
    optimum=minimize(objective,x0,jac=gradient,method="SLSQP",bounds=bounds,
        constraints=[dict(type="ineq",fun=force_fun,jac=force_jac),
                     dict(type="ineq",fun=shape_fun,jac=shape_jac)],
        options=dict(ftol=1e-10,maxiter=100))
    checked=solve(optimum.x)
    e=case["geometry_error"]+case["O"]@checked["displacement"].ravel()
    result=dict(method="ASDShellQ4 corotational local FD refinement",
        start="Checkerboard feasible optimum; not equal-start timing comparison",
        mesh=[12,8],step_m=step_m,trust_radius_mm=trust_radius_mm,
        adjustment_mm=optimum.x.tolist(),success=bool(optimum.success),
        optimizer_message=str(optimum.message),iterations=int(optimum.nit),
        analyses=forward,seconds=time.perf_counter()-started,
        baseline_objective=objective(x0),predicted_objective=objective(optimum.x),
        verified=dict(status=checked["status"],
            **metrics(case,checked["displacement"].ravel(),e),
            max_mount_reaction_N=float(np.max(np.linalg.norm(
                checked["mount_reaction"],axis=1))),
            mount_reaction_N=checked["mount_reaction"].tolist()))
    (OUT/"automotive_corotational_design.json").write_text(json.dumps(result,indent=2)+"\n")
    print("verified",result["verified"]["max_abs_error_mm"],
          result["verified"]["max_mount_reaction_N"],flush=True)
    return result


if __name__=="__main__":run()
