"""ShellMITC4 finite-difference basis design with the same elastic brackets."""
from __future__ import annotations

import json
import time
from pathlib import Path
import numpy as np
from scipy.optimize import minimize

from experiments.engineering_inverse.automotive_panel import build_case,metrics,PARAMS
from experiments.engineering_inverse.mount_springs import mitc4_forward

OUT=Path(__file__).resolve().parent/"results"


def run(nx=12,ny=8):
    case=build_case(nx,ny)
    start=time.perf_counter()
    zero=mitc4_forward(case,np.zeros(8))
    if zero["status"]!=0:raise RuntimeError("zero-input MITC4 failed")
    step=1e-4
    q_basis=[];reaction_basis=[]
    for k in range(8):
        d=np.zeros(8);d[k]=step
        s=mitc4_forward(case,d)
        if s["status"]!=0:raise RuntimeError(f"MITC4 basis {k} failed")
        q_basis.append((s["displacement"].ravel()-zero["displacement"].ravel())/step)
        reaction_basis.append((s["mount_reaction"]-zero["mount_reaction"])/step)
    Q=np.column_stack(q_basis)
    F=np.stack(reaction_basis,axis=-1)
    basis_seconds=time.perf_counter()-start
    edge=case["geometry_error"]+case["O"]@zero["displacement"].ravel()
    J=case["O"]@Q
    weights=np.tile([8.,1.],nx+1)
    scale=.001
    def data(x):
        d=x*scale
        return (edge+J@d)/scale,zero["mount_reaction"]+np.einsum("pcj,j->pc",F,d)
    def fun(x):
        error,_=data(x)
        return float((weights*error)@error/len(error)+.0005*np.dot(x,x))
    def jac(x):
        error,_=data(x)
        return 2*J.T@(weights*error)/len(error)+.001*x
    def force_con(x):return PARAMS["reaction_limit_N"]-np.linalg.norm(data(x)[1],axis=1)
    def force_jac(x):
        force=data(x)[1]
        return -.001*np.einsum("pc,pcj->pj",force,F)/np.maximum(
            np.linalg.norm(force,axis=1)[:,None],1e-12)
    def shape_con(x):
        error,_=data(x);limit=PARAMS["max_edge_tolerance_mm"]
        return np.r_[limit-error,limit+error]
    def shape_jac(x):return np.vstack((-J,J))
    start=time.perf_counter()
    sol=minimize(fun,np.zeros(8),jac=jac,method="SLSQP",
        bounds=[(-1.5,1.5)]*8,
        constraints=[dict(type="ineq",fun=force_con,jac=force_jac),
                     dict(type="ineq",fun=shape_con,jac=shape_jac)],
        options=dict(ftol=1e-11,maxiter=200))
    optimize_seconds=time.perf_counter()-start
    d=sol.x*scale
    checks={}
    for name,(mx,my) in (("design_mesh",(nx,ny)),("reference_mesh",(20,16))):
        check_case=build_case(mx,my)
        shell=mitc4_forward(check_case,d)
        corot=mitc4_forward(check_case,d,corotational=True)
        checks[name]={}
        for method,s in (("mitc4_linear",shell),("asd_corotational",corot)):
            q=s["displacement"].ravel()
            error=check_case["geometry_error"]+check_case["O"]@q
            checks[name][method]=dict(status=s["status"],
                **metrics(check_case,q,error),
                max_mount_reaction_N=float(np.max(np.linalg.norm(
                    s["mount_reaction"],axis=1))))
    payload=dict(method="MITC4 linear + finite-difference response basis",
        mesh=[nx,ny],forward_analyses=9,
        adjustment_mm=sol.x.tolist(),optimizer_success=bool(sol.success),
        optimizer_message=str(sol.message),iterations=int(sol.nit),
        objective=float(sol.fun),basis_seconds=basis_seconds,
        optimization_seconds=optimize_seconds,total_seconds=basis_seconds+optimize_seconds,
        checks=checks)
    (OUT/"automotive_mitc4_spring_design.json").write_text(json.dumps(payload,indent=2)+"\n")
    return payload


if __name__=="__main__":
    p=run();print(p["optimizer_success"],p["adjustment_mm"],p["checks"])
