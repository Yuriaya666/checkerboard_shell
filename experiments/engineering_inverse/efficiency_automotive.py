"""Equal-start nonlinear assembly optimization and measured workflow costs.

The two FD variants use the same central 8-variable difference scheme.  A
single cached state supplies objective and both constraint groups.  This
avoids charging either finite-element implementation for redundant SLSQP
callbacks, while all physical parameters remain those in config.json.
"""
from __future__ import annotations

import argparse
import json
import resource
import time
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize,linprog

from experiments.engineering_inverse.automotive_panel import build_case, metrics, PARAMS
from experiments.engineering_inverse.automotive_spring_design import response_sensitivity
from experiments.engineering_inverse.mount_springs import (
    bracket_matrices, checkerboard_forward, mitc4_forward)

OUT=Path(__file__).resolve().parent/"results"
METHODS=("cb_consistent","cb_fd","asd_fd")


def run(method="asd_fd",nx=12,ny=8,start="zero",max_iterations=80,
        difference_step_mm=.05,hessian_assembly_strategy="cached"):
    if method not in METHODS:raise ValueError(method)
    wall_start=time.perf_counter()
    case=build_case(nx,ny)
    model_setup_seconds=time.perf_counter()-wall_start
    Kspring,D,W=bracket_matrices(case) if method=="cb_consistent" else (None,None,None)
    weights=np.tile([8.,1.],nx+1)
    x0=np.zeros(8)
    if start=="cb_feasible":
        source=json.loads((OUT/"automotive_spring_design_12x8.json").read_text())
        x0=np.asarray(source["adjustment_mm"],dtype=float)
    elif start!="zero":raise ValueError(start)
    counts=dict(objective_calls=0,objective_jac_calls=0,
                reaction_constraint_calls=0,reaction_jac_calls=0,
                shape_constraint_calls=0,shape_jac_calls=0,
                equilibrium_solves=0,failed_equilibria=0,retries=0,
                newton_iterations=0,linearizations=0,
                tangent_assemblies=None,numeric_factorizations=None,
                finite_difference_gradients=0)
    times=dict(model_setup_seconds=model_setup_seconds,forward_seconds=0.,
               model_construction_seconds=0.,nonlinear_analysis_seconds=0.,
               extraction_seconds=0.,implicit_gradient_seconds=0.,
               finite_difference_seconds=0.,tangent_assembly_seconds=0.,
               numeric_factorization_seconds=0.)
    if method.startswith("cb"):
        counts["tangent_assemblies"]=0
        counts["numeric_factorizations"]=0
    cache={};history=[]
    suffix="_legacy" if hessian_assembly_strategy=="legacy" else ""
    checkpoint=OUT/f"automotive_efficiency_{method}_{nx}x{ny}_{start}{suffix}.partial.json"

    def calculate(x):
        x=np.asarray(x,dtype=float)
        t=time.perf_counter()
        if method.startswith("cb"):
            # Deterministic state map: no dependence on previously evaluated
            # optimizer trial. This matches the reference solver's cold builds.
            result=checkerboard_forward(case,x*.001,initial=None,
                hessian_assembly_strategy=hessian_assembly_strategy)
            for field,key in (("model_construction_seconds","model_construction_seconds"),
                              ("nonlinear_analysis_seconds","nonlinear_analysis_seconds"),
                              ("extraction_seconds","extraction_seconds")):
                times[key]+=result[field]
            counts["tangent_assemblies"]+=result["tangent_assemblies_total"]
            counts["numeric_factorizations"]+=result["numeric_factorizations_total"]
            times["tangent_assembly_seconds"]+=result["tangent_assembly_seconds"]
            times["numeric_factorization_seconds"]+=result["numeric_factorization_seconds"]
            counts["newton_iterations"]+=int(result["iterations"])
        else:
            result=None
            # Fixed-input diagnostics showed almost every nonzero one-step
            # trial exhausted Newton iterations before converging at two
            # increments. Start with the measured reliable continuation.
            for steps in (2,4,8,16):
                result=mitc4_forward(case,x*.001,corotational=True,
                                     analysis_steps=steps,norm_type=0,
                                     gauge_x_fix=True)
                times["model_construction_seconds"]+=result["model_construction_seconds"]
                times["nonlinear_analysis_seconds"]+=result["analysis_seconds"]
                times["extraction_seconds"]+=result["extraction_seconds"]
                counts["newton_iterations"]+=int(result["newton_iterations_total"])
                counts["equilibrium_solves"]+=1
                if result["status"]==0:break
                counts["failed_equilibria"]+=1;counts["retries"]+=1
        if method.startswith("cb"):
            counts["newton_iterations"]+=int(result["newton_iterations_total"])-int(result["iterations"])
            counts["equilibrium_solves"]+=1+int(result["continuation_steps"])
            counts["retries"]+=int(result["continuation_steps"]>0)
        times["forward_seconds"]+=time.perf_counter()-t
        if (method.startswith("cb") and result["status"]!="CONVERGED") or (
                method=="asd_fd" and result["status"]!=0):
            counts["failed_equilibria"]+=1
            raise RuntimeError(f"{method} equilibrium failed at {x.tolist()}: {result['status']}")
        q=result["displacement"].ravel()
        error=(case["geometry_error"]+case["O"]@q)*1000.
        F=result["mount_reaction"]
        objective=float((weights*error)@error/len(error)+.0005*x@x)
        return dict(solution=result,error_mm=error,reaction=F,objective=objective,
                    constraint=np.r_[PARAMS["reaction_limit_N"]-np.linalg.norm(F,axis=1),
                        PARAMS["max_edge_tolerance_mm"]-error,
                        PARAMS["max_edge_tolerance_mm"]+error])

    def state(x,record=True):
        key=tuple(np.asarray(x,dtype=float))
        if key!=cache.get("key"):
            data=calculate(x)
            cache.clear();cache.update(key=key,data=data)
            if record:
                history.append(dict(adjustment_mm=list(key),objective=data["objective"],
                    max_edge_error_mm=float(np.max(abs(data["error_mm"]))),
                    max_mount_reaction_N=float(np.max(np.linalg.norm(data["reaction"],axis=1))),
                    elapsed_seconds=time.perf_counter()-wall_start))
                if len(history)%5==0:
                    checkpoint.write_text(json.dumps(dict(method=method,start=start,
                        elapsed_seconds=time.perf_counter()-wall_start,
                        counts=counts,latest=history[-1]),indent=2)+"\n")
                print("evaluation",method,len(history),round(history[-1]["objective"],8),
                      round(history[-1]["max_edge_error_mm"],6),
                      round(history[-1]["elapsed_seconds"],2),flush=True)
        return cache["data"]

    def derivatives(x):
        key=tuple(np.asarray(x,dtype=float))
        state(x)
        if cache.get("derivative_key")==key:return cache["derivative"]
        if method=="cb_consistent":
            t=time.perf_counter()
            J,dF=response_sensitivity(case,cache["data"]["solution"],Kspring,D,W,
                hessian_assembly_strategy=hessian_assembly_strategy)
            times["implicit_gradient_seconds"]+=time.perf_counter()-t
            counts["linearizations"]+=1
            J*=1.0  # mm output per mm input; displacement derivative is m/m.
            dF*=.001
            objective_jac=2*J.T@(weights*cache["data"]["error_mm"])/len(weights)+.001*np.asarray(x)
            F=cache["data"]["reaction"]
            reaction_jac=-np.einsum("pc,pcj->pj",F,dF)/np.maximum(
                np.linalg.norm(F,axis=1)[:,None],1e-12)
            shape_jac=np.vstack((-J,J))
            cjac=np.vstack((reaction_jac,shape_jac))
        else:
            t=time.perf_counter()
            objective_jac=np.empty(8);cjac=np.empty((4+2*len(cache["data"]["error_mm"]),8))
            for i in range(8):
                xp=np.asarray(x,dtype=float).copy();xp[i]+=difference_step_mm
                xm=np.asarray(x,dtype=float).copy();xm[i]-=difference_step_mm
                plus=calculate(xp);minus=calculate(xm)
                objective_jac[i]=(plus["objective"]-minus["objective"])/(2*difference_step_mm)
                cjac[:,i]=(plus["constraint"]-minus["constraint"])/(2*difference_step_mm)
            times["finite_difference_seconds"]+=time.perf_counter()-t
            counts["finite_difference_gradients"]+=1
        cache.update(derivative_key=key,derivative=(objective_jac,cjac))
        return cache["derivative"]

    def fun(x):counts["objective_calls"]+=1;return state(x)["objective"]
    def jac(x):counts["objective_jac_calls"]+=1;return derivatives(x)[0]
    def rf(x):counts["reaction_constraint_calls"]+=1;return state(x)["constraint"][:4]
    def rj(x):counts["reaction_jac_calls"]+=1;return derivatives(x)[1][:4]
    def sf(x):counts["shape_constraint_calls"]+=1;return state(x)["constraint"][4:]
    def sj(x):counts["shape_jac_calls"]+=1;return derivatives(x)[1][4:]
    constraints=[dict(type="ineq",fun=rf,jac=rj),dict(type="ineq",fun=sf,jac=sj)]
    optimum=minimize(fun,x0,jac=jac,method="SLSQP",bounds=[(-1.5,1.5)]*8,
        constraints=constraints,options=dict(ftol=1e-10,maxiter=max_iterations))
    diagnostic=None
    if not optimum.success and optimum.nit<=1:
        initial_state=state(x0)
        _,initial_jac=derivatives(x0)
        linearized=linprog(np.zeros(8),A_ub=-initial_jac,
            b_ub=initial_state["constraint"],
            bounds=[(-1.5,1.5)]*8,method="highs")
        diagnostic=dict(initial_min_constraint=float(np.min(initial_state["constraint"])),
            initial_jac_max=float(np.max(abs(initial_jac))),
            initial_jac_rank=int(np.linalg.matrix_rank(initial_jac)),
            initial_linearized_feasible=bool(linearized.success),
            initial_linearized_message=linearized.message)
    final=state(optimum.x)
    r=final["solution"]
    quality=dict(**metrics(case,r["displacement"].ravel(),final["error_mm"]*.001),
        objective=final["objective"],max_mount_reaction_N=float(np.max(np.linalg.norm(
            r["mount_reaction"],axis=1))),mount_reaction_N=r["mount_reaction"].tolist(),
        engineering_feasible=bool(np.min(final["constraint"])>=-1e-5 and
            np.max(np.abs(optimum.x))<=1.5+1e-8),
        force_balance_N=r["force_balance"].tolist(),
        residual_inf=(float(r["residual_inf"]) if method.startswith("cb") else None))
    payload=dict(method=method,mesh=[nx,ny],start=start,
        hessian_assembly_strategy=hessian_assembly_strategy,
        state_initialization="reference-state cold solve for every new design input",
        common_x_locator=True,
        start_adjustment_mm=x0.tolist(),difference_step_mm=difference_step_mm,
        adjustment_mm=optimum.x.tolist(),optimizer_success=bool(optimum.success),
        optimizer_message=str(optimum.message),iterations=int(optimum.nit),
        failure_diagnostic=diagnostic,
        total_seconds=time.perf_counter()-wall_start,
        peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
        state_dofs=(3 if method.startswith("cb") else 6)*case["mesh"].num_vertices,
        free_shell_dofs=(3 if method.startswith("cb") else 6)*case["mesh"].num_vertices-
            3*len(case["datum"])-(1 if method.startswith("cb") else 0),
        counts=counts,times=times,final=quality,history=history)
    path=OUT/f"automotive_efficiency_{method}_{nx}x{ny}_{start}{suffix}.json"
    path.write_text(json.dumps(payload,indent=2)+"\n")
    if checkpoint.exists():checkpoint.unlink()
    print(method,start,optimum.success,optimum.nit,payload["total_seconds"],
          quality["max_abs_error_mm"],quality["max_mount_reaction_N"],flush=True)
    return payload


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("method",choices=METHODS)
    parser.add_argument("--start",choices=("zero","cb_feasible"),default="zero")
    parser.add_argument("--nx",type=int,default=12)
    parser.add_argument("--ny",type=int,default=8)
    parser.add_argument("--max-iterations",type=int,default=80)
    parser.add_argument("--difference-step-mm",type=float,default=.05)
    parser.add_argument("--hessian-strategy",choices=("cached","legacy"),
                        default="cached")
    args=parser.parse_args()
    run(args.method,args.nx,args.ny,args.start,args.max_iterations,
        args.difference_step_mm,args.hessian_strategy)
