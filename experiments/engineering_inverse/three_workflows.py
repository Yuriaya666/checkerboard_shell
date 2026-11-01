"""Three equal-layout linear aircraft design workflows, twist case.

Each row uses the same geometry, pressure, strokes, target, and force limits.
The MITC4 state solve is OpenSees ShellMITC4 with frame beam lines; the
Checkerboard state solve is its production linearized tangent plus position
beam lines. Full nonlinear acceptance is reported separately.
"""

import json
import time

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from scipy.optimize import minimize

from experiments.engineering_inverse.aircraft_skin import (
    OUT,PARAMS,build_case,layout_operator)
from experiments.engineering_inverse.common import (
    checkerboard_displacement_operator,mitc4_prescribed_response)


def run():
    case=build_case(12,24)
    source=json.loads((OUT/"aircraft_nonlinear_final_12x24.json").read_text())
    selected=source["stations"]
    base=checkerboard_displacement_operator(case["mesh"],case["material"],
        case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"])
    op=layout_operator(case,selected,{"K":base["K"],"residual":base["residual"]})
    groups=op["groups"]
    target=case["targets"][1]
    n=len(groups);scale=.002
    K=op["K"];free=op["free"];fixed=op["fixed"]
    lu=spla.splu(sp.csc_matrix(K[np.ix_(free,free)]))
    B=np.zeros((K.shape[0],n))
    for k,g in enumerate(groups):B[g,k]=1
    counter={"analyses":0}
    def cb_forward(d):
        counter["analyses"]+=1
        q=B@d
        q[free]=lu.solve(-op["residual"][free]-K[np.ix_(free,fixed)]@q[fixed])
        return case["O"]@q
    def cb_exact(d):return op["response0"]+op["response_map"]@d
    G=op["actuator_reaction_map"]
    b=op["actuator_reaction0"]
    def constraints(d):return np.r_[80-b-G@d,80+b+G@d]
    def constraints_jac(d):return np.vstack((-G,G))
    entries=[]
    for method in ("Checkerboard + exact derivative","Checkerboard + central FD"):
        counter["analyses"]=0
        exact=method.endswith("exact derivative")
        response=cb_exact if exact else cb_forward
        evaluations={"objective":0,"gradient":0}
        def objective(d):
            evaluations["objective"]+=1
            r=(response(d)-target)/scale
            return float(r@r/len(r)+1e-4*d@d/scale**2)
        def jac(d):
            evaluations["gradient"]+=1
            if exact:
                r=response(d)-target
                return 2*op["response_map"].T@r/(len(r)*scale**2)+2e-4*d/scale**2
            eps=1e-6
            grad=np.empty(n)
            for i in range(n):
                step=np.zeros(n);step[i]=eps
                grad[i]=(objective(d+step)-objective(d-step))/(2*eps)
            return grad
        start=time.perf_counter()
        result=minimize(objective,np.zeros(n),jac=jac,method="SLSQP",
            bounds=[(-PARAMS["stroke_limit_m"],PARAMS["stroke_limit_m"])]*n,
            constraints=[dict(type="ineq",fun=constraints,jac=constraints_jac)],
            options=dict(ftol=1e-11,maxiter=200))
        elapsed=time.perf_counter()-start
        error=response(result.x)-target
        entries.append(dict(method=method,success=bool(result.success),
            iterations=int(result.nit),objective=float(result.fun),
            rms_error_mm=float(np.sqrt(np.mean(error**2))*1e3),
            max_error_mm=float(np.max(np.abs(error))*1e3),
            strokes_m=result.x.tolist(),reaction_N=(b+G@result.x).tolist(),
            forward_solves=counter["analyses"],
            objective_evaluations=evaluations["objective"],
            gradient_evaluations=evaluations["gradient"],
            design_seconds=elapsed,
            common_assembly_seconds=base["assembly_seconds"]))
    fixed_list=list(fixed)
    lookup={int(d):i for i,d in enumerate(fixed_list)}
    baseline_reaction=[];reaction_columns=[]
    start=time.perf_counter()
    _,_,r0=mitc4_prescribed_response(case["mesh"],case["material"],
        fixed_list,groups,np.zeros(n),case["O"],case["load"],beam_lines=case["beam_lines"])
    baseline_reaction=np.array([r0[[lookup[d] for d in g]].sum() for g in groups])
    for i in range(n):
        probe=np.zeros(n);probe[i]=.001
        _,_,ri=mitc4_prescribed_response(case["mesh"],case["material"],
            fixed_list,groups,probe,case["O"],case["load"],beam_lines=case["beam_lines"])
        reaction_columns.append([(ri[[lookup[d] for d in g]].sum()-baseline_reaction[k])/.001
                                 for k,g in enumerate(groups)])
    reaction_map=np.column_stack(reaction_columns)
    basis_seconds=time.perf_counter()-start
    mt_counter={"analyses":n+1,"objective":0,"gradient":0}
    cache={}
    def mt_forward(d):
        key=tuple(float(x) for x in d)
        if key not in cache:
            response,_,_=mitc4_prescribed_response(case["mesh"],case["material"],
                fixed_list,groups,d,case["O"],case["load"],beam_lines=case["beam_lines"])
            cache[key]=response;mt_counter["analyses"]+=1
        return cache[key]
    def mt_objective(d):
        mt_counter["objective"]+=1
        r=(mt_forward(d)-target)/scale
        return float(r@r/len(r)+1e-4*d@d/scale**2)
    def mt_jac(d):
        mt_counter["gradient"]+=1
        eps=1e-6
        grad=np.empty(n)
        for i in range(n):
            step=np.zeros(n);step[i]=eps
            grad[i]=(mt_objective(d+step)-mt_objective(d-step))/(2*eps)
        return grad
    def mt_constraints(d):return np.r_[80-baseline_reaction-reaction_map@d,
                                        80+baseline_reaction+reaction_map@d]
    def mt_constraints_jac(d):return np.vstack((-reaction_map,reaction_map))
    start=time.perf_counter()
    result=minimize(mt_objective,np.zeros(n),jac=mt_jac,method="SLSQP",
        bounds=[(-PARAMS["stroke_limit_m"],PARAMS["stroke_limit_m"])]*n,
        constraints=[dict(type="ineq",fun=mt_constraints,jac=mt_constraints_jac)],
        options=dict(ftol=1e-11,maxiter=200))
    elapsed=time.perf_counter()-start
    error=mt_forward(result.x)-target
    entries.append(dict(method="MITC4 + central FD",success=bool(result.success),
        iterations=int(result.nit),objective=float(result.fun),
        rms_error_mm=float(np.sqrt(np.mean(error**2))*1e3),
        max_error_mm=float(np.max(np.abs(error))*1e3),
        strokes_m=result.x.tolist(),
        reaction_N=(baseline_reaction+reaction_map@result.x).tolist(),
        forward_solves=mt_counter["analyses"],
        objective_evaluations=mt_counter["objective"],
        gradient_evaluations=mt_counter["gradient"],
        design_seconds=elapsed,basis_seconds=basis_seconds))
    data=dict(case="aircraft twist, fixed verified three-station layout",
              stations=selected,methods=entries,
              notice="linear screening comparison only; MITC4 uses full rebuild per forward, CB reuses factorization")
    (OUT/"three_workflows.json").write_text(json.dumps(data,indent=2)+"\n")
    return data


if __name__=="__main__":
    for row in run()["methods"]:print(row["method"],row["rms_error_mm"],row["design_seconds"],row["forward_solves"])
