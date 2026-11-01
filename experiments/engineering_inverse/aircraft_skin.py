"""Distributed-rib upper-skin demonstrator with stroke layouts.

This is a scaled representative test article, not the published FishBAC CAD.
It contains a bonded thin metal skin, two spanwise beam stringers, five compliant
chordwise rib lines, a spar clamp, and selectable stroke-driven attachments.
"""

from __future__ import annotations

import itertools
import json
import time
from pathlib import Path

import numpy as np
import scipy.optimize as opt
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch

from checkerboard_shell import (BoundaryConditions, CheckerboardMesh,
                                CheckerboardShell, IsotropicMaterial, NodalLoad,
                                NewtonSettings)
from experiments.engineering_inverse.common import (
    checkerboard_displacement_operator, mitc4_prescribed_response,
    positional_beam_stiffness)
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent

HERE = Path(__file__).resolve().parent
OUT = HERE/"results"
OUT.mkdir(exist_ok=True)
PARAMS=json.loads((HERE/"config.json").read_text())["aircraft"]


def naca2412_upper(xc):
    xc = np.asarray(xc)
    m,p,t = .02,.4,.12
    camber = np.where(xc<=p, m/p**2*(2*p*xc-xc**2),
        m/(1-p)**2*((1-2*p)+2*p*xc-xc**2))
    thickness = 5*t*(.2969*np.sqrt(xc)-.1260*xc-.3516*xc**2+
                      .2843*xc**3-.1015*xc**4)
    return camber+thickness


def beam_line_strain_estimate(case, displacement):
    """Maximum linearized axial and outer-fiber bending strain of beam lines."""
    xyz=case["mesh"].reference_vertices.numpy().reshape(-1,3)
    q=np.asarray(displacement).reshape(-1,3)
    maxima=[]
    for line,E,area,inertia in case["beam_lines"]:
        axial=[]; bending=[]
        depth=np.sqrt(12*inertia/area)
        for a,b in zip(line[:-1],line[1:]):
            delta=xyz[b]-xyz[a];length=np.linalg.norm(delta)
            axial.append(abs(np.dot(q[b]-q[a],delta/length))/length)
        for a,b,c in zip(line[:-2],line[1:-1],line[2:]):
            length=.5*(np.linalg.norm(xyz[b]-xyz[a])+
                       np.linalg.norm(xyz[c]-xyz[b]))
            curvature=np.linalg.norm(q[a]-2*q[b]+q[c])/length**2
            bending.append(curvature*depth/2)
        maxima.append(dict(max_axial_strain=float(max(axial)),
                           max_outer_bending_strain=float(max(bending)),
                           E_Pa=E,area_m2=area,inertia_m4=inertia))
    return maxima


def build_case(nx=12, ny=24):
    chord, span = PARAMS["chord_m"],PARAMS["span_m"]
    aft = PARAMS["aft_m"]
    base = CheckerboardMesh.rectangle(aft, span, nx, ny)
    reference = base.reference_vertices.clone()
    x = reference[...,0].numpy()+chord-aft
    reference[...,2] = torch.as_tensor(chord*naca2412_upper(x/chord))
    mesh = CheckerboardMesh(reference,base.hx,base.hy,base.quads,base.boundaries)
    node = lambda i,j: j*(nx+1)+i
    station_fractions=tuple(PARAMS["candidate_span_fractions"])
    stations = [int(round(f*ny)) for f in station_fractions]
    lines = []
    for i in (nx//3,2*nx//3):
        line = [node(i,j) for j in range(ny+1)]
        E=70e9; width=.006; depth=.001
        lines.append((line,E,width*depth,width*depth**3/12))
    for j in stations:
        line = [node(i,j) for i in range(nx+1)]
        E=70e9; width=.008; depth=.001
        lines.append((line,E,width*depth,width*depth**3/12))
    pressure = np.empty((ny,nx))
    for j in range(ny):
        for i in range(nx):
            xi=(i+.5)/nx; eta=(j+.5)/ny
            pressure[j,i]=PARAMS["pressure_peak_Pa"]*(1-.4*xi)*(1+.15*np.cos(2*np.pi*eta))
    load = NodalLoad.dead_pressure(mesh,pressure)
    observations = [node(i,j) for j in range(ny+1) for i in range(2,nx+1)]
    O = np.zeros((len(observations),3*mesh.num_vertices))
    for row,n in enumerate(observations): O[row,3*n+2]=1
    targets = np.empty((2,len(observations)))
    for row,n in enumerate(observations):
        i=n%(nx+1); j=n//(nx+1); xi=i/nx; eta=j/ny
        shape=1.5*xi**2-.5*xi**3
        targets[0,row]=PARAMS["uniform_trailing_edge_increment_m"]*shape
        targets[1,row]=PARAMS["uniform_trailing_edge_increment_m"]*shape*(.4+1.2*eta)
    spar = [3*node(i,j)+c for j in range(ny+1)
            for i in (0,1) for c in range(3)]
    patches = []
    for station_index,j in enumerate(stations):
        y_center=station_fractions[station_index]*span
        patches.append([node(i,jj) for jj in range(ny+1)
                        for i in range(nx+1)
                        if .83 <= i/nx <= 1.0
                        and abs(jj*span/ny-y_center) <= .05])
    return dict(mesh=mesh,material=IsotropicMaterial(PARAMS["skin_E_Pa"],
                                  PARAMS["skin_poisson"],PARAMS["skin_thickness_m"]),
                load=load,beam_lines=lines,station_rows=stations,
                patches=patches,spar=spar,O=O,targets=targets,
                station_fractions=station_fractions,
                observations=observations,pressure=pressure,
                chord=chord,span=span,aft=aft)


def layout_operator(case, selected, prepared):
    if "compliance" in prepared:
        return _layout_from_compliance(case, selected, prepared)
    patches=[case["patches"][s] for s in selected]
    groups=[[3*n+2 for n in patch] for patch in patches]
    fixed=case["spar"]+[d for group in groups for d in group]
    op = checkerboard_displacement_operator(case["mesh"],case["material"],
        fixed,groups,case["O"],case["load"],prepared=prepared)
    # Net stroke reaction for each physical attachment, not per FE node.
    local={d:k for k,d in enumerate(op["fixed"])}
    reaction0=[]; reaction_map=[]
    for group in groups:
        indices=[local[d] for d in group]
        reaction0.append(op["reaction0"][indices].sum())
        reaction_map.append(op["reaction_map"][indices].sum(axis=0))
    op["actuator_reaction0"]=np.asarray(reaction0)
    op["actuator_reaction_map"]=np.asarray(reaction_map).reshape(len(groups),len(groups))
    op["groups"]=groups
    return op


def prepare_layout_compliance(case, base):
    """Reuse one spar-clamped factorization across every actuator subset.

    For candidate DOFs C, H=K_UU^{-1}E_C and Cmat=E_C^T H. A selected
    layout enforces q_C=B d via Cmat_SS lambda=B d-q_C^0. The resulting
    state and reaction are exactly those of the eliminated linear system.
    """
    candidate=sorted({3*n+2 for patch in case["patches"] for n in patch})
    free=base["free"]
    index=np.searchsorted(free,candidate)
    if not np.array_equal(free[index],candidate):
        raise ValueError("candidate actuator DOF is in the spar constraint")
    rhs=np.zeros((len(free),len(candidate)))
    rhs[index,np.arange(len(candidate))]=1.0
    influence=base["factorized"](rhs)
    compliance=influence[index,:]
    return dict(base=base,candidate=np.asarray(candidate),
                candidate_lookup={d:i for i,d in enumerate(candidate)},
                influence=influence,compliance=compliance,
                spar_influence=(base["cross_action"](influence)
                    if "cross_action" in base else
                    base["K"][np.ix_(base["fixed"],free)]@influence))


def _layout_from_compliance(case, selected, prepared):
    base=prepared["base"]
    groups=[[3*n+2 for n in case["patches"][s]] for s in selected]
    selected_dofs=sorted(set(d for group in groups for d in group))
    indices=np.array([prepared["candidate_lookup"][d] for d in selected_dofs])
    C=prepared["compliance"][np.ix_(indices,indices)]
    B=np.zeros((len(indices),len(groups)))
    lookup={d:i for i,d in enumerate(selected_dofs)}
    for col,group in enumerate(groups):
        B[[lookup[d] for d in group],col]=1.0
    base_at_controls=base["q0"][selected_dofs]
    lambda0=np.linalg.solve(C,-base_at_controls)
    lambda_map=np.linalg.solve(C,B)
    H=prepared["influence"][:,indices]
    q0=base["q0"].copy()
    q0[base["free"]]+=H@lambda0
    Q=np.zeros((len(q0),len(groups)))
    Q[base["free"]]=H@lambda_map
    fixed=np.array(sorted(set(base["fixed"]).union(selected_dofs)))
    free=np.setdiff1d(base["free"],selected_dofs)
    # Base spar reaction plus the change induced by candidate multipliers.
    spar0=base["reaction0"]+prepared["spar_influence"][:,indices]@lambda0
    spar_map=prepared["spar_influence"][:,indices]@lambda_map
    reaction0=np.zeros(len(fixed));reaction_map=np.zeros((len(fixed),len(groups)))
    fixed_lookup={d:i for i,d in enumerate(fixed)}
    reaction0[[fixed_lookup[d] for d in base["fixed"]]]=spar0
    reaction_map[[fixed_lookup[d] for d in base["fixed"]]]=spar_map
    reaction0[[fixed_lookup[d] for d in selected_dofs]]=lambda0
    reaction_map[[fixed_lookup[d] for d in selected_dofs]]=lambda_map
    net0=np.array([lambda0[[lookup[d] for d in group]].sum() for group in groups])
    netmap=np.array([lambda_map[[lookup[d] for d in group]].sum(axis=0)
                     for group in groups])
    return dict(q0=q0,Q=Q,response0=case["O"]@q0,
                response_map=case["O"]@Q,reaction0=reaction0,
                reaction_map=reaction_map,actuator_reaction0=net0,
                actuator_reaction_map=netmap,groups=groups,
                fixed=fixed,free=free,K=base.get("K"),residual=base["residual"],
                free_dofs=len(free),total_dofs=len(q0),assembly_seconds=0.0)


def optimize_controls(op,target,stroke_limit=None,force_limit=None):
    if stroke_limit is None: stroke_limit=PARAMS["stroke_limit_m"]
    if force_limit is None: force_limit=PARAMS["reaction_limit_N"]
    S=op["response_map"]
    r0=op["response0"]-target
    n=S.shape[1]
    scale=.003
    regularizer=.0001
    def fun(d):
        r=(r0+S@d)/scale
        return float(r@r/len(r)+regularizer*d@d/scale**2)
    def jac(d):
        r=(r0+S@d)/scale
        return 2*S.T@r/(len(r)*scale)+2*regularizer*d/scale**2
    A=op["actuator_reaction_map"]
    b=op["actuator_reaction0"]
    constraints=[]
    if n:
        constraints=[dict(type="ineq",fun=lambda d: force_limit-b-A@d,
                          jac=lambda d:-A),
                     dict(type="ineq",fun=lambda d: force_limit+b+A@d,
                          jac=lambda d:A)]
    start=time.perf_counter()
    sol=opt.minimize(fun,np.zeros(n),jac=jac,method="SLSQP",
                     bounds=[(-stroke_limit,stroke_limit)]*n,
                     constraints=constraints,
                     options=dict(ftol=1e-11,maxiter=300))
    elapsed=time.perf_counter()-start
    err=r0+S@sol.x
    reaction=b+A@sol.x
    return dict(strokes_m=sol.x.tolist(),reactions_N=reaction.tolist(),
                rms_error_mm=float(np.sqrt(np.mean(err**2))*1e3),
                max_error_mm=float(np.max(np.abs(err))*1e3),
                objective=float(sol.fun),iterations=int(sol.nit),
                success=bool(sol.success),time_seconds=elapsed,
                force_ok=bool(np.max(np.abs(reaction))<=force_limit+1e-5),
                stroke_ok=bool(np.max(np.abs(sol.x))<=stroke_limit+1e-8))


def nonlinear_optimize_controls(case, selected, initial_strokes, target,
                                max_iterations=8, initial_displacement=None):
    """Sequential exact-tangent design correction for the coupled shell/lines."""
    groups=[[3*n+2 for n in case["patches"][s]] for s in selected]
    fixed=np.array(sorted(set(case["spar"]+[d for g in groups for d in g])))
    free=np.setdiff1d(np.arange(3*case["mesh"].num_vertices),fixed)
    mask=torch.zeros(case["mesh"].shape,dtype=torch.bool)
    mask.reshape(-1)[fixed]=True
    Kbeam=sum((positional_beam_stiffness(case["mesh"].reference_vertices.numpy(),
                                         *line) for line in case["beam_lines"]),
              np.zeros((3*case["mesh"].num_vertices,)*2))
    prescribed=torch.zeros(case["mesh"].shape,dtype=torch.float64)
    model=CheckerboardShell(case["mesh"],case["material"],
        BoundaryConditions(mask,prescribed),case["load"],
        NewtonSettings(residual_tolerance=1e-5,maximum_newton_iterations=30),
        sp.csr_matrix(Kbeam))
    B=np.zeros((3*case["mesh"].num_vertices,len(groups)))
    for col,g in enumerate(groups): B[g,col]=1
    state=(None if initial_displacement is None else
           torch.as_tensor(initial_displacement,dtype=torch.float64).reshape(case["mesh"].shape))
    history=[]
    work={"equilibrium_solves":0,"newton_iterations":0}
    def evaluate(d,with_tangent=False):
        nonlocal state
        values=torch.zeros(case["mesh"].shape,dtype=torch.float64)
        values.reshape(-1)[:]=torch.as_tensor(B@d,dtype=torch.float64)
        model.boundary_conditions=BoundaryConditions(mask,values)
        result=model.solve(initial_displacement=state)
        work["equilibrium_solves"]+=1
        work["newton_iterations"]+=result.iterations
        if result.status!="CONVERGED":
            # A failed direct Newton step is a solver outcome, not proof that
            # this actuator layout is mechanically infeasible. Ramp both
            # pressure and stroke from the stress-free reference instead.
            fraction=0.0;step=.25
            continuation_state=torch.zeros(case["mesh"].shape,dtype=torch.float64)
            while fraction<1.0-1e-12:
                next_fraction=min(1.0,fraction+step)
                model.boundary_conditions=BoundaryConditions(
                    mask,values*next_fraction)
                trial=model.solve(initial_displacement=continuation_state,
                                  load_factor=next_fraction)
                work["equilibrium_solves"]+=1
                work["newton_iterations"]+=trial.iterations
                if trial.status=="CONVERGED":
                    continuation_state=(trial.deformed_vertices-
                        case["mesh"].reference_vertices).detach()
                    fraction=next_fraction
                    result=trial
                    step=min(step*1.5,.5)
                else:
                    step*=.5
                    if step<1/128:
                        model.boundary_conditions=BoundaryConditions(mask,values)
                        raise RuntimeError(
                            f"nonlinear equilibrium: direct={result.status}, "
                            f"continuation={trial.status} at {fraction:.4f}")
            model.boundary_conditions=BoundaryConditions(mask,values)
        state=(result.deformed_vertices-case["mesh"].reference_vertices).detach()
        q=state.numpy().reshape(-1)
        r=case["O"]@q-target
        force=model.full_residual_and_reactions(result.coordinates)[1].numpy().reshape(-1)
        reactions=np.array([force[g].sum() for g in groups])
        out=dict(rms_error_mm=float(np.sqrt(np.mean(r**2))*1e3),
                 max_error_mm=float(np.max(np.abs(r))*1e3),
                 max_error_node=int(case["observations"][int(np.argmax(np.abs(r)))]),
                 max_error_signed_mm=float(r[np.argmax(np.abs(r))]*1e3),
                 reactions_N=reactions.tolist(),
                 newton_iterations=result.iterations,
                 residual_inf=result.residual_inf,
                 strain_energy_J=result.energies["total"])
        if not with_tangent: return out,q
        tangent,_,_=assemble_paper_low_rank_tangent(
            state,case["mesh"].reference_vertices,case["mesh"].hx,
            case["mesh"].hy,model.targets,case["material"],variant="CB-Full")
        K=sp.csr_matrix(tangent.local+sp.csr_matrix(Kbeam))
        K=K+tangent.coefficient*sp.csr_matrix(np.outer(tangent.vector,tangent.vector))
        lu=spla.splu(K[free][:,free].tocsc())
        Q=B.copy()
        Q[free]=lu.solve(-K[free,:]@B)
        J=case["O"]@Q
        return out,(r,J)
    d=np.asarray(initial_strokes,dtype=float).copy()
    start=time.perf_counter()
    for iteration in range(max_iterations):
        out,derivative=evaluate(d,True)
        r,J=derivative
        history.append(dict(iteration=iteration,strokes_m=d.tolist(),**out))
        regularizer=1e-4
        h=J.T@J+regularizer*np.eye(len(d))
        gradient=J.T@r+regularizer*d
        step=np.linalg.solve(h,-gradient)
        if np.linalg.norm(step)<1e-7: break
        d=np.clip(d+step,-PARAMS["stroke_limit_m"],PARAMS["stroke_limit_m"])
    final,final_q=evaluate(d,False)
    final.update(strokes_m=d.tolist(),history=history,
                 displacement_m=final_q.tolist(),
                 beam_line_strain_estimate=beam_line_strain_estimate(case,final_q),
                 equilibrium_solves=work["equilibrium_solves"],
                 newton_total_iterations=work["newton_iterations"],
                 time_seconds=time.perf_counter()-start,
                 success=bool(final["rms_error_mm"]<=PARAMS["rms_tolerance_mm"] and
                              final["max_error_mm"]<=PARAMS["max_tolerance_mm"] and
                              np.max(np.abs(final["reactions_N"]))<=PARAMS["reaction_limit_N"]))
    return final


def run(nx=12,ny=24,with_mitc4=True):
    case=build_case(nx,ny)
    assembled=checkerboard_displacement_operator(case["mesh"],case["material"],
        case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"])
    preparation_start=time.perf_counter()
    prepared=prepare_layout_compliance(case,assembled)
    preparation_seconds=time.perf_counter()-preparation_start
    search=[]
    start=time.perf_counter()
    for count in range(1,6):
        for subset in itertools.combinations(range(5),count):
            op=layout_operator(case,subset,prepared)
            modes=[optimize_controls(op,target,
                PARAMS["stroke_limit_m"],PARAMS["reaction_limit_N"])
                for target in case["targets"]]
            feasible=all(m["success"] and m["force_ok"] and m["stroke_ok"]
                and m["rms_error_mm"]<=PARAMS["rms_tolerance_mm"]
                and m["max_error_mm"]<=PARAMS["max_tolerance_mm"]
                for m in modes)
            search.append(dict(stations=list(subset),count=count,
                               cases=modes,feasible=feasible,
                               worst_rms_mm=max(m["rms_error_mm"] for m in modes),
                               worst_max_mm=max(m["max_error_mm"] for m in modes)))
    search_time=time.perf_counter()-start
    feasible=[r for r in search if r["feasible"]]
    winner=min(feasible,key=lambda r:(r["count"],r["worst_rms_mm"])) if feasible else None
    payload=dict(case="scaled_distributed_rib_metal_upper_skin",
        mesh=[nx,ny],geometry=dict(chord_m=case["chord"],span_m=case["span"],
                                   aft_fraction=.30,airfoil="NACA 2412 upper"),
        material=dict(skin_E_Pa=70e9,skin_nu=.33,skin_t_m=.0005,
                      stringer_E_Pa=70e9,rib_E_Pa=70e9,
                      material_status="research assumptions; thin metallic compliant skin, not published FishBAC properties"),
        load=dict(pressure_range_Pa=[float(case["pressure"].min()),
                                      float(case["pressure"].max())],
                  type="reference-normal dead pressure"),
        actuation=dict(control="vertical stroke at 2x2 node attachment",
                       stroke_limit_mm=PARAMS["stroke_limit_m"]*1000,
                       net_reaction_limit_N=PARAMS["reaction_limit_N"],
                       candidates_span_fraction=list(case["station_fractions"])),
        target=dict(case_1="2 mm uniform TE increment with clamped-root cubic chord law",
                    case_2="span-gradient TE increment 0.8 to 3.2 mm",
                    rms_tolerance_mm=PARAMS["rms_tolerance_mm"],
                    max_tolerance_mm=PARAMS["max_tolerance_mm"]),
        checkerboard=dict(free_dofs_base=assembled["free_dofs"],
                          total_dofs=assembled["total_dofs"],
                          shell_assembly_seconds=assembled["assembly_seconds"],
                          compliance_preparation_seconds=preparation_seconds,
                          layout_search_seconds=search_time,
                          total_seconds=assembled["assembly_seconds"]+preparation_seconds+search_time,
                          layouts_evaluated=len(search),search=search,winner=winner))
    if with_mitc4:
        mt_search=[]; mt_analyses=0
        start=time.perf_counter()
        for cb_row in search:
            selected=cb_row["stations"]
            groups=[[3*n+2 for n in case["patches"][s]] for s in selected]
            fixed=sorted(set(case["spar"]+[d for g in groups for d in g]))
            zero=np.zeros(len(groups))
            mt0,_,mr0=mitc4_prescribed_response(case["mesh"],case["material"],
                fixed,groups,zero,case["O"],case["load"],beam_lines=case["beam_lines"])
            mt_analyses+=1
            probes=[]; rprobes=[]
            for i in range(len(groups)):
                values=np.zeros(len(groups)); values[i]=.001
                response,_,reaction=mitc4_prescribed_response(
                    case["mesh"],case["material"],fixed,groups,values,case["O"],
                    case["load"],beam_lines=case["beam_lines"])
                mt_analyses+=1
                probes.append((response-mt0)/.001)
                rprobes.append((reaction-mr0)/.001)
            M=np.column_stack(probes); RM=np.column_stack(rprobes)
            local={d:k for k,d in enumerate(fixed)}
            net0=np.array([mr0[[local[d] for d in g]].sum() for g in groups])
            netM=np.array([RM[[local[d] for d in g]].sum(axis=0) for g in groups])
            mt_op=dict(response0=mt0,response_map=M,
                       actuator_reaction0=net0,actuator_reaction_map=netM)
            modes=[optimize_controls(mt_op,target,
                PARAMS["stroke_limit_m"],PARAMS["reaction_limit_N"])
                for target in case["targets"]]
            feasible=all(m["success"] and m["force_ok"] and m["stroke_ok"]
                and m["rms_error_mm"]<=PARAMS["rms_tolerance_mm"]
                and m["max_error_mm"]<=PARAMS["max_tolerance_mm"]
                for m in modes)
            mt_search.append(dict(stations=selected,count=len(selected),
                cases=modes,feasible=feasible,
                worst_rms_mm=max(m["rms_error_mm"] for m in modes),
                worst_max_mm=max(m["max_error_mm"] for m in modes)))
        mt_time=time.perf_counter()-start
        mt_feasible=[r for r in mt_search if r["feasible"]]
        mt_winner=min(mt_feasible,key=lambda r:(r["count"],r["worst_rms_mm"])) if mt_feasible else None
        mt_same=next((r for r in mt_search if winner is not None and
                      r["stations"]==winner["stations"]),None)
        payload["mitc4"]=dict(method="OpenSees ShellMITC4 + elasticBeamColumn",
            free_dofs_base=6*case["mesh"].num_vertices-len(case["spar"]),
            search_seconds=mt_time,forward_analyses=mt_analyses,
            layouts_evaluated=len(mt_search),search=mt_search,
            winner=mt_winner,same_checkerboard_layout=mt_same)
    path=OUT/f"aircraft_{nx}x{ny}.json"
    path.write_text(json.dumps(payload,indent=2)+"\n")
    return payload


if __name__=="__main__":
    data=run()
    print(json.dumps({"winner":data["checkerboard"]["winner"],
                      "mitc4":data.get("mitc4")},indent=2))
