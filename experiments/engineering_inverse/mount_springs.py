"""Area-integrated elastic automotive brackets shared by both shell solvers.

Each bracket is idealized as two orthogonal spring-steel flexure leaves. The
assumed leaves are each 30 mm long, 30 mm wide and 1 mm thick; each direction
has transverse stiffness 3*E*I/L**3. This is an explicitly hypothetical
mounting adapter, not a production-vehicle drawing or calibrated stiffness.
Bilinear integration over a fixed flange area distributes the effective
spring energy, so total stiffness and centroid are mesh independent.
"""
from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp
import torch

from checkerboard_shell import BoundaryConditions, CheckerboardShell, NodalLoad, NewtonSettings
from experiments.engineering_inverse.automotive_panel import build_case, PARAMS
from experiments.cmame_round4.mitc4 import _ops

OUT=Path(__file__).resolve().parent/"results"
E_BRACKET=PARAMS["bracket_E_Pa"]
BRACKET_LENGTH=PARAMS["bracket_cantilever_length_m"]
BRACKET_WIDTH=PARAMS["bracket_cantilever_width_m"]
BRACKET_THICKNESS=PARAMS["bracket_cantilever_thickness_m"]
K_BRACKET=3*E_BRACKET*(BRACKET_WIDTH*BRACKET_THICKNESS**3/12)/BRACKET_LENGTH**3


def bracket_outer_fiber_stress(reaction):
    """Euler--Bernoulli root stress of either orthogonal leaf, in Pa."""
    force=np.asarray(reaction,dtype=float)
    return 6*BRACKET_LENGTH*np.abs(force[...,1:3]) / (
        BRACKET_WIDTH*BRACKET_THICKNESS**2)


def patch_weights(case):
    nx,ny=case["mesh"].nx,case["mesh"].ny
    width,height=case["width_m"],case["height_m"]
    sx,sy=width/nx,height/ny
    half_x=PARAMS["mount_half_width_m"]
    half_y=PARAMS["mount_half_height_m"]
    results=[]
    for fraction in PARAMS["mount_x_fractions"]:
        center_x=fraction*width
        center_y=PARAMS["mount_y_fraction"]*height
        weights=np.zeros(case["mesh"].num_vertices)
        for j in range(ny):
            y0,y1=j*sy,(j+1)*sy
            lo_y=max(y0,center_y-half_y);hi_y=min(y1,center_y+half_y)
            if hi_y<=lo_y:continue
            v0,v1=(lo_y-y0)/sy,(hi_y-y0)/sy
            iy0=(v1-v0)-.5*(v1*v1-v0*v0)
            iy1=.5*(v1*v1-v0*v0)
            for i in range(nx):
                x0,x1=i*sx,(i+1)*sx
                lo_x=max(x0,center_x-half_x);hi_x=min(x1,center_x+half_x)
                if hi_x<=lo_x:continue
                u0,u1=(lo_x-x0)/sx,(hi_x-x0)/sx
                ix0=(u1-u0)-.5*(u1*u1-u0*u0)
                ix1=.5*(u1*u1-u0*u0)
                for di,wi in ((0,ix0),(1,ix1)):
                    for dj,wj in ((0,iy0),(1,iy1)):
                        weights[(j+dj)*(nx+1)+i+di]+=sx*sy*wi*wj
        weights/=4*half_x*half_y
        if abs(weights.sum()-1)>1e-10:raise AssertionError("patch integration failed")
        results.append(weights)
    return np.stack(results)


def bracket_matrices(case, stiffness=K_BRACKET):
    W=patch_weights(case)
    ndof=3*case["mesh"].num_vertices
    D=np.zeros((ndof,8))
    diagonal=np.zeros(ndof)
    for p,weights in enumerate(W):
        for node,w in enumerate(weights):
            if w==0:continue
            for c in (1,2):
                dof=3*node+c
                value=stiffness*w
                diagonal[dof]+=value
                D[dof,2*p+c-1]+=value
    return sp.diags(diagonal,format="csr"),D,W


def checkerboard_forward(case, controls, stiffness=K_BRACKET,
                         tolerance=1e-6, initial=None,
                         hessian_assembly_strategy="cached"):
    wall_start=time.perf_counter()
    Kspring,D,W=bracket_matrices(case,stiffness)
    fixed=np.zeros(case["mesh"].shape,dtype=bool)
    fixed.reshape(-1)[[3*n+c for n in case["datum"] for c in range(3)]]=True
    fixed.reshape(-1)[case["gauge_x_dof"]]=True
    load=torch.as_tensor((D@np.asarray(controls)).reshape(case["mesh"].shape),dtype=torch.float64)
    model=CheckerboardShell(case["mesh"],case["material"],
        BoundaryConditions(torch.as_tensor(fixed)),NodalLoad(load),
        NewtonSettings(residual_tolerance=tolerance,maximum_newton_iterations=35),
        Kspring,hessian_assembly_strategy=hessian_assembly_strategy)
    model_construction_seconds=time.perf_counter()-wall_start
    solve_start=time.perf_counter()
    solved=model.solve(initial_displacement=initial)
    newton_iterations_total=int(solved.iterations)
    tangent_assemblies_total=len(solved.history)
    tangent_assembly_seconds=sum(h.get("assembly_seconds",0.) for h in solved.history)
    numeric_factorization_seconds=sum(h.get("factorization_seconds",0.) for h in solved.history)
    numeric_factorizations_total=sum("factorization_seconds" in h for h in solved.history)
    continuation_steps=0
    if solved.status!="CONVERGED":
        fraction=0.0;step=.25
        continuation_state=torch.zeros(case["mesh"].shape,dtype=torch.float64)
        while fraction<1.0-1e-12:
            trial_fraction=min(1.0,fraction+step)
            trial=model.solve(initial_displacement=continuation_state,
                              load_factor=trial_fraction)
            newton_iterations_total+=int(trial.iterations)
            tangent_assemblies_total+=len(trial.history)
            tangent_assembly_seconds+=sum(h.get("assembly_seconds",0.) for h in trial.history)
            numeric_factorization_seconds+=sum(h.get("factorization_seconds",0.) for h in trial.history)
            numeric_factorizations_total+=sum("factorization_seconds" in h for h in trial.history)
            continuation_steps+=1
            if trial.status=="CONVERGED":
                continuation_state=(trial.deformed_vertices-
                                    case["mesh"].reference_vertices).detach()
                solved=trial;fraction=trial_fraction;step=min(step*1.5,.5)
            else:
                step*=.5
                if step<1/128:
                    solved=trial
                    break
    nonlinear_analysis_seconds=time.perf_counter()-solve_start
    post_start=time.perf_counter()
    q=(solved.deformed_vertices-case["mesh"].reference_vertices).detach().numpy().reshape(-1,3)
    mount=np.zeros((4,3))
    xyz=case["mesh"].reference_vertices.detach().numpy().reshape(-1,3)
    nodal_mount=np.zeros_like(q)
    for p in range(4):
        for c in (1,2):
            mount[p,c]=stiffness*np.dot(W[p],float(controls[2*p+c-1])-q[:,c])
            nodal_mount[:,c]+=stiffness*W[p]*(float(controls[2*p+c-1])-q[:,c])
    residual,reaction=model.full_residual_and_reactions(solved.coordinates)
    datum=reaction.detach().numpy().reshape(-1,3)[case["datum"]]
    gauge_force=float(reaction.detach().numpy().reshape(-1)[case["gauge_x_dof"]])
    spring_energy=.5*sum(stiffness*np.dot(W[p],
        (q[:,c]-controls[2*p+c-1])**2)
        for p in range(4) for c in (1,2))
    return dict(status=solved.status,iterations=solved.iterations,
        newton_iterations_total=newton_iterations_total,
        tangent_assemblies_total=tangent_assemblies_total,
        tangent_assembly_seconds=tangent_assembly_seconds,
        numeric_factorization_seconds=numeric_factorization_seconds,
        numeric_factorizations_total=numeric_factorizations_total,
        final_shell_local_nnz=int(solved.history[-1]["local_tangent_nnz"]),
        rank_one_storage=int(solved.history[-1]["rank_one_storage"]),
        continuation_steps=continuation_steps,
        residual_inf=solved.residual_inf,displacement=q,
        mount_reaction=mount,datum_reaction=datum,
        gauge_reaction_N=gauge_force,
        shell_energy=solved.energies["total"]-float(.5*q.ravel()@(Kspring@q.ravel())),
        spring_energy=spring_energy,
        force_balance=mount.sum(axis=0)+datum.sum(axis=0)+
                      np.array([gauge_force,0.,0.]),
        reference_moment_balance=np.cross(xyz,nodal_mount).sum(axis=0)+
                                 np.cross(xyz[case["datum"]],datum).sum(axis=0)+
                                 np.cross(xyz[case["gauge_x_dof"]//3],
                                          [gauge_force,0.,0.]),
        current_moment_balance=np.cross(xyz+q,nodal_mount).sum(axis=0)+
                               np.cross(xyz[case["datum"]]+q[case["datum"]],datum).sum(axis=0)+
                               np.cross(xyz[case["gauge_x_dof"]//3]+q[case["gauge_x_dof"]//3],
                                        [gauge_force,0.,0.]),
        model=model,coordinates=solved.coordinates,
        spring_matrix=Kspring,design_force=D,weights=W,
        model_construction_seconds=model_construction_seconds,
        nonlinear_analysis_seconds=nonlinear_analysis_seconds,
        extraction_seconds=time.perf_counter()-post_start)


def mitc4_forward(case, controls, stiffness=K_BRACKET,
                  corotational=False, asd_linear=False, analysis_steps=None,
                  norm_type=2,residual_tolerance=1e-6,
                  gauge_x_fix=False):
    wall_start=time.perf_counter()
    ops=_ops();ops.wipe()
    ops.model("basic","-ndm",3,"-ndf",6)
    xyz=case["mesh"].reference_vertices.detach().numpy().reshape(-1,3)
    n=len(xyz)
    for node,point in enumerate(xyz,1):ops.node(node,*map(float,point))
    ops.section("ElasticMembranePlateSection",1,case["material"].young,
                case["material"].poisson,case["material"].thickness,0.0)
    for tag,quad in enumerate(case["mesh"].quads,1):
        args=("ASDShellQ4" if (corotational or asd_linear) else "ShellMITC4",tag,
              *[int(v)+1 for v in quad],1)
        ops.element(*args,*(["-corotational"] if corotational else []))
    for node in case["datum"]:ops.fix(int(node)+1,1,1,1,0,0,0)
    if gauge_x_fix:
        ops.fix(case["gauge_x_dof"]//3+1,1,0,0,0,0,0)
    W=patch_weights(case)
    ops.timeSeries("Linear",1);ops.pattern("Plain",1,1)
    anchor_by_patch=[[] for _ in range(4)]
    anchor_reference={}
    next_node=n+1;next_tag=len(case["mesh"].quads)+1
    for p,weights in enumerate(W):
        for node,w in enumerate(weights):
            if w<=1e-14:continue
            anchor=next_node;next_node+=1
            ops.node(anchor,*map(float,xyz[node]))
            anchor_reference[anchor]=xyz[node]
            ops.fix(anchor,1,0,0,1,1,1)
            ops.sp(anchor,2,float(controls[2*p]))
            ops.sp(anchor,3,float(controls[2*p+1]))
            ops.uniaxialMaterial("Elastic",next_tag,float(stiffness*w))
            mat_y=next_tag;next_tag+=1
            ops.uniaxialMaterial("Elastic",next_tag,float(stiffness*w))
            mat_z=next_tag;next_tag+=1
            ops.element("zeroLength",next_tag,int(node)+1,anchor,
                        "-mat",mat_y,mat_z,"-dir",2,3)
            next_tag+=1
            anchor_by_patch[p].append(anchor)
    ops.system("SparseGeneral");ops.numberer("RCM")
    steps=(16 if corotational else 1) if analysis_steps is None else int(analysis_steps)
    ops.constraints("Transformation");ops.integrator("LoadControl",1.0/steps)
    ops.test("NormUnbalance",float(residual_tolerance),80,0,int(norm_type))
    ops.algorithm("NewtonLineSearch" if corotational else "Linear")
    ops.analysis("Static")
    construction_seconds=time.perf_counter()-wall_start
    solve_start=time.perf_counter()
    status=0
    newton_iterations_total=0
    completed_steps=0
    for _ in range(steps):
        status=ops.analyze(1)
        newton_iterations_total+=int(ops.testIter())
        completed_steps+=1
        if status!=0:break
    norms=ops.testNorm()
    last_convergence_norm=float(norms[-1] if isinstance(norms,(list,tuple)) and norms
                                else norms if not isinstance(norms,(list,tuple)) else np.nan)
    analysis_seconds=time.perf_counter()-solve_start
    extraction_start=time.perf_counter()
    q=np.array([ops.nodeDisp(i+1,c+1) for i in range(n) for c in range(3)]).reshape(-1,3)
    ops.reactions()
    mount=np.zeros((4,3))
    for p,anchors in enumerate(anchor_by_patch):
        for anchor in anchors:
            for c in (1,2):mount[p,c]+=ops.nodeReaction(anchor,c+1)
    datum=np.array([[ops.nodeReaction(int(node)+1,c+1) for c in range(3)]
                    for node in case["datum"]])
    gauge_reaction=(float(ops.nodeReaction(case["gauge_x_dof"]//3+1,1))
                    if gauge_x_fix else 0.)
    moment=np.zeros(3)
    for p,anchors in enumerate(anchor_by_patch):
        for anchor in anchors:
            force=np.array([ops.nodeReaction(anchor,c+1) for c in range(3)])
            moment+=np.cross(anchor_reference[anchor],force)
    moment+=np.cross(xyz[case["datum"]],datum).sum(axis=0)
    moment+=np.cross(xyz[case["gauge_x_dof"]//3],[gauge_reaction,0.,0.])
    current_moment=np.zeros(3)
    for p,anchors in enumerate(anchor_by_patch):
        for anchor in anchors:
            force=np.array([ops.nodeReaction(anchor,c+1) for c in range(3)])
            position=anchor_reference[anchor].copy()
            position[1]+=float(controls[2*p]);position[2]+=float(controls[2*p+1])
            current_moment+=np.cross(position,force)
    current_moment+=np.cross(xyz[case["datum"]]+q[case["datum"]],datum).sum(axis=0)
    current_moment+=np.cross(xyz[case["gauge_x_dof"]//3]+q[case["gauge_x_dof"]//3],
                             [gauge_reaction,0.,0.])
    spring_energy=.5*sum(stiffness*np.dot(W[p],
        (q[:,c]-controls[2*p+c-1])**2)
        for p in range(4) for c in (1,2))
    return dict(status=status,displacement=q,mount_reaction=mount,
                datum_reaction=datum,gauge_reaction_N=gauge_reaction,
                spring_energy=spring_energy,
                force_balance=mount.sum(axis=0)+datum.sum(axis=0)+
                              np.array([gauge_reaction,0.,0.]),
                reference_moment_balance=moment,
                current_moment_balance=current_moment,
                anchor_counts=[len(a) for a in anchor_by_patch],
                model_construction_seconds=construction_seconds,
                analysis_seconds=analysis_seconds,
                extraction_seconds=time.perf_counter()-extraction_start,
                newton_iterations_total=newton_iterations_total,
                analysis_steps=steps,completed_analysis_steps=completed_steps,
                last_convergence_norm=last_convergence_norm,
                convergence_norm_type=int(norm_type),
                residual_tolerance=float(residual_tolerance),
                system_size=int(ops.systemSize()))


def audit():
    report={"stiffness_N_per_m":K_BRACKET,"bracket_geometry_m":{
        "length":BRACKET_LENGTH,"width":BRACKET_WIDTH,
        "thickness":BRACKET_THICKNESS},"meshes":{}}
    for nx,ny in ((12,8),(16,12),(20,16)):
        case=build_case(nx,ny)
        W=patch_weights(case)
        material_xy=case["mesh"].reference_vertices.detach().numpy().reshape(-1,3)[:,:2]
        rows={}
        for name,d in {"zero":np.zeros(8),
                       "one_z":np.array([0,1e-4,0,0,0,0,0,0]),
                       "opposite_z":np.array([0,1e-4,0,-1e-4,0,0,0,0]),
                       "four_z":np.tile([0,1e-4],4)}.items():
            cb=checkerboard_forward(case,d)
            mt=mitc4_forward(case,d)
            rows[name]={"input_m":d.tolist(),
                "checkerboard":{"status":cb["status"],"residual_inf":cb["residual_inf"],
                    "reaction_N":cb["mount_reaction"].tolist(),
                    "energy_J":cb["shell_energy"]+cb["spring_energy"],
                    "balance_N":cb["force_balance"].tolist()},
                "mitc4":{"status":mt["status"],"reaction_N":mt["mount_reaction"].tolist(),
                    "spring_energy_J":mt["spring_energy"],"balance_N":mt["force_balance"].tolist()}}
            print(nx,ny,name,np.max(np.linalg.norm(cb["mount_reaction"],axis=1)),
                  np.max(np.linalg.norm(mt["mount_reaction"],axis=1)),flush=True)
        report["meshes"][f"{nx}x{ny}"]={"weight_sums":W.sum(axis=1).tolist(),
            "weight_centers_xy":(W@material_xy).tolist(),"cases":rows}
        (OUT/"automotive_spring_audit.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


if __name__=="__main__":audit()
