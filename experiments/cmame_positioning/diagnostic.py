"""Compare CB and ShellMITC4 on the same Pq midpoint geometry.

This is a small *diagnostic* nonuniform-curvature patch, not a published
Checkerboard design object or a benchmark for shell-element accuracy.
Both methods receive the same physical midpoint coordinates, translation
supports and nodal loads. MITC4 is allowed to inherit the same q->Pq design
map; it is not charged for a fictional re-meshing operation.
"""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path
import time

import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import torch

from checkerboard_shell.geometry.core import torch_geometry
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent


HERE = Path(__file__).resolve().parent
RESULTS = HERE / "results"


def grid(n: int, parameter: float) -> torch.Tensor:
    x = torch.linspace(0, 1, n+1, dtype=torch.float64)
    yy, xx = torch.meshgrid(x, x, indexing="ij")
    zz = (0.12*torch.sin(math.pi*xx)*torch.sin(math.pi*yy)
          + 0.035*torch.sin(3*math.pi*xx)*torch.sin(2*math.pi*yy)
          + parameter*0.07*torch.sin(2*math.pi*xx)*torch.sin(math.pi*yy))
    return torch.stack((xx, yy, zz), -1)


def midpoint_map(n: int):
    vertex = lambda j,i: j*(n+1)+i
    nh = n*(n+1)
    horizontal = lambda j,i: j*n+i
    vertical = lambda j,i: nh+j*(n+1)+i
    edges = []
    for j in range(n+1):
        for i in range(n):
            edges.append((vertex(j,i),vertex(j,i+1)))
    for j in range(n):
        for i in range(n+1):
            edges.append((vertex(j,i),vertex(j+1,i)))
    rows = np.repeat(np.arange(len(edges)),2)
    cols = np.asarray(edges).reshape(-1)
    scalar = sp.coo_matrix((np.full(len(rows),.5),(rows,cols)),
                           shape=(len(edges),(n+1)**2)).tocsr()
    P = sp.kron(scalar,sp.eye(3),format="csr")
    faces=[]
    for j in range(n):
        for i in range(n):
            faces.append((horizontal(j,i),vertical(j,i+1),
                          horizontal(j+1,i),vertical(j,i)))
    for j in range(1,n):
        for i in range(1,n):
            faces.append((vertical(j-1,i),horizontal(j,i),
                          vertical(j,i),horizontal(j,i-1)))
    return P, np.asarray(faces,dtype=int), scalar, (horizontal,vertical)


def physical_fields(reference: torch.Tensor, P: sp.csr_matrix, faces, n: int):
    y=(P@reference.numpy().reshape(-1)).reshape(-1,3)
    first, normal = observed_first_normal(y, faces, n)
    return y, first, normal


def observed_first_normal(y: np.ndarray, faces, n: int):
    b=y[faces[:n*n]]
    ell=math.sqrt(2)/n
    t1=2*(b[:,1]-b[:,0])/ell
    t2=2*(b[:,3]-b[:,0])/ell
    first=np.stack((np.stack((np.sum(t1*t1,axis=1),np.sum(t1*t2,axis=1)),-1),
                    np.stack((np.sum(t2*t1,axis=1),np.sum(t2*t2,axis=1)),-1)),-2)
    normal=np.cross(t1,t2); normal/=np.linalg.norm(normal,axis=-1)[:,None]
    return first, normal


def observed_change(y0: np.ndarray, y1: np.ndarray, faces, n: int):
    """B/W observation changes; W on a warped B net is only a diagnostic."""
    a0,_=observed_first_normal(y0,faces,n)
    a1,_=observed_first_normal(y1,faces,n)
    b0=midpoint_second_form(y0,faces,n)
    b1=midpoint_second_form(y1,faces,n)
    quad=y1[faces[:n*n]]
    closure=quad[:,0]-quad[:,1]+quad[:,2]-quad[:,3]
    return {"B_metric_change_rms":float(np.sqrt(np.mean((a1-a0)**2))),
            "W_observed_curvature_change_rms":float(np.sqrt(np.mean((b1-b0)**2))),
            "B_parallelogram_closure_max_m":float(np.max(np.linalg.norm(closure,axis=1)))}


def midpoint_second_form(y: np.ndarray, faces: np.ndarray, n: int):
    """Compute the W observation directly from the shared midpoint geometry."""
    b=y[faces[:n*n]]
    ell=math.sqrt(2)/n
    t1=(2*(b[:,1]-b[:,0])/ell).reshape(n,n,3)
    t2=(2*(b[:,3]-b[:,0])/ell).reshape(n,n,3)
    normal=np.cross(t1,t2);normal/=np.linalg.norm(normal,axis=-1,keepdims=True)
    d1=(normal[1:,1:]-normal[:-1,:-1])/ell
    d2=(normal[1:,:-1]-normal[:-1,1:])/ell
    t1w=.5*(t1[1:,1:]+t1[:-1,:-1])
    t2w=.5*(t2[1:,:-1]+t2[:-1,1:])
    b11=-np.sum(d1*t1w,axis=-1)
    b12=-np.sum(d1*t2w,axis=-1)
    b21=-np.sum(d2*t1w,axis=-1)
    b22=-np.sum(d2*t2w,axis=-1)
    return np.stack((np.stack((b11,.5*(b12+b21)),-1),
                     np.stack((.5*(b12+b21),b22),-1)),-2)


def loaded_nodes(y: np.ndarray, n: int):
    support=np.flatnonzero(y[:,0]<=0.25+1e-12)
    loaded=np.flatnonzero(y[:,0]>=0.75-1e-12)
    force=np.zeros_like(y)
    force[loaded,2]=-0.001/len(loaded)
    return support,loaded,force


def cb_solve(reference, P, faces, n, support, loaded, force):
    targets=make_targets(reference,1/n,1/n)
    material=IsotropicMaterial(2e6,.3,.008)
    start=time.perf_counter()
    tangent,_,audit=assemble_paper_low_rank_tangent(
        torch.zeros_like(reference),reference,1/n,1/n,targets,material)
    assembly_s=time.perf_counter()-start
    fixed=np.asarray([3*int(k)+c for k in support for c in range(3)])
    C0=P[fixed].toarray()
    # Use exact row-space reduction; dependent support constraints must not
    # be passed into a Lagrange-multiplier solve.
    Q,R,pivot=la.qr(C0.T,pivoting=True,mode="economic")
    rank=int(np.sum(np.abs(np.diag(R))>1e-10))
    C=C0[pivot[:rank]]
    parity=np.fromfunction(lambda j,i:(-1.)**(i+j),(n+1,n+1)).reshape(-1)
    gauge=np.zeros((3,3*(n+1)**2))
    for c in range(3):gauge[c,c::3]=parity/len(parity)
    C=sp.csr_matrix(np.vstack((C,gauge)))
    rhs=P.T@force.reshape(-1)
    start=time.perf_counter()
    u,audit_solve=tangent.solve(rhs,constraints=C)
    solve_s=time.perf_counter()-start
    midpoint_u=(P@u).reshape(-1,3)
    K_u=tangent.action(u)
    free_residual=K_u-rhs
    # Least-squares reaction in the independent support/gauge basis.
    reaction_basis=la.solve((C@C.T).toarray(),
                            np.asarray(C@(-free_residual)),assume_a="sym")
    physical_reaction=np.zeros_like(force.reshape(-1))
    # Reduced C rows represent the selected original physical constraints.
    physical_reaction[fixed[pivot[:rank]]]=-reaction_basis[:rank]
    y=(P@reference.numpy().reshape(-1)).reshape(-1,3)
    current=reference+torch.as_tensor(u.reshape(reference.shape),dtype=reference.dtype)
    components=paper_components(current,reference,1/n,1/n,targets,material)
    internal={key:float(components[key].detach()) for key in
              ("membrane_b","membrane_xm","membrane","bending",
               "bending_compatibility","total")}
    total_force=force+physical_reaction.reshape(-1,3)
    return {"tip_z_m":float(midpoint_u[loaded,2].mean()),
            "strain_energy_J":float(.5*rhs@u),
            "energy_components_J":internal,
            **observed_change(y,y+midpoint_u,faces,n),
            "reaction_z_N":float(physical_reaction.reshape(-1,3)[:,2].sum()),
            "loaded_z_N":float(force[:,2].sum()),
            "force_balance_inf_N":float(np.max(np.abs(total_force.sum(axis=0)))),
            "moment_balance_inf_Nm":float(np.max(np.abs(np.cross(y,total_force).sum(axis=0)))),
            "free_residual_inf_N":float(np.max(np.abs(free_residual+C.T@reaction_basis))),
            "raw_states":reference.numel(),
            "active_physical_state_dim":reference.numel()-3,
            "independent_support_rows":rank,
            "gauge_rows":3,
            "tangent_nnz":int(tangent.local.nnz),
            "assembly_s":assembly_s,"solve_s":solve_s,
            "total_s":assembly_s+solve_s}


def mitc4_solve(y,faces,support,loaded,force,element="ShellMITC4"):
    import openseespy.opensees as ops
    n=int(round((1+math.sqrt(2*len(faces)-1))/2))
    start=time.perf_counter()
    ops.wipe();ops.model("basic","-ndm",3,"-ndf",6)
    for k,(x,yy,z) in enumerate(y,1):ops.node(k,float(x),float(yy),float(z))
    ops.section("ElasticMembranePlateSection",1,2e6,.3,.008,0.0)
    for k,face in enumerate(faces,1):
        ops.element(element,k,*[int(v)+1 for v in face],1)
    for k in support:ops.fix(int(k)+1,1,1,1,0,0,0)
    ops.timeSeries("Linear",1);ops.pattern("Plain",1,1)
    for k in loaded:
        ops.load(int(k)+1,*[float(v) for v in force[k]],0.,0.,0.)
    setup_s=time.perf_counter()-start
    ops.constraints("Plain");ops.numberer("RCM");ops.system("SparseGeneral")
    ops.test("NormUnbalance",1e-12,20);ops.algorithm("Linear")
    ops.integrator("LoadControl",1.0);ops.analysis("Static")
    start=time.perf_counter()
    status=ops.analyze(1)
    solve_s=time.perf_counter()-start
    if status:raise RuntimeError(f"MITC4 analyze returned {status}")
    displacement=np.asarray([ops.nodeDisp(k+1)[:3] for k in range(len(y))])
    ops.reactions()
    reactions=np.asarray([ops.nodeReaction(int(k)+1)[:3] for k in support])
    total_force=force.copy();total_force[support]+=reactions
    return {"tip_z_m":float(displacement[loaded,2].mean()),
            "strain_energy_J":float(.5*np.sum(force*displacement)),
            **observed_change(y,y+displacement,faces,n),
            "reaction_z_N":float(reactions[:,2].sum()),
            "loaded_z_N":float(force[:,2].sum()),
            "force_balance_inf_N":float(np.max(np.abs(total_force.sum(axis=0)))),
            "moment_balance_inf_Nm":float(np.max(np.abs(np.cross(y,total_force).sum(axis=0)))),
            "raw_states":int(6*len(y)),
            "active_states":int(6*len(y)-3*len(support)),
            "support_nodes":len(support),
            "model_setup_s":setup_s,"solve_s":solve_s,
            "total_s":setup_s+solve_s}


def evaluate(n:int, parameter:float):
    q=grid(n,parameter)
    P,faces,scalar,_=midpoint_map(n)
    start=time.perf_counter()
    y,first,normal=physical_fields(q,P,faces,n)
    map_s=time.perf_counter()-start
    cb_geometry=torch_geometry(q,1/n,1/n)
    metric_delta=float(np.max(np.abs(first-cb_geometry.first_form.numpy().reshape(-1,2,2))))
    normal_delta=float(np.max(np.abs(normal-cb_geometry.normals_b.numpy().reshape(-1,3))))
    curvature_delta=float(np.max(np.abs(midpoint_second_form(y,faces,n)-
                                          cb_geometry.second_form.numpy())))
    support,loaded,force=loaded_nodes(y,n)
    cb=cb_solve(q,P,faces,n,support,loaded,force)
    mitc=mitc4_solve(y,faces,support,loaded,force)
    dq=(grid(n,parameter+1e-5)-grid(n,parameter-1e-5)).numpy().reshape(-1)/2e-5
    dy_fd=((P@grid(n,parameter+1e-5).numpy().reshape(-1))-
           (P@grid(n,parameter-1e-5).numpy().reshape(-1)))/2e-5
    chain_delta=float(np.max(np.abs(P@dq-dy_fd)))
    gauge=np.fromfunction(lambda j,i:(-1.)**(i+j),(n+1,n+1)).reshape(-1)
    g=np.repeat(gauge[:,None],3,axis=1).reshape(-1)
    kernel_residual=float(np.max(np.abs(P@g)))
    # The dense-rank check is limited to this small diagnostic model.
    rank=int(np.linalg.matrix_rank(P.toarray(),tol=1e-10))
    return {"n":n,"parameter":parameter,"surface":"synthetic nonuniform-curvature diagnostic",
            "raw_vertices":(n+1)**2,"physical_midpoint_nodes":len(y),
            "P_shape":list(P.shape),"P_rank":rank,
            "P_kernel_dim":P.shape[1]-rank,
            "independent_midpoint_compatibility_constraints":P.shape[0]-rank,
            "gauge_P_inf":kernel_residual,"B_metric_max_abs_difference":metric_delta,
            "B_normal_max_abs_difference":normal_delta,
            "W_second_form_max_abs_difference":curvature_delta,
            "design_chain_rule_max_abs_difference":chain_delta,
            "geometry_map_s":map_s,"CB":cb,"MITC4_Pq":mitc}


def main():
    RESULTS.mkdir(exist_ok=True)
    rows=[]
    for n in (8,12):
        for parameter in (-.1,0.,.1):
            row=evaluate(n,parameter)
            rows.append(row)
            (RESULTS/"diagnostic.json").write_text(json.dumps(rows,indent=2)+"\n")
            print(n,parameter,row["CB"]["tip_z_m"],row["MITC4_Pq"]["tip_z_m"],flush=True)
    flat=[]
    for r in rows:
        for method in ("CB","MITC4_Pq"):
            flat.append({k:r[k] for k in ("n","parameter","surface","P_rank","P_kernel_dim",
                                             "independent_midpoint_compatibility_constraints",
                                             "B_metric_max_abs_difference","B_normal_max_abs_difference",
                                             "W_second_form_max_abs_difference",
                                             "design_chain_rule_max_abs_difference","geometry_map_s")}
                        |{"method":method,**r[method]})
    with (RESULTS/"diagnostic.csv").open("w",newline="") as f:
        writer=csv.DictWriter(f,fieldnames=list(dict.fromkeys(k for row in flat for k in row)))
        writer.writeheader();writer.writerows(flat)
    for n in (8,12):
        subset=[r for r in rows if r["n"]==n]
        for method in ("CB","MITC4_Pq"):
            sensitivity=(subset[2][method]["tip_z_m"]-subset[0][method]["tip_z_m"])/.2
            print("design sensitivity",n,method,sensitivity,flush=True)


if __name__=="__main__":main()
