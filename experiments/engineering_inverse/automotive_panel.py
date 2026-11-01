"""Representative front-fender patch: elastic mounting-point adjustment.

The surface is parametric, not a production vehicle CAD or measured stamping.
All coordinates use metres, forces newtons, and Young's modulus pascals.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import scipy.optimize as opt
import torch

from checkerboard_shell import (BoundaryConditions, CheckerboardMesh,
                                CheckerboardShell, IsotropicMaterial,
                                NewtonSettings)
from experiments.engineering_inverse.common import (
    checkerboard_displacement_operator, mitc4_prescribed_response)

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"
OUT.mkdir(exist_ok=True)
PARAMS=json.loads((HERE/"config.json").read_text())["automotive"]


def build_case(nx=12, ny=8):
    width, height = PARAMS["width_m"],PARAMS["height_m"]
    base = CheckerboardMesh.rectangle(width, height, nx, ny)
    nominal = base.reference_vertices.detach().numpy().copy()
    xi = nominal[..., 0] / width
    eta = nominal[..., 1] / height
    # Curved wheel-arch side; the top two rows stand for a continuous flange.
    nominal[..., 1] += .035 * np.sin(np.pi*xi) * eta**2
    nominal[..., 2] = .026*np.sin(np.pi*xi)*np.sin(.8*np.pi*eta) + .023*eta**2
    manufactured = nominal.copy()
    # Prescribed smooth manufacturing error, not scan or production data.
    manufactured[..., 1] += PARAMS["manufacturing_y_amplitude_m"]*np.sin(np.pi*xi)*eta**2
    manufactured[..., 2] += (PARAMS["manufacturing_z_sine_amplitude_m"]*np.sin(2*np.pi*xi)
                             + PARAMS["manufacturing_z_cosine_amplitude_m"]*np.cos(np.pi*xi))*eta**2
    mesh = CheckerboardMesh(torch.as_tensor(manufactured), base.hx, base.hy,
                            base.quads, base.boundaries)
    node = lambda i, j: j*(nx+1)+i
    datum = [node(0, 0), node(nx, 0)]
    # Both datum nodes lie on the same checkerboard sublattice for even nx.
    # Fix one x-coordinate on the other sublattice to remove the exact
    # relative-sublattice translation gauge of the position-only mesh.
    gauge_x_dof = 3*node(1, 0)
    # Four wheel-arch flange brackets at fixed physical fractions. The
    # geometry has no CAD mounting holes; the selected nodes discretize them.
    mount_j=int(round(PARAMS["mount_y_fraction"]*ny))
    mount_fractions=tuple(PARAMS["mount_x_fractions"])
    mounts = [node(int(round(f*nx)),mount_j) for f in mount_fractions]
    mount_patches=[]
    for fraction in mount_fractions:
        x_center=fraction*width; y_center=PARAMS["mount_y_fraction"]*height
        patch=[node(i,j) for j in range(ny+1) for i in range(nx+1)
               if abs(i*width/nx-x_center)<=PARAMS["mount_half_width_m"]+1e-12
               and abs(j*height/ny-y_center)<=PARAMS["mount_half_height_m"]+1e-12]
        if not patch:raise ValueError("mount patch has no nodes")
        mount_patches.append(patch)
    fixed_nodes=sorted(set(datum+[n for patch in mount_patches for n in patch]))
    fixed = [3*n+c for n in fixed_nodes for c in range(3)]
    design = [[3*n+c for n in patch] for patch in mount_patches for c in (1,2)]
    obs = np.zeros((2*(nx+1), 3*mesh.num_vertices))
    geom_error = np.zeros(2*(nx+1))
    directions = []
    for i in range(nx+1):
        ip, im = min(i+1,nx), max(i-1,0)
        tangent = nominal[ny, ip]-nominal[ny, im]
        tangent /= np.linalg.norm(tangent)
        across = nominal[ny]-nominal[ny-1]
        across = across[i]; across /= np.linalg.norm(across)
        normal = np.cross(tangent, across); normal /= np.linalg.norm(normal)
        gap_direction = np.cross(normal, tangent)
        directions.append([gap_direction.tolist(), normal.tolist()])
        for component, direction in enumerate((gap_direction, normal)):
            row = 2*i+component
            obs[row, 3*node(i,ny):3*node(i,ny)+3] = direction
            geom_error[row] = direction @ (manufactured[ny,i]-nominal[ny,i])
    return dict(mesh=mesh, nominal=nominal, manufactured=manufactured,
                fixed=fixed, design=design, mounts=mounts,
                mount_patches=mount_patches,datum=datum,gauge_x_dof=gauge_x_dof,
                O=obs, geometry_error=geom_error, directions=directions,
                material=IsotropicMaterial(PARAMS["steel_E_Pa"],
                    PARAMS["poisson"],PARAMS["thickness_m"]),
                target_gap_m=PARAMS["target_gap_m"],
                width_m=width, height_m=height)


def optimize_map(response0, response_map, geometry_error, bounds=None,
                 reaction0=None, reaction_map=None, reaction_limit=None,
                 reaction_groups=None):
    if bounds is None: bounds=PARAMS["adjustment_limit_m"]
    if reaction_limit is None: reaction_limit=PARAMS["reaction_limit_N"]
    scale = .001
    n = response_map.shape[1]
    weights = np.tile([8.0, 1.0], len(geometry_error)//2)
    def fun(x):
        error = (geometry_error+response0+response_map@(scale*x))/scale
        return float((weights*error)@error/len(error)+.0005*(x@x))
    def grad(x):
        error = (geometry_error+response0+response_map@(scale*x))/scale
        return 2*response_map.T@(weights*error)/len(error)+.001*x
    start = time.perf_counter()
    constraints=[]
    if reaction0 is not None:
        A=scale*reaction_map
        if reaction_groups is None:
            constraints=[dict(type="ineq",fun=lambda x: reaction_limit-reaction0-A@x,
                              jac=lambda x:-A),
                         dict(type="ineq",fun=lambda x: reaction_limit+reaction0+A@x,
                              jac=lambda x:A)]
        else:
            def resultant(x):
                r=reaction0+A@x
                return np.array([np.linalg.norm([r[component].sum()
                                      for component in group])
                                 for group in reaction_groups])
            def resultant_jac(x):
                r=reaction0+A@x
                rows=[]
                for group in reaction_groups:
                    vector=np.array([r[component].sum() for component in group])
                    derivative=np.array([A[component].sum(axis=0)
                                         for component in group])
                    rows.append(vector@derivative/max(np.linalg.norm(vector),1e-12))
                return np.array(rows)
            constraints=[dict(type="ineq",fun=lambda x: reaction_limit-resultant(x),
                              jac=lambda x:-resultant_jac(x))]
    solution = opt.minimize(fun, np.zeros(n), jac=grad,
                            method="SLSQP" if constraints else "L-BFGS-B",
                            bounds=[(-bounds/scale,bounds/scale)]*n,
                            constraints=constraints,
                            options=dict(ftol=1e-12, maxiter=300))
    solution.x *= scale
    return solution, time.perf_counter()-start


def metrics(case, q, error):
    gap = case["target_gap_m"] - error[0::2]
    flush = error[1::2]
    return dict(gap_mm=(gap*1e3).tolist(), flush_mm=(flush*1e3).tolist(),
                gap_error_rms_mm=float(np.sqrt(np.mean(error[0::2]**2))*1e3),
                flush_rms_mm=float(np.sqrt(np.mean(flush**2))*1e3),
                max_abs_error_mm=float(np.max(np.abs(error))*1e3),
                max_elastic_displacement_mm=float(np.max(np.linalg.norm(q.reshape(-1,3),axis=1))*1e3))


def membrane_principal_strain(case,q):
    reference=case["manufactured"]
    deformed=reference+q.reshape(reference.shape)
    values=[]
    for j in range(reference.shape[0]-1):
        row=[]
        for i in range(reference.shape[1]-1):
            def tangents(vertices):
                p=vertices[j:j+2,i:i+2]
                return np.column_stack((.5*(p[0,1]+p[1,1]-p[0,0]-p[1,0]),
                                        .5*(p[1,0]+p[1,1]-p[0,0]-p[0,1])))
            A=tangents(reference); a=tangents(deformed)
            G=A.T@A
            eig,V=np.linalg.eigh(G)
            invroot=(V/np.sqrt(eig))@V.T
            strain=.5*invroot@(a.T@a-G)@invroot
            row.append(float(np.max(np.abs(np.linalg.eigvalsh(strain)))))
        values.append(row)
    return values


def mount_reactions(case, reaction, fixed_lookup=None):
    """Sum nodal reaction vectors over each physical bracket footprint."""
    if fixed_lookup is None:
        nodal=np.asarray(reaction).reshape(-1,3)
        return np.array([nodal[patch].sum(axis=0)
                         for patch in case["mount_patches"]])
    return np.array([[sum(reaction[fixed_lookup[3*node+c]] for node in patch)
                      for c in range(3)] for patch in case["mount_patches"]])


def prescribed_field(case, controls):
    values=torch.zeros(case["mesh"].shape,dtype=torch.float64)
    flat=values.reshape(-1)
    for group,value in zip(case["design"],controls):
        flat[group]=float(value)
    return values


def run(nx=12, ny=8, with_mitc4=True, nonlinear_check=True):
    case = build_case(nx, ny)
    cb = checkerboard_displacement_operator(case["mesh"], case["material"],
        case["fixed"], case["design"], case["O"])
    fixed_lookup={int(dof):k for k,dof in enumerate(cb["fixed"])}
    cb_groups=[[[fixed_lookup[3*node+c] for node in patch]
                for c in range(3)]
               for patch in case["mount_patches"]]
    solution, opt_time = optimize_map(cb["response0"], cb["response_map"],
        case["geometry_error"],reaction0=cb["reaction0"],
        reaction_map=cb["reaction_map"],reaction_groups=cb_groups)
    d = solution.x
    q = cb["q0"]+cb["Q"]@d
    error0 = case["geometry_error"]+cb["response0"]
    error = error0+cb["response_map"]@d
    reaction = cb["reaction0"]+cb["reaction_map"]@d
    mount_reaction=mount_reactions(case,reaction,fixed_lookup)
    energy = .5*q@cb["K"]@q
    result = dict(method="Checkerboard + exact design derivative",
        free_dofs=cb["free_dofs"], total_dofs=cb["total_dofs"],
        assembly_seconds=cb["assembly_seconds"], optimization_seconds=opt_time,
        total_seconds=cb["assembly_seconds"]+opt_time,
        iterations=int(solution.nit), success=bool(solution.success),
        objective=float(solution.fun),
        adjustment_mm=(d*1e3).tolist(),
        maximum_mount_reaction_N=float(np.max(np.linalg.norm(mount_reaction,axis=1))),
        maximum_constraint_reaction_component_N=float(np.max(np.abs(reaction))),
        mount_reaction_N=mount_reaction.tolist(),
        strain_energy_J=float(energy),
        max_abs_membrane_principal_strain=float(np.max(membrane_principal_strain(case,q))),
        membrane_principal_strain=membrane_principal_strain(case,q),
        initial=metrics(case, cb["q0"], error0),
        final=metrics(case, q, error))
    payload = dict(case="representative_automotive_fender_patch",
        mesh=[nx,ny], material=dict(E_Pa=case["material"].young,
            nu=case["material"].poisson, thickness_m=case["material"].thickness),
        geometry_source="parametric research assumption",
        manufacturing_error_source="prescribed smooth research field",
        gap_target_mm=PARAMS["target_gap_m"]*1000,
        gap_definition="adjacent nominal edge offset along local in-surface edge normal",
        flush_definition="offset along nominal surface normal",
        mount_nodes=case["mounts"], datum_nodes=case["datum"],
        mount_patch_nodes=case["mount_patches"],
        mount_patch_footprint_m=[2*PARAMS["mount_half_width_m"],
                                 2*PARAMS["mount_half_height_m"]],
        permitted_adjustment="mount local/global y and z; x fixed; +/-1.5 mm",
        reaction_limit_N=PARAMS["reaction_limit_N"],
        nominal_vertices_m=case["nominal"].reshape(-1,3).tolist(),
        manufactured_vertices_m=case["manufactured"].reshape(-1,3).tolist(),
        checkerboard=result)
    # Best rigid motion compatible with the two upper-flange datums is a
    # rotation about their connecting line. This is a diagnostic baseline.
    xyz=case["manufactured"].reshape(-1,3)
    axis=xyz[case["datum"][1]]-xyz[case["datum"][0]]
    axis/=np.linalg.norm(axis)
    rigid=np.cross(np.broadcast_to(axis,xyz.shape),xyz-xyz[case["datum"][0]]).reshape(-1)
    rigid_response=case["O"]@rigid
    weights=np.tile([8.0,1.0],len(error0)//2)
    angle=-np.dot(weights*rigid_response,error0)/np.dot(weights*rigid_response,rigid_response)
    angle=float(np.clip(angle,-.0015/case["height_m"],.0015/case["height_m"]))
    rigid_q=angle*rigid
    payload["rigid_datum_baseline"]=dict(rotation_about_datum_axis_rad=angle,
        final=metrics(case,rigid_q,error0+angle*rigid_response))
    if nonlinear_check:
        mask=torch.zeros(case["mesh"].shape,dtype=torch.bool)
        mask.reshape(-1)[case["fixed"]]=True
        prescribed=prescribed_field(case,d)
        model=CheckerboardShell(case["mesh"],case["material"],
            BoundaryConditions(mask,prescribed),solver_settings=NewtonSettings(
                residual_tolerance=1e-6,maximum_newton_iterations=30))
        solved=model.solve()
        nq=(solved.deformed_vertices-case["mesh"].reference_vertices).numpy().reshape(-1)
        nerror=case["geometry_error"]+case["O"]@nq
        reaction_full=model.full_residual_and_reactions(solved.coordinates)[1].numpy().reshape(-1)
        mount_reaction_nl=mount_reactions(case,reaction_full)
        payload["checkerboard_nonlinear_check"]=dict(status=solved.status,
            iterations=solved.iterations,residual_inf=solved.residual_inf,
            final=metrics(case,nq,nerror),
            displacement_m=nq.tolist(),
            maximum_mount_reaction_N=float(np.max(np.linalg.norm(mount_reaction_nl,axis=1))),
            mount_reaction_N=mount_reaction_nl.tolist(),
            maximum_constraint_reaction_component_N=float(np.max(np.abs(reaction_full))),
            max_abs_membrane_principal_strain=float(np.max(membrane_principal_strain(case,nq))),
            strain_energy_J=solved.energies["total"])
        if (solved.status!="CONVERGED" or
            payload["checkerboard_nonlinear_check"]["maximum_mount_reaction_N"]
                >PARAMS["reaction_limit_N"]+1e-5):
            for factor in (.9,.8,.7,.6,.5):
                model.boundary_conditions=BoundaryConditions(mask,prescribed_field(case,d*factor))
                trial=model.solve()
                if trial.status!="CONVERGED":continue
                tq=(trial.deformed_vertices-case["mesh"].reference_vertices).numpy().reshape(-1)
                te=case["geometry_error"]+case["O"]@tq
                tr=mount_reactions(case,model.full_residual_and_reactions(
                    trial.coordinates)[1].numpy().reshape(-1))
                max_reaction=float(np.max(np.linalg.norm(tr,axis=1)))
                if max_reaction<=PARAMS["reaction_limit_N"]+1e-5:
                    payload["checkerboard_nonlinear_feasible"]=dict(
                        adjustment_mm=(d*factor*1e3).tolist(),
                        scale_of_linear_optimum=factor,status=trial.status,
                        iterations=trial.iterations,residual_inf=trial.residual_inf,
                        final=metrics(case,tq,te),
                        displacement_m=tq.tolist(),
                        maximum_mount_reaction_N=max_reaction,
                        mount_reaction_N=tr.tolist(),
                        max_abs_membrane_principal_strain=float(
                            np.max(membrane_principal_strain(case,tq))),
                        strain_energy_J=trial.energies["total"])
                    break
    if with_mitc4:
        start = time.perf_counter()
        zero = np.zeros(len(case["design"]))
        mt0, mq0, reaction0 = mitc4_prescribed_response(case["mesh"],case["material"],
            case["fixed"],case["design"],zero,case["O"])
        responses=[]; reactions=[]
        for i in range(len(zero)):
            response,_,reaction=mitc4_prescribed_response(
                case["mesh"],case["material"],case["fixed"],case["design"],
                np.eye(len(zero))[i]*.001,case["O"])
            responses.append((response-mt0)/.001)
            reactions.append(reaction)
        M=np.column_stack(responses)
        R=np.column_stack([(reaction-reaction0)/.001 for reaction in reactions])
        basis_time = time.perf_counter()-start
        mt_lookup={int(dof):k for k,dof in enumerate(case["fixed"])}
        mt_groups=[[[mt_lookup[3*node+c] for node in patch]
                    for c in range(3)]
                   for patch in case["mount_patches"]]
        ms, mo = optimize_map(mt0,M,case["geometry_error"],
                              reaction0=reaction0,reaction_map=R,
                              reaction_groups=mt_groups)
        mr, mq, mreact = mitc4_prescribed_response(case["mesh"],case["material"],
            case["fixed"],case["design"],ms.x,case["O"])
        merr = case["geometry_error"]+mr
        mt_mount=mount_reactions(case,mreact,mt_lookup)
        payload["mitc4"] = dict(method="OpenSees ShellMITC4 + sampled linear response",
            free_dofs=6*case["mesh"].num_vertices-len(case["fixed"]),
            basis_seconds=basis_time, optimization_seconds=mo,
            total_seconds=basis_time+mo, iterations=int(ms.nit),
            success=bool(ms.success), objective=float(ms.fun),
            optimizer_message=str(ms.message),
            adjustment_mm=(ms.x*1e3).tolist(),
            maximum_mount_reaction_N=float(np.max(np.linalg.norm(mt_mount,axis=1))),
            maximum_constraint_reaction_component_N=float(np.max(np.abs(mreact))),
            mount_reaction_N=mt_mount.tolist(),
            final=metrics(case,mq,merr))
    path = OUT/f"automotive_{nx}x{ny}.json"
    path.write_text(json.dumps(payload,indent=2)+"\n")
    return payload


if __name__ == "__main__":
    summary = run()
    print(json.dumps({k:summary[k] for k in ("checkerboard","mitc4")},indent=2))
