"""Small-displacement design operators built on the frozen shell tangent."""

from __future__ import annotations

import time
import numpy as np
import scipy.sparse.linalg as spla
import torch

from checkerboard_shell import (BoundaryConditions, CheckerboardMesh,
                                CheckerboardShell, NodalLoad)


def strip_mesh(parent, node_grid):
    """A bonded reinforcement strip sharing parent translational nodes."""
    nodes = np.asarray(node_grid, dtype=int)
    ny, nx = nodes.shape[0]-1, nodes.shape[1]-1
    if nx < 3 or ny < 3:
        raise ValueError("Checkerboard reinforcement needs at least 3x3 cells")
    quads = []
    for j in range(ny):
        for i in range(nx):
            a = j*(nx+1)+i
            quads.append((a,a+1,a+nx+2,a+nx+1))
    sub = CheckerboardMesh(parent.reference_vertices.reshape(-1,3)[nodes].clone(),
                           parent.hx, parent.hy, np.asarray(quads), {})
    return sub


def positional_beam_stiffness(vertices, node_line, young, area, inertia):
    """Small-strain axial and transverse Euler-Bernoulli line stiffness.

    Translational-only centerline discretization. Interior curvatures use a
    three-node second difference; the beam shares shell positions exactly.
    """
    xyz = np.asarray(vertices).reshape(-1,3)
    line = list(map(int,node_line))
    if len(line) < 3:
        raise ValueError("beam line requires at least three nodes")
    K = np.zeros((3*len(xyz), 3*len(xyz)))
    eye = np.eye(3)
    for a,b in zip(line[:-1],line[1:]):
        delta = xyz[b]-xyz[a]; length = np.linalg.norm(delta)
        tangent = delta/length
        block = young*area/length*np.outer(tangent,tangent)
        for u,v,sign in ((a,a,1),(b,b,1),(a,b,-1),(b,a,-1)):
            K[3*u:3*u+3,3*v:3*v+3] += sign*block
    for a,b,c in zip(line[:-2],line[1:-1],line[2:]):
        delta = xyz[c]-xyz[a]; tangent = delta/np.linalg.norm(delta)
        length = .5*(np.linalg.norm(xyz[b]-xyz[a])+
                     np.linalg.norm(xyz[c]-xyz[b]))
        projector = eye-np.outer(tangent,tangent)
        block = young*inertia/length**3*projector
        for u,wu in ((a,1),(b,-2),(c,1)):
            for v,wv in ((a,1),(b,-2),(c,1)):
                K[3*u:3*u+3,3*v:3*v+3] += wu*wv*block
    return K


def checkerboard_displacement_operator(mesh, material, fixed_dofs, design_dofs,
                                       observations, load=None,
                                       reinforcements=(), beam_lines=(),
                                       prepared=None):
    """Linearize at the manufactured reference and retain exact K_fc.

    ``observations`` is a dense matrix acting on all 3V displacement DOFs.
    Returns the response and reaction matrices for prescribed displacements.
    """
    start = time.perf_counter()
    shell_tangent_seconds=0.0
    local_nnz=None
    rank_one_storage=None
    if prepared is None:
        bc = BoundaryConditions(torch.zeros(mesh.shape, dtype=torch.bool))
        model = CheckerboardShell(mesh, material, bc, NodalLoad.zero(mesh) if load is None else load)
        zero = torch.zeros(3*mesh.num_vertices, dtype=torch.float64)
        tangent, residual, _ = model.low_rank_tangent(zero)
        shell_tangent_seconds=time.perf_counter()-start
        local_nnz=int(tangent.local.nnz)
        rank_one_storage=tangent.rank_one_storage
        K = tangent.local + tangent.coefficient * np.outer(tangent.vector, tangent.vector)
        K = np.asarray(K)
    else:
        K = prepared["K"]
        residual = prepared["residual"]
    all_dofs = np.arange(3*mesh.num_vertices)
    fixed = np.array(sorted(set(fixed_dofs)), dtype=int)
    free = np.setdiff1d(all_dofs, fixed)
    groups = [([int(d)] if np.isscalar(d) else list(map(int,d)))
              for d in design_dofs]
    design = np.array([d for group in groups for d in group], dtype=int)
    if not np.all(np.isin(design, fixed)):
        raise ValueError("design DOF must be constrained")
    for node_grid, extra_material in (() if prepared is not None else reinforcements):
        sub = strip_mesh(mesh, node_grid)
        sub_model = CheckerboardShell(
            sub, extra_material,
            BoundaryConditions(torch.zeros(sub.shape, dtype=torch.bool)),
            NodalLoad.zero(sub))
        sub_zero = torch.zeros(3*sub.num_vertices, dtype=torch.float64)
        sub_tangent, sub_residual, _ = sub_model.low_rank_tangent(sub_zero)
        sub_K = sub_tangent.local + sub_tangent.coefficient*np.outer(
            sub_tangent.vector, sub_tangent.vector)
        ids = np.array([3*n+c for n in np.asarray(node_grid).ravel()
                        for c in range(3)])
        K[np.ix_(ids,ids)] += np.asarray(sub_K)
        residual[ids] += sub_residual
    if prepared is None:
        for line, young, area, inertia in beam_lines:
            K += positional_beam_stiffness(mesh.reference_vertices.detach().numpy(),
                                           line, young, area, inertia)
    system_build_seconds=time.perf_counter()-start
    Kff = K[np.ix_(free, free)]
    factor_start=time.perf_counter()
    solve = spla.factorized(__import__('scipy').sparse.csc_matrix(Kff))
    factor_seconds=time.perf_counter()-factor_start
    force = np.zeros(len(all_dofs)) if load is None else load.values.detach().numpy().reshape(-1)
    baseline = solve(-residual[free])
    controls = np.zeros((len(all_dofs), len(groups)))
    for col, group in enumerate(groups):
        controls[group, col] = 1.0
    H = np.empty((len(free), len(groups)))
    for col in range(len(groups)):
        H[:, col] = solve(-K[np.ix_(free, groups[col])].sum(axis=1))
    q0 = np.zeros(len(all_dofs)); q0[free] = baseline
    Q = controls.copy(); Q[free, :] = H
    reaction0 = (K @ q0 - force)[fixed]
    reaction_map = (K[fixed, :] @ Q)
    return {
        "q0": q0, "Q": Q, "response0": observations @ q0,
        "response_map": observations @ Q,
        "reaction0": reaction0, "reaction_map": reaction_map,
        "fixed": fixed, "free": free, "K": K,
        "residual": residual,
        "factorized": solve,
        "assembly_seconds": time.perf_counter()-start,
        "shell_tangent_seconds": shell_tangent_seconds,
        "system_build_seconds": system_build_seconds,
        "factor_seconds": factor_seconds,
        "shell_local_nnz": local_nnz,
        "rank_one_storage": rank_one_storage,
        "free_dofs": len(free), "total_dofs": len(all_dofs),
    }


def mitc4_prescribed_response(mesh, material, fixed_dofs, design_dofs, values,
                              observation_matrix, load=None,
                              reinforcements=(), beam_lines=()):
    """OpenSees ShellMITC4 with the same translational mounts and geometry."""
    from experiments.cmame_round4.mitc4 import _ops
    ops = _ops(); ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    xyz = mesh.reference_vertices.detach().numpy().reshape(-1, 3)
    for node, point in enumerate(xyz, 1):
        ops.node(node, *map(float, point))
    ops.section("ElasticMembranePlateSection", 1, material.young,
                material.poisson, material.thickness, 0.0)
    for tag, quad in enumerate(mesh.quads, 1):
        ops.element("ShellMITC4", tag, *map(lambda v: int(v)+1, quad), 1)
    tag = len(mesh.quads)+1
    for section, (node_grid, extra_material) in enumerate(reinforcements, 2):
        ops.section("ElasticMembranePlateSection", section, extra_material.young,
                    extra_material.poisson, extra_material.thickness, 0.0)
        grid = np.asarray(node_grid, dtype=int)
        for j in range(grid.shape[0]-1):
            for i in range(grid.shape[1]-1):
                quad = (grid[j,i],grid[j,i+1],grid[j+1,i+1],grid[j+1,i])
                ops.element("ShellMITC4", tag,
                            *[int(v)+1 for v in quad], section)
                tag += 1
    # The beam baseline uses ordinary 3D Euler-Bernoulli frame elements.
    # Rotations are shared with the MITC4 shell nodes; the Checkerboard line
    # uses slopes derived from positions. Their continuum EA and EI match.
    for line_number, (line, young, area, inertia) in enumerate(beam_lines, 1):
        ops.geomTransf("Linear", line_number, 0.0, 0.0, 1.0)
        shear = young/(2.0*(1.0+0.3))
        J = 2.0*inertia
        for a,b in zip(line[:-1],line[1:]):
            ops.element("elasticBeamColumn", tag, int(a)+1, int(b)+1,
                        area, young, shear, J, inertia, inertia, line_number)
            tag += 1
    groups = [([int(d)] if np.isscalar(d) else list(map(int,d)))
              for d in design_dofs]
    design_lookup = {int(d): float(v) for group, v in zip(groups, values)
                     for d in group}
    fixed_lookup = set(map(int, fixed_dofs))
    for node in range(mesh.num_vertices):
        fixity = [int((3*node+c) in fixed_lookup and
                      (3*node+c) not in design_lookup) for c in range(3)]
        if any(fixity):
            ops.fix(node+1, *fixity, 0, 0, 0)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for dof, value in design_lookup.items():
        ops.sp(dof//3+1, dof%3+1, value)
    if load is not None:
        forces = load.values.detach().numpy().reshape(-1,3)
        for node, force in enumerate(forces, 1):
            if np.any(force):
                ops.load(node, *map(float, force), 0.0, 0.0, 0.0)
    ops.system("SparseGeneral")
    ops.numberer("RCM")
    ops.constraints("Transformation")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    if ops.analyze(1) != 0:
        raise RuntimeError("MITC4 prescribed-displacement analysis failed")
    q = np.array([ops.nodeDisp(n+1, c+1)
                  for n in range(mesh.num_vertices) for c in range(3)])
    ops.reactions()
    reaction = np.array([ops.nodeReaction(int(d//3+1), int(d%3+1))
                         for d in fixed_dofs])
    return observation_matrix @ q, q, reaction
