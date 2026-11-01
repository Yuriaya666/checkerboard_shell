"""Exact local-sparse plus rank-one aircraft linear operator.

Keeps the validated shell tangent as its sparse local block and one global
vector.  A bordered SuperLU solve handles the rank-one term without creating
the dense outer product used by the historical screening implementation.
"""
from __future__ import annotations

import time
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch

from checkerboard_shell import BoundaryConditions, CheckerboardShell
from checkerboard_shell.mechanics.low_rank import LowRankTangent
from experiments.engineering_inverse.common import positional_beam_stiffness


def sparse_base_operator(case, hessian_assembly_strategy="cached"):
    started=time.perf_counter()
    mesh=case["mesh"]
    model=CheckerboardShell(mesh,case["material"],
        BoundaryConditions(torch.zeros(mesh.shape,dtype=torch.bool)),case["load"],
        hessian_assembly_strategy=hessian_assembly_strategy)
    tangent,residual,audit=model.low_rank_tangent(
        torch.zeros(3*mesh.num_vertices,dtype=torch.float64))
    shell_tangent_seconds=time.perf_counter()-started
    beam_start=time.perf_counter()
    beam=sum((positional_beam_stiffness(mesh.reference_vertices.numpy(),*line)
              for line in case["beam_lines"]),
             np.zeros(tangent.shape,dtype=float))
    local=(tangent.local+sp.csr_matrix(beam)).tocsr()
    del beam
    build_seconds=time.perf_counter()-beam_start
    fixed=np.array(sorted(set(case["spar"])),dtype=int)
    free=np.setdiff1d(np.arange(3*mesh.num_vertices),fixed)
    reduced=LowRankTangent(local[free][:,free].tocsr(),
                           tangent.vector[free],tangent.coefficient)
    augmented=reduced.augmented_matrix()
    factor_start=time.perf_counter()
    lu=spla.splu(augmented)
    factor_seconds=time.perf_counter()-factor_start
    def solve(rhs):
        rhs=np.asarray(rhs,dtype=float)
        vector=rhs.ndim==1
        if vector:rhs=rhs[:,None]
        extended=np.vstack((rhs,np.zeros((1,rhs.shape[1]))))
        result=lu.solve(extended)[:len(free)]
        return result[:,0] if vector else result
    q0=np.zeros(tangent.shape[0])
    q0[free]=solve(-residual[free])
    external=case["load"].values.numpy().reshape(-1)
    reaction0=(local[fixed,:]@q0 + tangent.coefficient*tangent.vector[fixed]*
               (tangent.vector@q0)-external[fixed])
    def cross_action(influence):
        return (local[fixed,:][:,free]@influence+
                tangent.coefficient*tangent.vector[fixed,None]*
                (tangent.vector[free]@influence)[None,:])
    return dict(q0=q0,free=free,fixed=fixed,reaction0=reaction0,
        residual=residual,factorized=solve,cross_action=cross_action,
        total_dofs=len(q0),free_dofs=len(free),
        shell_local_nnz=int(tangent.local.nnz),
        full_local_nnz=int(local.nnz),
        augmented_nnz=int(augmented.nnz),
        rank_one_storage=tangent.rank_one_storage,
        shell_tangent_seconds=shell_tangent_seconds,
        shell_tangent_audit=audit,
        beam_build_seconds=build_seconds,
        factor_seconds=factor_seconds,
        assembly_seconds=time.perf_counter()-started)
