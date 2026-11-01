"""Exact production tangent at nonlinear final states.

The archived final states (CB-shell tag ``v1.0-arxiv``) are reused only as
*states*.  The tangent is the current production ``K_local + alpha g g^T``
reduced to the benchmark's kinematic subspace and compared with a dense AD
Hessian of the current production energy.
"""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import torch

from checkerboard_shell.mechanics.energy import make_targets
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.cmame_round4.cantilever import (
    CantileverProblem, MAXIMUM_MOMENT,
)
from experiments.cmame_round4.common import DATA, energy_parts, write_json
from experiments.cmame_round4.hemisphere import HemisphereProblem
from experiments.cmame_round4.open_cylinder import Problem


CB = Path("/home/song/projects/CB-shell/outputs/checkerboard_nonlinear_benchmarks")


def _relative(matrix, reference):
    return float(torch.linalg.matrix_norm(matrix-reference)
                 / max(float(torch.linalg.matrix_norm(reference)), 1e-30))


def _low_rank_tangent(problem, positions, transform):
    reference = problem.reference
    displacement = (positions - reference).reshape(-1)
    tangent, residual, audit = assemble_paper_low_rank_tangent(
        displacement.reshape_as(reference), reference,
        problem.hx, problem.hy, problem.targets,
        IsotropicMaterial(*problem.material_args),
        variant="CB-Full")
    T = transform.detach().cpu().numpy()
    local = T.T @ tangent.local @ T
    vector = T.T @ tangent.vector
    reduced = local + tangent.coefficient * np.outer(vector, vector)
    return torch.as_tensor(reduced, dtype=reference.dtype), tangent, audit


def check_open_cylinder() -> dict:
    archived = json.loads((CB/"open_cylinder_pullout/mesh_8x8.json").read_text())
    positions = torch.tensor(archived["steps"][-1]["positions"], dtype=torch.float64)
    problem = Problem(8, 8)
    problem.material_args = (10.5e6, 0.3125, 0.094)
    values = problem.transform.T @ (positions.reshape(-1) - problem.base)
    dense = torch.autograd.functional.hessian(
        lambda x: problem.potential(x, 1.0), values, vectorize=True)
    dense = 0.5*(dense+dense.T)
    reduced, tangent, audit = _low_rank_tangent(
        problem, positions, problem.transform)
    return {
        "case": "open_cylinder_8x8_final",
        "reduced_dimension": int(values.numel()),
        "relative_tangent_difference": _relative(reduced, dense),
        "rank_one_coefficient": float(tangent.coefficient),
        "rank_one_storage": int(tangent.rank_one_storage),
        "X_M": float(audit["X_M"]),
    }


def check_hemisphere() -> dict:
    archived = json.loads((CB/"hemisphere/mesh_12.json").read_text())
    positions = torch.tensor(archived["steps"][-1]["positions"], dtype=torch.float64)
    problem = HemisphereProblem(12)
    problem.material_args = (6.825e7, 0.3, 0.04)
    values = problem.transform.T @ (positions.reshape(-1) - problem.base)
    dense = torch.autograd.functional.hessian(
        lambda x: problem.potential(x, 1.0), values, vectorize=True)
    dense = 0.5*(dense+dense.T)
    reduced, tangent, audit = _low_rank_tangent(
        problem, positions, problem.transform)
    return {
        "case": "hemisphere_12_final",
        "reduced_dimension": int(values.numel()),
        "relative_tangent_difference": _relative(reduced, dense),
        "rank_one_coefficient": float(tangent.coefficient),
        "rank_one_storage": int(tangent.rank_one_storage),
        "X_M": float(audit["X_M"]),
    }


def check_cantilever() -> dict:
    archived = json.loads((CB/"cantilever/mesh_32x4.json").read_text())
    positions = torch.tensor(archived["steps"][-1]["positions"], dtype=torch.float64)
    problem = CantileverProblem(32, 4)
    problem.material_args = (1.2e6, 0.0, 0.1)
    values = torch.stack((
        positions[2, 2:33, 0]-problem.reference[2, 2:33, 0],
        positions[2, 2:33, 2]-problem.reference[2, 2:33, 2],
    ), -1).reshape(-1)
    # Compare the internal energy only; the end-moment work is an external
    # potential and its Hessian is checked separately by the load path.
    dense = torch.autograd.functional.hessian(
        lambda x: problem.parts(x)["total"], values, vectorize=True)
    dense = 0.5*(dense+dense.T)
    # Linear map from the reduced cantilever coordinates to full positions.
    columns = []
    for index in range(values.numel()):
        basis = torch.zeros_like(values)
        basis[index] = 1.0
        columns.append((problem.current(basis)-problem.reference).reshape(-1))
    transform = torch.stack(columns, dim=1)
    reduced, tangent, audit = _low_rank_tangent(
        problem, positions, transform)
    return {
        "case": "cantilever_32x4_final",
        "reduced_dimension": int(values.numel()),
        "relative_tangent_difference": _relative(reduced, dense),
        "rank_one_coefficient": float(tangent.coefficient),
        "rank_one_storage": int(tangent.rank_one_storage),
        "X_M": float(audit["X_M"]),
        "end_moment": float(MAXIMUM_MOMENT),
    }


def run() -> dict:
    checks = [check_open_cylinder(), check_hemisphere(), check_cantilever()]
    payload = {
        "checks": checks,
        "maximum_relative_tangent_difference": max(
            check["relative_tangent_difference"] for check in checks),
        "note": (
            "States are the archived final configurations; the tangent and "
            "energy are the current production implementation."),
    }
    write_json(DATA/"tangent_final_states.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
