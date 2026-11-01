"""Shared machinery for the aircraft-skin and automotive-panel studies.

Both cases are force-inverse problems on a thin curved panel.  The Checkerboard
model uses the production 3V tangent; the MITC4 baseline uses OpenSees
ShellMITC4 with central finite-difference design gradients.
"""

from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import minimize

from checkerboard_shell import (
    BoundaryConditions, IsotropicMaterial, NodalLoad, NewtonSettings,
    CheckerboardShell,
)
from experiments.cmame_round5.common import (
    cantilever_mesh, force_matrix, linear_solve, node_id, selection_matrix,
    tangent_at_reference,
)

ALPHA = 1.0e-10
OPTIONS = {"maxiter": 30, "ftol": 1e-14, "gtol": 1e-10, "maxls": 60}


def panel_mesh(nx, ny, width, height, heights):
    return cantilever_mesh(nx, ny, width, height, heights=heights)


def gaussian_weights(mesh, nodes, *, x0, y0, sx, sy, total, x_range):
    vertices = mesh.reference_vertices.detach().cpu().numpy()
    x = np.array([vertices[node // (mesh.nx + 1), node % (mesh.nx + 1), 0]
                  for node in nodes])
    y = np.array([vertices[node // (mesh.nx + 1), node % (mesh.nx + 1), 1]
                  for node in nodes])
    weights = np.exp(-((x - x0) / sx) ** 2 - ((y - y0) / sy) ** 2)
    weights *= (x >= x_range[0]) & (x <= x_range[1])
    if weights.sum() <= 0:
        raise ValueError("empty Gaussian patch")
    return total * weights / weights.sum()


def checkerboard_problem(*, nx, ny, width, height, heights, material,
                         fixed_nodes, force_nodes, obs_nodes, f_true):
    started = time.perf_counter()
    mesh = panel_mesh(nx, ny, width, height, heights)
    boundary = BoundaryConditions.fixed_nodes(mesh, fixed_nodes)
    model = CheckerboardShell(
        mesh, material, boundary, NodalLoad.zero(mesh),
        NewtonSettings(residual_tolerance=1.0e-12,
                       maximum_newton_iterations=8))
    tangent, residual, _ = tangent_at_reference(model)
    P = force_matrix(model, force_nodes)
    S = selection_matrix(model, obs_nodes)
    q_true, _ = linear_solve(tangent, P @ f_true)
    target = S @ q_true
    build = time.perf_counter() - started
    return {
        "mesh": mesh, "model": model, "tangent": tangent, "P": P, "S": S,
        "f_true": f_true, "target": target, "build_seconds": build,
    }


def checkerboard_design(problem, *, max_iterations=30):
    counter = {"solves": 0}

    def solve(rhs):
        counter["solves"] += 1
        return linear_solve(problem["tangent"], rhs)[0]

    def objective(forces):
        q = solve(problem["P"] @ forces)
        residual = problem["S"] @ q - problem["target"]
        return float(residual @ residual + ALPHA * forces @ forces), q, residual

    def gradient(forces):
        q = solve(problem["P"] @ forces)
        residual = problem["S"] @ q - problem["target"]
        adjoint = solve(2.0 * problem["S"].T @ residual)
        return problem["P"].T @ adjoint + 2.0 * ALPHA * forces

    def function(values):
        value, _, _ = objective(values)
        return value, gradient(values)

    started = time.perf_counter()
    result = minimize(
        function, np.zeros_like(problem["f_true"]), jac=True,
        method="L-BFGS-B", options=OPTIONS)
    runtime = time.perf_counter() - started
    value, q, residual = objective(result.x)
    target_error = float(np.linalg.norm(residual)
                         / max(np.linalg.norm(problem["target"]), 1e-30))
    force_error = float(np.linalg.norm(result.x - problem["f_true"])
                        / max(np.linalg.norm(problem["f_true"]), 1e-30))
    return {
        "method": "Checkerboard + adjoint",
        "f_opt": result.x.tolist(),
        "optimized": (problem["S"] @ q).tolist(),
        "iterations": int(result.nit),
        "success": bool(result.success),
        "message": str(result.message),
        "objective_initial": float(objective(np.zeros_like(result.x))[0]),
        "objective_optimal": float(value),
        "target_error": target_error,
        "force_error": force_error,
        "solves": int(counter["solves"]),
        "optimization_seconds": float(runtime),
        "total_seconds": float(problem["build_seconds"] + runtime),
    }


def _ops():
    from experiments.cmame_round4.mitc4 import _ops as ops_factory
    return ops_factory()


def mitc4_response(*, nx, ny, width, height, heights, thickness, young,
                  poisson, fixed_nodes, force_nodes, forces, obs_nodes):
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for j in range(ny + 1):
        for i in range(nx + 1):
            ops.node(j * (nx + 1) + i + 1, i * width / nx, j * height / ny,
                     float(heights[j, i]))
    ops.section("ElasticMembranePlateSection", 1, young, poisson,
                thickness, 0.0)
    element = 1
    for j in range(ny):
        for i in range(nx):
            ops.element("ShellMITC4", element,
                        j * (nx + 1) + i + 1,
                        j * (nx + 1) + i + 2,
                        (j + 1) * (nx + 1) + i + 2,
                        (j + 1) * (nx + 1) + i + 1, 1)
            element += 1
    for node in fixed_nodes:
        ops.fix(node + 1, 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for node, force in zip(force_nodes, forces):
        ops.load(node + 1, 0.0, 0.0, float(force), 0.0, 0.0, 0.0)
    ops.system("SparseGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    status = ops.analyze(1)
    if status != 0:
        raise RuntimeError(f"MITC4 analysis failed: {status}")
    return np.array([ops.nodeDisp(node + 1, 3) for node in obs_nodes])


def mitc4_design(*, nx, ny, width, height, heights, thickness, young, poisson,
                 fixed_nodes, force_nodes, obs_nodes, target, f_true,
                 max_iterations=30):
    analyses = {"count": 0}
    runtime = {"seconds": 0.0}

    def response(forces):
        started = time.perf_counter()
        values = mitc4_response(
            nx=nx, ny=ny, width=width, height=height, heights=heights,
            thickness=thickness, young=young, poisson=poisson,
            fixed_nodes=fixed_nodes, force_nodes=force_nodes,
            forces=forces, obs_nodes=obs_nodes)
        runtime["seconds"] += time.perf_counter() - started
        analyses["count"] += 1
        return values

    def objective(forces):
        residual = response(forces) - target
        return float(residual @ residual + ALPHA * forces @ forces)

    def gradient(forces, epsilon=1.0e-6):
        result = np.zeros_like(forces)
        for index in range(forces.size):
            direction = np.zeros_like(forces)
            direction[index] = epsilon
            result[index] = (objective(forces + direction)
                             - objective(forces - direction)) / (2.0 * epsilon)
        return result

    initial = objective(np.zeros_like(f_true))
    started = time.perf_counter()
    result = minimize(
        lambda values: (objective(values), gradient(values)),
        np.zeros_like(f_true), jac=True, method="L-BFGS-B",
        options=OPTIONS)
    wall = time.perf_counter() - started
    response_opt = response(result.x)
    residual = response_opt - target
    target_error = float(np.linalg.norm(residual)
                         / max(np.linalg.norm(target), 1e-30))
    force_error = float(np.linalg.norm(result.x - f_true)
                        / max(np.linalg.norm(f_true), 1e-30))
    return {
        "method": "MITC4 + finite difference",
        "f_opt": result.x.tolist(),
        "optimized": response_opt.tolist(),
        "iterations": int(result.nit),
        "success": bool(result.success),
        "message": str(result.message),
        "objective_initial": initial,
        "objective_optimal": float(result.fun),
        "target_error": target_error,
        "force_error": force_error,
        "analyses": int(analyses["count"]),
        "optimization_seconds": float(wall),
        "analysis_seconds": float(runtime["seconds"]),
        "total_seconds": float(wall),
    }


def write_json(path, payload):
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n",
                    encoding="utf-8")
