"""Actual-units aircraft morphing-skin study (NACA 2412, aluminium, pressure).

Geometry: aft 30% of a NACA 2412 chord, 0.6 m chord, 0.3 m span, 2 mm
aluminium skin.  Load: 5 kPa pressure plus nine actuator forces at three rib
stations.  Target: an external 5 mm trailing-edge camber increase.  The
Checkerboard adjoint and an OpenSees ShellMITC4 finite-difference baseline
solve the same bounded force-design problem.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
import torch
from scipy.optimize import minimize

from checkerboard_shell import (
    BoundaryConditions, CheckerboardShell, IsotropicMaterial, NodalLoad,
    NewtonSettings,
)
from experiments.cmame_round5.common import (
    cantilever_mesh, force_matrix, linear_solve, node_id, selection_matrix,
    tangent_at_reference,
)

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
FIGURES = HERE/"figures"
DATA.mkdir(parents=True, exist_ok=True)
FIGURES.mkdir(parents=True, exist_ok=True)

# Real-scale parameters (SI).
CHORD = 0.60
SPAN = 0.30
X0 = 0.70 * CHORD
X1 = CHORD
THICKNESS = 0.002
YOUNG = 73.0e9
POISSON = 0.33
PRESSURE = 5000.0          # Pa
TARGET_DEFORMATION = 0.005  # m at the trailing edge
FORCE_LIMIT = 500.0        # N per actuator

NX, NY = 12, 8
FORCE_IJ = [(i, j) for j in (1, 4, 7) for i in (2, 6, 10)]
OBS_IJ = [(i, j) for j in range(NY + 1) for i in range(2, NX + 1)]
ALPHA = 1.0e-12
BOUNDS = [(-FORCE_LIMIT, FORCE_LIMIT)] * len(FORCE_IJ)
OPTIONS = {"maxiter": 60, "ftol": 1e-14, "gtol": 1e-8, "maxls": 80}


def naca_camber(xc, m=0.02, p=0.4):
    xc = np.asarray(xc, dtype=float)
    result = np.empty_like(xc)
    fore = xc <= p
    result[fore] = m / p**2 * (2.0 * p * xc[fore] - xc[fore]**2)
    aft = ~fore
    result[aft] = m / (1.0 - p)**2 * (
        (1.0 - 2.0 * p) + 2.0 * p * xc[aft] - xc[aft]**2)
    return result


def build_case():
    base = cantilever_mesh(NX, NY, X1 - X0, SPAN)
    vertices = base.reference_vertices.detach().cpu().numpy()
    x = vertices[..., 0] + X0
    y = vertices[..., 1]
    z0 = CHORD * naca_camber(x / CHORD)
    mesh = cantilever_mesh(NX, NY, X1 - X0, SPAN, heights=z0)
    fixed_nodes = [node_id(mesh, i, j) for i in (0, 1)
                   for j in range(NY + 1)]
    boundary = BoundaryConditions.fixed_nodes(mesh, fixed_nodes)

    hx = (X1 - X0) / NX
    hy = SPAN / NY
    area = np.full((NY + 1, NX + 1), hx * hy)
    area[0, :] *= 0.5
    area[-1, :] *= 0.5
    area[:, 0] *= 0.5
    area[:, -1] *= 0.5
    area[0, 0] *= 0.5
    area[0, -1] *= 0.5
    area[-1, 0] *= 0.5
    area[-1, -1] *= 0.5
    pressure = PRESSURE * area
    load_values = np.zeros((NY + 1, NX + 1, 3))
    load_values[..., 2] = pressure
    material = IsotropicMaterial(YOUNG, POISSON, THICKNESS)
    model = CheckerboardShell(
        mesh, material, boundary, NodalLoad(torch.as_tensor(load_values)),
        NewtonSettings(residual_tolerance=1.0e-10,
                       maximum_newton_iterations=8))
    tangent, residual, _ = tangent_at_reference(model)
    force_nodes = [node_id(mesh, i, j) for i, j in FORCE_IJ]
    obs_nodes = [node_id(mesh, i, j) for i, j in OBS_IJ]
    P = force_matrix(model, force_nodes)
    S = selection_matrix(model, obs_nodes)

    target_z = z0 + TARGET_DEFORMATION * ((x - X0) / (X1 - X0))**2
    target = np.array([target_z[j, i] - z0[j, i] for i, j in OBS_IJ])
    return {
        "mesh": mesh, "model": model, "tangent": tangent,
        "residual": residual, "P": P, "S": S, "target": target,
        "force_nodes": force_nodes, "obs_nodes": obs_nodes,
        "pressure": pressure, "z0": z0, "x": x, "y": y,
        "fixed_nodes": fixed_nodes,
    }


def checkerboard_design(problem):
    counter = {"solves": 0}

    def solve(rhs):
        counter["solves"] += 1
        return linear_solve(problem["tangent"], rhs)[0]

    def objective(forces):
        q = solve(problem["P"] @ forces - problem["residual"])
        residual = problem["S"] @ q - problem["target"]
        return float(residual @ residual + ALPHA * forces @ forces), q, residual

    def gradient(forces):
        q = solve(problem["P"] @ forces - problem["residual"])
        residual = problem["S"] @ q - problem["target"]
        adjoint = solve(2.0 * problem["S"].T @ residual)
        return problem["P"].T @ adjoint + 2.0 * ALPHA * forces

    initial = objective(np.zeros(len(problem["force_nodes"])))[0]
    started = time.perf_counter()
    result = minimize(
        lambda values: (objective(values)[0], gradient(values)),
        np.zeros(len(problem["force_nodes"])), jac=True,
        method="L-BFGS-B", bounds=BOUNDS, options=OPTIONS)
    runtime = time.perf_counter() - started
    value, q, residual = objective(result.x)
    error = float(np.linalg.norm(residual)
                  / max(np.linalg.norm(problem["target"]), 1e-30))
    return {
        "method": "Checkerboard + adjoint",
        "iterations": int(result.nit),
        "success": bool(result.success),
        "objective_initial": float(initial),
        "objective_optimal": float(value),
        "target_error": error,
        "target_error_max": float(np.max(np.abs(residual))),
        "force_max": float(np.max(np.abs(result.x))),
        "force_rms": float(np.sqrt(np.mean(result.x**2))),
        "solves": int(counter["solves"]),
        "optimization_seconds": float(runtime),
        "total_seconds": float(runtime),
        "f_opt": result.x.tolist(),
    }


def _ops():
    from experiments.cmame_round4.mitc4 import _ops as factory
    return factory()


def mitc4_response(problem, forces):
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    mesh = problem["mesh"]
    z = problem["z0"]
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.node(j * (NX + 1) + i + 1, i * (X1 - X0) / NX,
                     j * SPAN / NY, float(z[j, i]))
    ops.section("ElasticMembranePlateSection", 1, YOUNG, POISSON,
                THICKNESS, 0.0)
    element = 1
    for j in range(NY):
        for i in range(NX):
            ops.element("ShellMITC4", element,
                        j * (NX + 1) + i + 1,
                        j * (NX + 1) + i + 2,
                        (j + 1) * (NX + 1) + i + 2,
                        (j + 1) * (NX + 1) + i + 1, 1)
            element += 1
    for node in problem["fixed_nodes"]:
        ops.fix(node + 1, 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(NY + 1):
        for i in range(NX + 1):
            ops.load(j * (NX + 1) + i + 1, 0.0, 0.0,
                     float(problem["pressure"][j, i]), 0.0, 0.0, 0.0)
    for node, force in zip(problem["force_nodes"], forces):
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
    return np.array([ops.nodeDisp(node + 1, 3) for node in problem["obs_nodes"]])


def mitc4_design(problem):
    counter = {"analyses": 0}
    runtime = {"seconds": 0.0}

    def response(forces):
        started = time.perf_counter()
        values = mitc4_response(problem, forces)
        runtime["seconds"] += time.perf_counter() - started
        counter["analyses"] += 1
        return values

    def objective(forces):
        residual = response(forces) - problem["target"]
        return float(residual @ residual + ALPHA * forces @ forces)

    def gradient(forces, epsilon=1.0e-3):
        result = np.zeros_like(forces)
        for index in range(forces.size):
            direction = np.zeros_like(forces)
            direction[index] = epsilon
            result[index] = (objective(forces + direction)
                             - objective(forces - direction)) / (2.0 * epsilon)
        return result

    initial = objective(np.zeros(len(problem["force_nodes"])))
    started = time.perf_counter()
    result = minimize(
        lambda values: (objective(values), gradient(values)),
        np.zeros(len(problem["force_nodes"])), jac=True,
        method="L-BFGS-B", bounds=BOUNDS, options=OPTIONS)
    wall = time.perf_counter() - started
    residual = response(result.x) - problem["target"]
    error = float(np.linalg.norm(residual)
                  / max(np.linalg.norm(problem["target"]), 1e-30))
    return {
        "method": "MITC4 + finite difference",
        "iterations": int(result.nit),
        "success": bool(result.success),
        "objective_initial": float(initial),
        "objective_optimal": float(result.fun),
        "target_error": error,
        "target_error_max": float(np.max(np.abs(residual))),
        "force_max": float(np.max(np.abs(result.x))),
        "force_rms": float(np.sqrt(np.mean(result.x**2))),
        "analyses": int(counter["analyses"]),
        "optimization_seconds": float(wall),
        "analysis_seconds": float(runtime["seconds"]),
        "total_seconds": float(wall),
    }


def run():
    problem = build_case()
    cb = checkerboard_design(problem)
    mt = mitc4_design(problem)
    mt["free_dofs"] = 6 * (NX + 1) * (NY + 1) - 6 * len(problem["fixed_nodes"])
    payload = {
        "case": "aircraft_morphing_skin_actual",
        "geometry": "NACA 2412 aft 30% chord skin",
        "chord_m": CHORD, "span_m": SPAN,
        "panel_x0_m": X0, "panel_x1_m": X1,
        "thickness_m": THICKNESS, "young_Pa": YOUNG,
        "poisson": POISSON, "pressure_Pa": PRESSURE,
        "target_trailing_edge_m": TARGET_DEFORMATION,
        "force_limit_N": FORCE_LIMIT,
        "design_variables": len(problem["force_nodes"]),
        "free_dofs": int(problem["model"].boundary_conditions.free_mask.sum()),
        "observations": len(problem["obs_nodes"]),
        "checkerboard": cb, "mitc4": mt,
        "target_definition": "external NACA 2412 camber + quadratic TE increase",
    }
    (DATA/"aircraft_actual.json").write_text(
        json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
