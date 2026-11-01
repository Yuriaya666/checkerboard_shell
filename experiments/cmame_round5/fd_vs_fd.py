"""Same-target finite-difference comparison: Checkerboard vs MITC4.

Both models use the same 12x8 flat plate, left-edge clamp, nine z-forces, 99
observation nodes, and the Checkerboard-generated target.  Only the
discretization/backend changes; the gradient method (central finite
difference) and optimizer settings are identical.
"""

from __future__ import annotations

import json
import time
from pathlib import Path

import numpy as np
from scipy.optimize import minimize

from experiments.cmame_round5.actuation_design import build_problem
from experiments.cmame_round5.common import linear_solve, write_json
from experiments.cmame_round6_applications.common import mitc4_response

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
NX, NY = 12, 8
WIDTH, HEIGHT = 1.0, 0.6
THICKNESS = 0.04
YOUNG, POISSON = 2.0e5, 0.3
EPSILON = 1.0e-6
OPTIONS = {"maxiter": 40, "ftol": 1e-14, "gtol": 1e-10, "maxls": 60}
ALPHA = 1.0e-10


def checkerboard_fd(problem):
    counter = {"solves": 0}

    def response(forces):
        counter["solves"] += 1
        q, _ = linear_solve(problem["tangent"], problem["P"] @ forces)
        return problem["S"] @ q

    def objective(forces):
        residual = response(forces) - problem["u_target"]
        return float(residual @ residual + ALPHA * forces @ forces)

    def gradient(forces):
        result = np.zeros_like(forces)
        for index in range(forces.size):
            direction = np.zeros_like(forces)
            direction[index] = EPSILON
            result[index] = (objective(forces + direction)
                             - objective(forces - direction)) / (2.0 * EPSILON)
        return result

    started = time.perf_counter()
    result = minimize(
        lambda values: (objective(values), gradient(values)),
        np.zeros_like(problem["f_true"]), jac=True, method="L-BFGS-B",
        options=OPTIONS)
    runtime = time.perf_counter() - started
    residual = response(result.x) - problem["u_target"]
    target_error = float(np.linalg.norm(residual)
                         / max(np.linalg.norm(problem["u_target"]), 1e-30))
    force_error = float(np.linalg.norm(result.x - problem["f_true"])
                        / max(np.linalg.norm(problem["f_true"]), 1e-30))
    return {
        "method": "Checkerboard + finite difference",
        "iterations": int(result.nit),
        "success": bool(result.success),
        "target_error": target_error,
        "target_error_max": float(np.max(np.abs(residual))),
        "force_error": force_error,
        "solves": int(counter["solves"]),
        "optimization_seconds": float(runtime),
    }


def mitc4_fd(problem, fixed_nodes):
    heights = np.zeros((NY + 1, NX + 1))
    counter = {"analyses": 0}
    runtime = {"seconds": 0.0}

    def response(forces):
        started = time.perf_counter()
        values = mitc4_response(
            nx=NX, ny=NY, width=WIDTH, height=HEIGHT, heights=heights,
            thickness=THICKNESS, young=YOUNG, poisson=POISSON,
            fixed_nodes=fixed_nodes, force_nodes=problem["force_nodes"],
            forces=forces, obs_nodes=problem["observation_nodes"])
        runtime["seconds"] += time.perf_counter() - started
        counter["analyses"] += 1
        return values

    def objective(forces):
        residual = response(forces) - problem["u_target"]
        return float(residual @ residual + ALPHA * forces @ forces)

    def gradient(forces):
        result = np.zeros_like(forces)
        for index in range(forces.size):
            direction = np.zeros_like(forces)
            direction[index] = EPSILON
            result[index] = (objective(forces + direction)
                             - objective(forces - direction)) / (2.0 * EPSILON)
        return result

    started = time.perf_counter()
    result = minimize(
        lambda values: (objective(values), gradient(values)),
        np.zeros_like(problem["f_true"]), jac=True, method="L-BFGS-B",
        options=OPTIONS)
    wall = time.perf_counter() - started
    residual = response(result.x) - problem["u_target"]
    target_error = float(np.linalg.norm(residual)
                         / max(np.linalg.norm(problem["u_target"]), 1e-30))
    force_error = float(np.linalg.norm(result.x - problem["f_true"])
                        / max(np.linalg.norm(problem["f_true"]), 1e-30))
    return {
        "method": "MITC4 + finite difference",
        "iterations": int(result.nit),
        "success": bool(result.success),
        "target_error": target_error,
        "target_error_max": float(np.max(np.abs(residual))),
        "force_error": force_error,
        "analyses": int(counter["analyses"]),
        "optimization_seconds": float(wall),
        "analysis_seconds": float(runtime["seconds"]),
    }


def run():
    started = time.perf_counter()
    problem = build_problem()
    build_seconds = time.perf_counter() - started
    fixed_nodes = [j * (NX + 1) + i for j in range(NY + 1) for i in (0, 1)]
    cb = checkerboard_fd(problem)
    mt = mitc4_fd(problem, fixed_nodes)
    cb["total_seconds"] = float(build_seconds + cb["optimization_seconds"])
    cb["build_seconds"] = float(build_seconds)
    cb["free_dofs"] = int(problem["P"].shape[0])
    mt["total_seconds"] = float(mt["optimization_seconds"])
    mt["free_dofs"] = 6 * (NX + 1) * (NY + 1) - 6 * len(fixed_nodes)
    payload = {
        "case": "same_target_fd_vs_fd",
        "geometry": "flat cantilever plate 12x8",
        "thickness": THICKNESS,
        "young": YOUNG,
        "poisson": POISSON,
        "force_nodes": problem["force_nodes"],
        "observation_nodes": problem["observation_nodes"],
        "design_variables": int(problem["f_true"].size),
        "target": "Checkerboard-generated target",
        "fd_step": EPSILON,
        "checkerboard": cb,
        "mitc4": mt,
    }
    write_json(DATA/"fd_vs_fd.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
