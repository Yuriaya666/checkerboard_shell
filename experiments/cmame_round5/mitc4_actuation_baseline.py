"""Head-to-head design-loop baseline: MITC4 + finite-difference gradients.

Same 9-actuator target-shape problem as ``actuation_design.py``, but the
forward model is OpenSeesPy ``ShellMITC4`` (6 DOF per node) and the design
gradient is central finite difference.  This measures the practical cost of a
black-box rotation-based shell in a design loop, versus the position-only
adjoint in ``actuation_design.py``.
"""

from __future__ import annotations

import json
import time

import numpy as np

from experiments.cmame_round4.mitc4 import _ops
from experiments.cmame_round5.actuation_design import (
    ALPHA, FORCE_NODES, OBS_NODES, build_problem,
)
from experiments.cmame_round5.common import DATA, LOGS, write_json


N = 12


def _tag(i, j):
    return j * (N + 1) + i + 1


def _build_and_solve(forces):
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    for j in range(N + 1):
        for i in range(N + 1):
            ops.node(_tag(i, j), i / N, j / N, 0.0)
    ops.section("ElasticMembranePlateSection", 1, 2.0e5, 0.3, 0.04, 0.0)
    element = 1
    for j in range(N):
        for i in range(N):
            ops.element("ShellMITC4", element, _tag(i, j), _tag(i + 1, j),
                        _tag(i + 1, j + 1), _tag(i, j + 1), 1)
            element += 1
    # Clamped left boundary: two layers, all six DOF (same physical BC as the
    # Checkerboard clamped_boundary).
    for j in range(N + 1):
        for i in (0, 1):
            ops.fix(_tag(i, j), 1, 1, 1, 1, 1, 1)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for value, (i, j) in zip(forces, FORCE_NODES):
        ops.load(_tag(i, j), 0.0, 0.0, float(value), 0.0, 0.0, 0.0)
    ops.system("SparseGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    status = ops.analyze(1)
    if status != 0:
        raise RuntimeError(f"MITC4 analysis failed: {status}")
    return np.array([ops.nodeDisp(_tag(i, j), 3) for i, j in OBS_NODES])


class Problem:
    def __init__(self):
        reference = build_problem()
        self.f_true = reference["f_true"]
        self.target = _build_and_solve(self.f_true)
        self.analyses = 0
        self.runtime = 0.0

    def response(self, forces):
        started = time.perf_counter()
        values = _build_and_solve(forces)
        self.runtime += time.perf_counter() - started
        self.analyses += 1
        return values

    def objective(self, forces):
        residual = self.response(forces) - self.target
        return float(residual @ residual + ALPHA * forces @ forces)

    def gradient(self, forces, epsilon=1.0e-6):
        gradient = np.zeros_like(forces)
        for index in range(forces.size):
            direction = np.zeros_like(forces)
            direction[index] = epsilon
            plus = self.objective(forces + direction)
            minus = self.objective(forces - direction)
            gradient[index] = (plus - minus) / (2.0 * epsilon)
        return gradient


def run() -> dict:
    from scipy.optimize import minimize
    problem = Problem()
    history = []

    def function(values):
        value = problem.objective(values)
        return value, problem.gradient(values)

    def callback(values):
        history.append({
            "iteration": len(history),
            "objective": float(problem.objective(values)),
        })

    started = time.perf_counter()
    result = minimize(
        function, np.zeros_like(problem.f_true), jac=True,
        method="L-BFGS-B", callback=callback,
        options={"maxiter": 13, "ftol": 1e-14, "gtol": 1e-10, "maxls": 40})
    wall = time.perf_counter() - started
    f_opt = result.x
    response = _build_and_solve(f_opt)
    target_error = float(np.linalg.norm(response - problem.target)
                         / max(np.linalg.norm(problem.target), 1e-30))
    force_error = float(np.linalg.norm(f_opt - problem.f_true)
                        / max(np.linalg.norm(problem.f_true), 1e-30))
    payload = {
        "case": "mitc4_actuation_fd_baseline",
        "method": "OpenSees ShellMITC4 + finite-difference gradient",
        "design_variables": int(problem.f_true.size),
        "nodes": (N + 1) ** 2,
        "elements": N * N,
        "state_dofs": 6 * (N + 1) ** 2,
        "free_dofs": 6 * (N + 1) ** 2 - 6 * 2 * (N + 1),
        "observations": len(OBS_NODES),
        "iterations": int(result.nit),
        "success": bool(result.success),
        "message": str(result.message),
        "target_displacement_relative_error": target_error,
        "force_relative_error": force_error,
        "objective_initial": float(history[0]["objective"]),
        "objective_optimal": float(result.fun),
        "analyses": int(problem.analyses),
        "runtime_seconds": float(wall),
        "analysis_runtime_seconds": float(problem.runtime),
        "history": history,
    }
    write_json(DATA/"mitc4_actuation_baseline.json", payload)
    write_json(LOGS/"mitc4_actuation_baseline_summary.json", {
        "target_error": target_error,
        "iterations": payload["iterations"],
        "analyses": payload["analyses"],
        "runtime_seconds": payload["runtime_seconds"],
    })
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
