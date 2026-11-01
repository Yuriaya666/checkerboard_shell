"""Fair same-method comparison: adjoint vs finite difference.

Same Checkerboard model, same target, same optimizer settings, same linear
solver.  The tangent/model assembly is timed once and included in both totals;
the optimization work is also reported separately.
"""

from __future__ import annotations

import json
import time

import numpy as np
from scipy.optimize import minimize

from experiments.cmame_round5.actuation_design import (
    ALPHA, _COUNTER, build_problem, gradient as adjoint_gradient,
    objective as adjoint_objective,
)
from experiments.cmame_round5.common import DATA, LOGS, linear_solve, write_json


OPTIONS = {"maxiter": 40, "ftol": 1e-14, "gtol": 1e-10, "maxls": 60}


def _run_adjoint(problem):
    _COUNTER["linear_solves"] = 0
    history = []

    def function(forces):
        value, _, _ = adjoint_objective(problem, forces)
        return value, adjoint_gradient(problem, forces)

    def callback(forces):
        value, _, _ = adjoint_objective(problem, forces)
        history.append({"iteration": len(history), "objective": float(value)})

    started = time.perf_counter()
    result = minimize(function, np.zeros_like(problem["f_true"]), jac=True,
                      method="L-BFGS-B", callback=callback, options=OPTIONS)
    runtime = time.perf_counter() - started
    q, _ = linear_solve(problem["tangent"],
                        problem["P"] @ result.x)
    error = float(np.linalg.norm(problem["S"] @ q - problem["u_target"])
                  / max(np.linalg.norm(problem["u_target"]), 1e-30))
    return result, history, runtime, _COUNTER["linear_solves"], error


def _run_fd(problem, epsilon=1.0e-6):
    counter = {"solves": 0}
    history = []

    def response(forces):
        counter["solves"] += 1
        q, _ = linear_solve(problem["tangent"], problem["P"] @ forces)
        return problem["S"] @ q

    def value(forces):
        r = response(forces) - problem["u_target"]
        return float(r @ r + ALPHA * forces @ forces)

    def gradient(forces):
        g = np.zeros_like(forces)
        for index in range(forces.size):
            direction = np.zeros_like(forces)
            direction[index] = epsilon
            g[index] = (value(forces + direction)
                        - value(forces - direction)) / (2.0 * epsilon)
        return g

    def function(forces):
        return value(forces), gradient(forces)

    def callback(forces):
        history.append({"iteration": len(history),
                        "objective": float(value(forces))})

    started = time.perf_counter()
    result = minimize(function, np.zeros_like(problem["f_true"]), jac=True,
                      method="L-BFGS-B", callback=callback, options=OPTIONS)
    runtime = time.perf_counter() - started
    error = float(np.linalg.norm(response(result.x) - problem["u_target"])
                  / max(np.linalg.norm(problem["u_target"]), 1e-30))
    return result, history, runtime, counter["solves"], error


def run() -> dict:
    started = time.perf_counter()
    problem = build_problem()
    build_time = time.perf_counter() - started
    _COUNTER["linear_solves"] = 0

    adj, adj_hist, adj_time, adj_solves, adj_error = _run_adjoint(problem)
    fd, fd_hist, fd_time, fd_solves, fd_error = _run_fd(problem)

    rows = [
        {
            "method": "adjoint",
            "iterations": int(adj.nit),
            "solves": int(adj_solves),
            "optimization_seconds": float(adj_time),
            "total_seconds": float(build_time + adj_time),
            "target_error": adj_error,
            "objective": float(adj.fun),
            "success": bool(adj.success),
        },
        {
            "method": "finite_difference",
            "iterations": int(fd.nit),
            "solves": int(fd_solves),
            "optimization_seconds": float(fd_time),
            "total_seconds": float(build_time + fd_time),
            "target_error": fd_error,
            "objective": float(fd.fun),
            "success": bool(fd.success),
        },
    ]
    payload = {
        "case": "fair_same_method_adjoint_vs_fd",
        "design_variables": int(problem["f_true"].size),
        "free_dofs": int(problem["P"].shape[0]),
        "observations": int(problem["S"].shape[0]),
        "build_seconds": float(build_time),
        "fd_step": 1.0e-6,
        "optimizer_options": OPTIONS,
        "rows": rows,
    }
    write_json(DATA/"fair_fd_adjoint.json", payload)
    write_json(LOGS/"fair_fd_adjoint_summary.json", payload)

    lines = [
        "# Same-method: adjoint vs finite difference",
        "",
        "Same Checkerboard model, same target, same L-BFGS settings, same linear",
        "solver.  Model/tangent assembly is timed once and included in `total`;",
        "`optimization` is the design loop only.",
        "",
        "| method | iterations | linear solves | optimization (s) | total (s) | target error |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['method']} | {row['iterations']} | {row['solves']} | "
            f"{row['optimization_seconds']:.2f} | {row['total_seconds']:.2f} | "
            f"{100.0*row['target_error']:.4f}% |")
    lines.extend([
        "",
        f"Model/tangent assembly: {build_time:.2f} s.",
        f"Linear-solve ratio: {rows[1]['solves']/rows[0]['solves']:.1f}x; "
        f"optimization-time ratio: {rows[1]['optimization_seconds']/rows[0]['optimization_seconds']:.1f}x.",
    ])
    (DATA.parent/"FAIR_FD_ADJOINT.md").write_text("\n".join(lines) + "\n",
                                                  encoding="utf-8")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
