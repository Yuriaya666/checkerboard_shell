"""Design demo 1: adjoint gradient check.

Verifies the exact design sensitivity of the production tangent against central
finite differences for the actuation shape-control problem.
"""

from __future__ import annotations

import json
import numpy as np

from experiments.cmame_round5.actuation_design import (
    build_problem, gradient, objective,
)
from experiments.cmame_round5.common import DATA, OUT, write_json


def run() -> dict:
    problem = build_problem()
    f0 = np.zeros_like(problem["f_true"])
    epsilon = 1.0e-7
    rows = []
    for index in (0, 3, 4, 7, 8):
        direction = np.zeros_like(f0)
        direction[index] = epsilon
        plus, _, _ = objective(problem, f0 + direction)
        minus, _, _ = objective(problem, f0 - direction)
        finite = (plus - minus) / (2.0 * epsilon)
        exact = gradient(problem, f0)[index]
        rows.append({
            "design_index": int(index),
            "adjoint": float(exact),
            "finite_difference": float(finite),
            "absolute_error": float(abs(exact - finite)),
            "relative_error": float(abs(exact - finite)
                                    / max(abs(finite), 1e-30)),
        })
    payload = {
        "case": "adjoint_gradient_check",
        "design_variables": int(f0.size),
        "free_dofs": int(problem["model"].boundary_conditions.free_mask.sum()),
        "observation_count": int(problem["S"].shape[0]),
        "epsilon": epsilon,
        "checks": rows,
        "maximum_relative_error": max(row["relative_error"] for row in rows),
        "note": (
            "Adjoint gradient uses the exact sparse-local plus rank-one "
            "production tangent; finite differences use two forward solves "
            "per design variable."),
    }
    write_json(DATA/"adjoint_check.json", payload)
    lines = [
        "# Adjoint gradient check",
        "",
        f"- design variables: {payload['design_variables']}",
        f"- free DOFs: {payload['free_dofs']}",
        f"- observations: {payload['observation_count']}",
        f"- maximum relative error: {payload['maximum_relative_error']:.3e}",
        "",
        "| index | adjoint | finite difference | relative error |",
        "|---:|---:|---:|---:|",
    ]
    for row in rows:
        lines.append(
            f"| {row['design_index']} | {row['adjoint']:.8e} | "
            f"{row['finite_difference']:.8e} | {row['relative_error']:.3e} |")
    (OUT/"ADJOINT_CHECK.md").write_text("\n".join(lines) + "\n",
                                        encoding="utf-8")
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
