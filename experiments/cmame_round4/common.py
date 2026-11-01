"""CMAME completion round: shared helpers.

This round does not change the Checkerboard formulation.  It closes the
experimental gaps identified in ``docs/CMAME_PROJECT_SYNTHESIS.md``:

* a full 3V simply supported plate solve (the manuscript driver is the scalar
  transverse branch);
* one adapter, one vertex/load/BC setup, and measured raw/active states for
  the MidedgeTan comparison rows;
* a dual roof reference (historical 0.3024 and modern 0.3005924566);
* reruns of the three finite-deformation benchmarks with the current
  production energy.
"""

from __future__ import annotations

import csv
import hashlib
import json
import os
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np
import torch

from checkerboard_shell import (
    BoundaryConditions, CheckerboardMesh, CheckerboardShell,
    IsotropicMaterial, NodalLoad, NewtonSettings,
)
from checkerboard_shell.mechanics.energy import (
    make_targets, paper_components as production_paper_components,
)


HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
OUT = HERE
DATA = OUT/"raw_data"
FIGURES = OUT/"figures"
LOGS = OUT/"logs"

MIDEDGETAN_BIN = Path(os.environ.get(
    "MIDEDGETAN_BIN", "/tmp/checkerboard_round2_midedgetan_distortion"))

REFERENCES = {
    "plate": 4.4360891e-4,
    "roof_historical": 0.3024,
    "roof_modern": 0.3005924566,
    "M3": 0.09352155,
}


def material(young, poisson, thickness):
    return IsotropicMaterial(young, poisson, thickness)


def energy_parts(vertices, reference, hx, hy, targets,
                 young, poisson, thickness):
    """Current production energy with the historical benchmark signature."""
    return production_paper_components(
        vertices, reference, hx, hy, targets,
        material(young, poisson, thickness))


def relative_error(value, reference):
    return abs(float(value) - float(reference)) / abs(float(reference))


def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n",
                    encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---------------------------------------------------------------------------
# Full 3V simply supported plate
# ---------------------------------------------------------------------------

def full_plate(n: int, *, solve: bool = True) -> dict:
    """Solve the complete 3V simply supported plate on the production energy.

    Boundary conditions: ``w = 0`` on all four edges; in-plane rigid motions
    are removed by three scalar pins, leaving the in-plane field otherwise
    free.  The transverse load is the same interior nodal load used by the
    recovered manuscript scalar branch, so the vertical response is directly
    comparable to that branch.
    """
    mesh = CheckerboardMesh.rectangle(
        width=1.0, height=1.0, nx=n, ny=n)
    reference = mesh.reference_vertices
    mask = torch.zeros_like(reference, dtype=torch.bool)
    mask[0, :, 2] = True
    mask[-1, :, 2] = True
    mask[:, 0, 2] = True
    mask[:, -1, 2] = True
    # Remove the three in-plane rigid modes.  The checkerboard gauge is also
    # annihilated because these pins fix its amplitude to zero.
    mask[0, 0, 0] = True
    mask[0, 0, 1] = True
    mask[n, 0, 1] = True
    boundary = BoundaryConditions(fixed_mask=mask)

    load_values = torch.zeros_like(reference)
    load_values[1:-1, 1:-1, 2] = 1.0e-8 / (n * n)
    load = NodalLoad(values=load_values, name="uniform_pressure_interior")

    model = CheckerboardShell(
        mesh, material(1.0, 0.3, 0.01), boundary, load,
        NewtonSettings(residual_tolerance=1.0e-14,
                       maximum_newton_iterations=5))
    row = {
        "benchmark": "plate_full_3v",
        "n": int(n),
        "raw_scalar_states": int(3 * reference.shape[0] * reference.shape[1]),
        "active_states": int(boundary.free_mask.sum()),
        "total_load_z": float(load_values[..., 2].sum()),
        "status": "NOT_SOLVED",
    }
    if not solve:
        return row
    started = time.perf_counter()
    result = model.solve(backend="sparse_low_rank")
    row["runtime_seconds"] = time.perf_counter() - started
    center = float(
        result.deformed_vertices[n // 2, n // 2, 2]
        - reference[n // 2, n // 2, 2])
    in_plane = (result.deformed_vertices[..., :2]
                - reference[..., :2]).abs().max()
    row.update({
        "status": result.status,
        "iterations": int(result.iterations),
        "residual_inf": float(result.residual_inf),
        "response": center,
        "reference": REFERENCES["plate"],
        "relative_error": relative_error(center, REFERENCES["plate"]),
        "max_in_plane_displacement": float(in_plane),
    })
    return row


# ---------------------------------------------------------------------------
# MidedgeTan adapter (verified official LibShell MidedgeAngleTanFormulation)
# ---------------------------------------------------------------------------

def midedgetan_row(benchmark: str, n: int) -> dict:
    """Run the provenance-recorded MidedgeTan adapter on one uniform grid."""
    from experiments.cmame_round2.run_distortion import (
        logical_grid, param_areas, reference_for,
    )

    if not MIDEDGETAN_BIN.exists():
        raise FileNotFoundError(
            f"MidedgeTan adapter not found at {MIDEDGETAN_BIN}; build it with "
            "experiments/cmame_round2/05_baselines/build_midedgetan_distortion.sh")
    xi, eta = logical_grid(n, 0.0, "smooth", benchmark)
    reference, _, _ = reference_for(benchmark, n, xi, eta)
    weights, signed = param_areas(xi, eta)
    if float(np.min(signed)) <= 0.0:
        raise ValueError("inverted logical cell in uniform grid")
    if benchmark == "roof":
        weights = weights * 25.0 * 2.0 * 25.0 * np.deg2rad(40.0)
    with tempfile.TemporaryDirectory(prefix="cmame_round4_mt_") as directory:
        vertices_path = Path(directory)/"vertices.txt"
        weights_path = Path(directory)/"weights.txt"
        np.savetxt(vertices_path, reference.numpy().reshape(-1, 3), fmt="%.17g")
        np.savetxt(weights_path, weights.reshape(-1), fmt="%.17g")
        completed = subprocess.run(
            [str(MIDEDGETAN_BIN), benchmark, str(n),
             str(vertices_path), str(weights_path)],
            text=True, capture_output=True)
    if completed.returncode != 0:
        return {
            "benchmark": benchmark, "n": int(n), "method": "MidedgeTan",
            "status": "FAILED", "reason": completed.stderr[-2000:],
        }
    payload = json.loads(completed.stdout.strip().splitlines()[-1])
    payload.update({
        "method": "MidedgeTan",
        "status": "SUCCESS",
        "binary_sha256": sha256(MIDEDGETAN_BIN),
    })
    return payload


# ---------------------------------------------------------------------------
# Checkerboard curved benchmarks (recovered manuscript drivers)
# ---------------------------------------------------------------------------

def checkerboard_curved_row(benchmark: str, n: int) -> dict:
    from experiments.cmame_round2.recovered_drivers import (
        solve_m3, solve_roof,
    )
    if benchmark == "roof":
        payload = solve_roof(n)
        payload["method"] = "Checkerboard"
        payload["response"] = float(payload["displacement"])
        payload["reference_historical"] = REFERENCES["roof_historical"]
        payload["reference_modern"] = REFERENCES["roof_modern"]
        payload["relative_error_historical"] = relative_error(
            payload["response"], REFERENCES["roof_historical"])
        payload["relative_error_modern"] = relative_error(
            payload["response"], REFERENCES["roof_modern"])
        return payload
    if benchmark == "M3":
        payload = solve_m3(n)
        payload["method"] = "Checkerboard"
        payload["response"] = float(payload["displacement"])
        payload["reference"] = REFERENCES["M3"]
        payload["relative_error"] = relative_error(
            payload["response"], REFERENCES["M3"])
        return payload
    raise ValueError(benchmark)
