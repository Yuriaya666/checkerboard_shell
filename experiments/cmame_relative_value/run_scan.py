"""Thickness/mesh study with the existing curved benchmark drivers.

Run with PYTHONPATH=src:.:/tmp/cb-openseespy and a Python containing torch,
scipy, matplotlib and OpenSeesPy. The MidedgeTan driver is built separately.
All three M3 solves use raw equator point loads and a raw equator observable.
Loads are scaled by (t/t0)^3 so that a pure-bending response is constant.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import os
from pathlib import Path
import subprocess
import tempfile
import time

import numpy as np
import torch

from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import paper_components
from experiments.cmame_round2.recovered_drivers import (
    m3_problem, roof_problem, solve_curved_sparse,
)
from experiments.cmame_round2.run_distortion import logical_grid, param_areas, reference_for
from experiments.cmame_round4.mitc4 import run_m3, run_roof


HERE = Path(__file__).resolve().parent
DATA = HERE / "results"
THICKNESS_RATIOS = (1.0, 0.5, 0.25, 0.125, 0.0625)
MESHES = {"roof": (8, 12, 16, 24), "M3": (6, 8, 12, 16)}
T0 = {"roof": 0.25, "M3": 0.04}
BIN = Path(os.environ.get("MIDEDGETAN_THICKNESS_BIN", "/tmp/checkerboard_midedgetan_thickness"))


def save_rows(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    keys = list(dict.fromkeys(key for row in rows for key in row))
    with path.open("w", newline="") as file:
        writer = csv.DictWriter(file, keys)
        writer.writeheader()
        writer.writerows(rows)
    path.with_suffix(".json").write_text(json.dumps(rows, indent=2) + "\n")


def midedge(benchmark: str, n: int, thickness: float) -> dict:
    xi, eta = logical_grid(n, 0.0, "smooth", benchmark)
    reference, _, _ = reference_for(benchmark, n, xi, eta)
    weights, signed = param_areas(xi, eta)
    assert signed.min() > 0
    if benchmark == "roof":
        weights *= 25.0 * 2.0 * 25.0 * math.radians(40.0)
    with tempfile.TemporaryDirectory() as temp:
        vertices_file, weights_file = Path(temp) / "v.txt", Path(temp) / "w.txt"
        np.savetxt(vertices_file, reference.numpy().reshape(-1, 3), fmt="%.17g")
        np.savetxt(weights_file, weights.reshape(-1), fmt="%.17g")
        process = subprocess.run(
            [str(BIN), benchmark, str(n), str(vertices_file), str(weights_file),
             repr(thickness)], capture_output=True, text=True, check=True)
    return json.loads(process.stdout.strip().splitlines()[-1])


def checkerboard(benchmark: str, n: int, thickness: float, scale: float) -> dict:
    problem = roof_problem(n) if benchmark == "roof" else m3_problem(n)
    old = problem["material"]
    problem["material"] = IsotropicMaterial(old.young, old.poisson, thickness)
    if benchmark == "roof":
        problem["force"] = problem["force"] * scale
    else:
        # Same raw equator point load as the MidedgeTan and MITC4 drivers.
        force = np.zeros(problem["reference"].numel())
        force[3 * n * (n + 1)] = scale
        force[3 * (n * (n + 1) + n) + 1] = -scale
        problem["force"] = force
    start = time.perf_counter()
    displacement, audit = solve_curved_sparse(problem)
    elapsed = time.perf_counter() - start
    if benchmark == "roof":
        response = -float(displacement[problem["output_raw"]])
        active = problem["active"].size
    else:
        response = -float(displacement[3 * (n * (n + 1) + n) + 1])
        active = problem["reference"].numel() - problem["constraints"].shape[0]
    amplitude = 1.0e-4
    positions = problem["reference"] + amplitude * torch.tensor(displacement).reshape_as(problem["reference"])
    with torch.no_grad():
        parts = paper_components(positions, problem["reference"], problem["hx"],
                                 problem["hy"], problem["targets"], problem["material"])
    em = float(parts["membrane"]) / amplitude**2
    eb = float(parts["bending"] + parts["bending_compatibility"]) / amplitude**2
    total = em + eb
    return {"response": response, "raw_states": int(problem["reference"].numel()),
            "active_states": int(active), "time_s": elapsed,
            "tangent_nnz": audit.get("local_tangent_nnz", audit.get("structural_nnz")),
            "membrane_energy_quadratic_reference_tangent": em,
            "bending_energy_quadratic_reference_tangent": eb,
            "membrane_to_bending": em / eb if eb else None,
            "xm_fraction": float(parts["membrane_xm"]) / amplitude**2 / total if total else None,
            "xw_fraction": float(parts["bending_compatibility"]) / amplitude**2 / total if total else None}


def one(benchmark: str, n: int, ratio: float, method: str) -> dict:
    thickness = T0[benchmark] * ratio
    scale = ratio ** 3
    started = time.perf_counter()
    if method == "CB-Full":
        payload = checkerboard(benchmark, n, thickness, scale)
    elif method == "MidedgeTan":
        payload = midedge(benchmark, n, thickness)
        payload = {"response": float(payload["observable"]) * scale,
                   "raw_states": payload["raw_scalar_states"],
                   "active_states": payload["active_states"],
                   "tangent_nnz": payload["tangent_nnz"]}
    else:
        payload = (run_roof(n, thickness=thickness, pressure=90.0 * scale)
                   if benchmark == "roof" else run_m3(n, thickness=thickness))
        payload = {"response": float(payload["response"]) *
                   (1.0 if benchmark == "roof" else scale),
                   "raw_states": payload["raw_scalar_states"],
                   "active_states": payload["active_states"]}
    payload.setdefault("time_s", time.perf_counter() - started)
    payload["total_wall_s"] = time.perf_counter() - started
    return {"benchmark": benchmark, "n": n, "t_over_t0": ratio,
            "thickness": thickness, "load_scale": scale,
            "method": method, **payload, "status": "SUCCESS"}


def run(methods: tuple[str, ...], *, references: bool = False) -> None:
    if "MidedgeTan" in methods and not BIN.exists():
        raise FileNotFoundError(f"build {BIN} with build_midedgetan_distortion.sh")
    rows: list[dict] = []
    path = DATA / ("reference_convergence.csv" if references else "thickness_scan.csv")
    if references:
        cases = ((b, n, r, m) for b in MESHES for r in THICKNESS_RATIOS
                 for m in methods for n in (48, 64, 96))
    else:
        cases = ((b, n, r, m) for b, meshes in MESHES.items()
                 for r in THICKNESS_RATIOS for n in meshes for m in methods)
    for b, n, ratio, method in cases:
        start = time.perf_counter()
        try:
            row = one(b, n, ratio, method)
        except Exception as error:
            row = {"benchmark": b, "n": n, "t_over_t0": ratio,
                   "method": method, "status": "FAILED", "reason": repr(error)}
        rows.append(row)
        save_rows(path, rows)
        print(f"{b:4s} t/t0={ratio:g} n={n:3d} {method:11s} "
              f"{row.get('response', float('nan')):.8g} "
              f"{time.perf_counter()-start:.2f}s {row['status']}", flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--methods", nargs="+", default=["CB-Full", "MidedgeTan", "MITC4"])
    parser.add_argument("--references", action="store_true")
    args = parser.parse_args()
    run(tuple(args.methods), references=args.references)
