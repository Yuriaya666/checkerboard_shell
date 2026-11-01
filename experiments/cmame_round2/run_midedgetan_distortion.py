"""Run the verified LibShell MidedgeTan on the Round-2 distorted vertices."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import subprocess
import tempfile

import numpy as np

from experiments.cmame_round2.run_distortion import (
    LEVELS, MESHES, SEED, logical_grid, param_areas, reference_for,
)
from experiments.cmame_round2.recovered_drivers import (
    M3_REFERENCE, PLATE_REFERENCE, ROOF_REFERENCE,
)


ROOT = Path(__file__).resolve().parent
BASELINE = ROOT/"05_baselines"
OUT = ROOT/"04_distortion_revalidated"/"raw_data"
BINARY = Path("/tmp/checkerboard_round2_midedgetan_distortion")
REFERENCES = {"plate": PLATE_REFERENCE, "roof": ROOF_REFERENCE, "M3": M3_REFERENCE}


def build() -> None:
    completed = subprocess.run(
        ["bash", str(BASELINE/"build_midedgetan_distortion.sh"), str(BINARY)],
        check=True, text=True, capture_output=True)
    (BASELINE/"build_midedgetan_distortion.log").write_text(
        completed.stdout+completed.stderr, encoding="utf-8")


def historical_uniform_values() -> dict[str, float]:
    path = ROOT/"00_driver_recovery"/"raw_data"/"midedgetan_verified_rerun.csv"
    rows = list(csv.DictReader(path.open()))
    return {
        "plate": float(next(r for r in rows if r["method"] == "MidedgeTan"
                            and r["mesh_cells"] == "16")["center_deflection"]),
        "roof": float(next(r for r in rows if r["benchmark"] == "Scordelis_Lo_roof"
                           and r["mesh_cells_each_direction"] == "16")[
                               "downward_free_edge_midspan_displacement"]),
        "M3": float(next(r for r in rows if r["method"] == "MidedgeTan_linear"
                         and r["mesh_cells"] == "12")["minus_u_B"]),
    }


def run() -> list[dict]:
    if not BINARY.exists():
        build()
    expected = historical_uniform_values()
    checkerboard = json.loads((OUT/"distortion.json").read_text())
    rows: list[dict] = []
    with tempfile.TemporaryDirectory(prefix="cb_round2_midedge_") as directory:
        temporary = Path(directory)
        for benchmark, n in MESHES.items():
            for kind in ("smooth", "random"):
                for delta in LEVELS:
                    xi, eta = logical_grid(n, delta, kind, benchmark)
                    reference, _, _ = reference_for(benchmark, n, xi, eta)
                    weights, _ = param_areas(xi, eta)
                    if benchmark == "roof":
                        weights *= 25.0*2.0*25.0*np.deg2rad(40.0)
                    vertices_path = temporary/"vertices.txt"
                    weights_path = temporary/"weights.txt"
                    np.savetxt(vertices_path, reference.numpy().reshape(-1, 3), fmt="%.17g")
                    np.savetxt(weights_path, weights.reshape(-1), fmt="%.17g")
                    completed = subprocess.run(
                        [str(BINARY), benchmark, str(n), str(vertices_path),
                         str(weights_path)], text=True, capture_output=True)
                    if completed.returncode == 0:
                        row = json.loads(completed.stdout.strip().splitlines()[-1])
                        row.update({
                            "status": "SUCCESS", "distortion_kind": kind,
                            "delta": delta, "seed": SEED,
                            "reference": REFERENCES[benchmark],
                            "relative_error": abs(row["observable"]-REFERENCES[benchmark])
                                              /REFERENCES[benchmark],
                        })
                        quality = next(r for r in checkerboard
                                       if r["benchmark"] == benchmark
                                       and r["distortion_kind"] == kind
                                       and r["delta"] == delta)
                        for key in ("minimum_scaled_jacobian", "maximum_aspect_ratio",
                                    "minimum_interior_angle_deg", "maximum_skew_angle_deg",
                                    "minimum_logical_cell_area", "maximum_logical_cell_area"):
                            row[key] = quality[key]
                        if delta == 0.0:
                            row["historical_uniform_value"] = expected[benchmark]
                            row["historical_uniform_relative_difference"] = (
                                abs(row["observable"]-expected[benchmark])
                                /abs(expected[benchmark]))
                            # ASCII vertex transport and a separately compiled
                            # executable change the M3 solve by about 2e-10.
                            if row["historical_uniform_relative_difference"] > 1e-8:
                                row["status"] = "FAILED_UNIFORM_REPRODUCTION_GATE"
                    else:
                        row = {
                            "benchmark": benchmark, "n": n, "method": "MidedgeTan",
                            "distortion_kind": kind, "delta": delta, "seed": SEED,
                            "status": "FAILED", "reason": completed.stderr[-2000:],
                        }
                    rows.append(row)
    keys: list[str] = []
    for row in rows:
        keys.extend(key for key in row if key not in keys)
    with (OUT/"midedgetan_distortion.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)
    (OUT/"midedgetan_distortion.json").write_text(
        json.dumps(rows, indent=2)+"\n", encoding="utf-8")
    provenance = {
        "binary": str(BINARY),
        "binary_sha256": hashlib.sha256(BINARY.read_bytes()).hexdigest(),
        "adapter_source": str(BASELINE/"midedgetan_distortion_driver.cpp"),
        "libshell_commit": "6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf",
        "uniform_reproduction_gate": expected,
        "failed_rows": sum(r["status"] != "SUCCESS" for r in rows),
    }
    (BASELINE/"midedgetan_distortion_provenance.json").write_text(
        json.dumps(provenance, indent=2)+"\n", encoding="utf-8")
    return rows


if __name__ == "__main__":
    result = run()
    print(json.dumps({"rows": len(result),
                      "failed": sum(r["status"] != "SUCCESS" for r in result)}, indent=2))
