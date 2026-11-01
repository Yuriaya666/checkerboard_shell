#!/usr/bin/env python3
"""Build derived tables, experiment manifests and scoped run logs from raw data."""

from __future__ import annotations

import csv
from pathlib import Path
import shutil
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]/"src"))
sys.path.insert(0, str(HERE))

from common import ROOT, write_csv, write_json


def read_csv(path):
    with path.open() as handle:
        return list(csv.DictReader(handle))


def main():
    distortion = read_csv(ROOT/"01_mesh_distortion"/"raw_data"/"all_cases.csv")
    scaling = read_csv(ROOT/"03_scaling"/"raw_data"/"scaling.csv")
    rank_one = read_csv(ROOT/"03_scaling"/"raw_data"/"low_rank_equivalence.csv")
    curved = read_csv(ROOT/"05_curved_spectrum"/"raw_data"/"spectra.csv")
    boundary = read_csv(ROOT/"06_boundary_modes"/"raw_data"/"spectra.csv")

    comparison = []
    for name in ("plate", "roof", "m3"):
        for method in ("CB-Full", "MidedgeTan", "MITC4"):
            matches = [row for row in distortion if row.get("study") == "convergence"
                       and row.get("benchmark") == name and row.get("n") == "12"
                       and row.get("method") == method]
            if matches:
                record = matches[0]
                comparison.append({key: record.get(key, "") for key in (
                    "benchmark", "method", "status", "reference_displacement",
                    "observable", "relative_error", "scalar_state_count",
                    "active_state_count", "condition_estimate", "reason")})
    write_csv(ROOT/"summary_comparison_table.csv", comparison)

    table_c = []
    for benchmark in ("plate", "m3"):
        for method in ("CB-0", "CB-M", "CB-Full", "MidedgeTan", "MITC4"):
            record = next((row for row in scaling if row.get("benchmark") == benchmark
                           and row.get("n") == "16" and row.get("method") == method), None)
            rank = next((row for row in rank_one if row.get("benchmark") == benchmark
                         and row.get("n") == "16"), {})
            table_c.append({
                "benchmark": benchmark, "method": method,
                "status": record["status"] if record else "UNAVAILABLE",
                "states": record.get("scalar_state_count", "") if record else "",
                "active_states": record.get("active_state_count", "") if record else "",
                "tangent_nnz": record.get("tangent_nnz", "") if record else "",
                "tangent_density": record.get("tangent_density", "") if record else "",
                "tangent_assembly_median_s": record.get("explicit_tangent_time_s", "") if record else "",
                "solve_s": record.get("factorization_solve_time_s", "") if record else "",
                "peak_process_rss_mb": record.get("peak_process_rss_mb", "") if record else "",
                "relative_error": record.get("relative_error", "") if record else "",
                "lowrank_augmented_nnz": rank.get("augmented_sparse_nnz", "")
                    if method == "CB-Full" else "",
            })
    write_csv(ROOT/"03_scaling"/"raw_data"/"Table_C1.csv", table_c)

    table_k = []
    for kind in ("sphere", "saddle"):
        for n in (8, 16, 32):
            for operator in ("J_B", "J_B+DXM", "J_bW", "J_bW+DXWrel", "full_tangent"):
                subset = [row for row in curved if row.get("curvature") == kind
                          and row.get("n") == str(n) and row.get("operator") == operator
                          and row.get("eigenvalue")]
                if not subset:
                    table_k.append({"curvature": kind, "n": n, "operator": operator,
                                    "status": "FAILED"})
                    continue
                near = sum(abs(float(row["relative_eigenvalue"])) <= 1e-9
                           for row in subset)
                positive = [float(row["eigenvalue"]) for row in subset
                            if float(row["eigenvalue"]) > 0]
                table_k.append({"curvature": kind, "n": n, "operator": operator,
                                "status": subset[0]["solver_status"],
                                "observed_near_nullity_first_20": near,
                                "smallest_observed_positive": min(positive) if positive else "",
                                "smallest_relative_eigenvalue": subset[0]["relative_eigenvalue"],
                                "global_nullity_not_claimed": True})
    write_csv(ROOT/"05_curved_spectrum"/"raw_data"/"Table_K1.csv", table_k)

    locations = {
        "01_mesh_distortion": ("Primary", ["Experiment 1", "distortion"]),
        "02_gauge_invariance": ("Primary", ["Experiment 2", "Gauge", "gauge"]),
        "03_scaling": ("Primary", ["Experiment 3", "Scaling", "rank-one"]),
        "04_beta_sensitivity": ("Primary", ["Experiment 4", "Beta", "beta"]),
        "05_curved_spectrum": ("Spectral", ["curved", "sphere_", "saddle_"]),
        "06_boundary_modes": ("Spectral", ["boundary", "flat_", "cylinder_"]),
    }
    sources = {
        "Primary": (ROOT/"primary_run.log").read_text().splitlines(),
        "Spectral": (ROOT/"spectral_run.log").read_text().splitlines(),
    }
    for name, (source, patterns) in locations.items():
        base = ROOT/name
        captured = [line for line in sources[source]
                    if any(pattern.lower() in line.lower() for pattern in patterns)]
        (base/"logs"/"run.log").write_text("\n".join(captured)+"\n")
        figures = sorted(str(path.relative_to(base)) for path in (base/"figures").iterdir()
                         if path.is_file())
        datasets = sorted(str(path.relative_to(base)) for path in (base/"raw_data").rglob("*")
                          if path.is_file())
        write_json(base/"artifact_manifest.json", {
            "figures": figures, "raw_data": datasets,
            "config": sorted(str(path.relative_to(base)) for path in (base/"config").iterdir()),
            "logs": ["logs/run.log"],
        })
    print("Built summary comparison, C1, K1, scoped logs and manifests.")


if __name__ == "__main__":
    main()

