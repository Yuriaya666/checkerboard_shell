"""Unified 12-row Checkerboard / MidedgeTan comparison for CMAME.

Differences from Round 3:

* the plate row is the complete 3V production solve, not the scalar
  transverse branch; the scalar branch is retained as a cross-check;
* every MidedgeTan row is rerun through the same adapter, on the same vertex
  table, load table, and uniform logical grid, and reports measured
  ``raw_scalar_states`` / ``active_states``;
* the roof row reports errors against both the historical 0.3024 and the
  modern 0.3005924566 reference.
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from experiments.cmame_round4.common import (
    DATA, FIGURES, LOGS, MIDEDGETAN_BIN, OUT, REFERENCES,
    checkerboard_curved_row, full_plate, midedgetan_row, relative_error,
    write_csv, write_json,
)


PLATE_MESHES = (8, 16, 32, 64)
ROOF_MESHES = (12, 16, 24, 32)
M3_MESHES = (6, 8, 12, 16)


def build_rows(*, plate_meshes=PLATE_MESHES) -> list[dict]:
    rows: list[dict] = []

    # --- plate: full 3V Checkerboard, scalar cross-check, MidedgeTan ---
    from experiments.cmame_round2.recovered_drivers import solve_plate
    for n in plate_meshes:
        full = full_plate(n)
        scalar = solve_plate(n)
        mt = midedgetan_row("plate", n)
        rows.append({
            "benchmark": "plate", "n": int(n), "method": "Checkerboard-full-3V",
            "response": full.get("response"),
            "reference": REFERENCES["plate"],
            "relative_error": full.get("relative_error"),
            "raw_scalar_states": full["raw_scalar_states"],
            "active_states": full["active_states"],
            "status": full["status"],
            "iterations": full.get("iterations"),
            "residual_inf": full.get("residual_inf"),
            "runtime_seconds": full.get("runtime_seconds"),
            "total_load_z": full.get("total_load_z"),
            "max_in_plane_displacement": full.get("max_in_plane_displacement"),
            "notes": "production paper_components; w=0 boundary; in-plane free except 3 rigid pins",
        })
        rows.append({
            "benchmark": "plate", "n": int(n), "method": "Checkerboard-scalar",
            "response": scalar["displacement"],
            "reference": REFERENCES["plate"],
            "relative_error": scalar["relative_error"],
            "raw_scalar_states": scalar["raw_scalar_states"],
            "active_states": scalar["active_states"],
            "status": "SUCCESS",
            "notes": "recovered manuscript scalar transverse branch (cross-check only)",
        })
        rows.append({
            "benchmark": "plate", "n": int(n), "method": "MidedgeTan",
            "response": mt.get("observable"),
            "reference": REFERENCES["plate"],
            "relative_error": relative_error(mt["observable"], REFERENCES["plate"])
            if mt.get("status") == "SUCCESS" else None,
            "raw_scalar_states": mt.get("raw_scalar_states"),
            "active_states": mt.get("active_states"),
            "tangent_nnz": mt.get("tangent_nnz"),
            "status": mt.get("status"),
            "notes": "official LibShell MidedgeAngleTanFormulation via uniform adapter",
        })

    # --- roof ---
    for n in ROOF_MESHES:
        cb = checkerboard_curved_row("roof", n)
        mt = midedgetan_row("roof", n)
        rows.append({
            "benchmark": "roof", "n": int(n), "method": "Checkerboard",
            "response": cb["response"],
            "reference": REFERENCES["roof_historical"],
            "relative_error": cb["relative_error_historical"],
            "reference_modern": REFERENCES["roof_modern"],
            "relative_error_modern": cb["relative_error_modern"],
            "raw_scalar_states": cb["raw_scalar_states"],
            "active_states": cb["active_states"],
            "status": "SUCCESS",
            "notes": "recovered manuscript driver; dual reference reported",
        })
        rows.append({
            "benchmark": "roof", "n": int(n), "method": "MidedgeTan",
            "response": mt.get("observable"),
            "reference": REFERENCES["roof_historical"],
            "relative_error": relative_error(mt["observable"], REFERENCES["roof_historical"])
            if mt.get("status") == "SUCCESS" else None,
            "reference_modern": REFERENCES["roof_modern"],
            "relative_error_modern": relative_error(mt["observable"], REFERENCES["roof_modern"])
            if mt.get("status") == "SUCCESS" else None,
            "raw_scalar_states": mt.get("raw_scalar_states"),
            "active_states": mt.get("active_states"),
            "tangent_nnz": mt.get("tangent_nnz"),
            "status": mt.get("status"),
            "notes": "official LibShell MidedgeAngleTanFormulation via uniform adapter",
        })

    # --- M3 ---
    for n in M3_MESHES:
        cb = checkerboard_curved_row("M3", n)
        mt = midedgetan_row("M3", n)
        rows.append({
            "benchmark": "M3", "n": int(n), "method": "Checkerboard",
            # Main comparison uses the same raw equator vertex as the
            # MidedgeTan driver; the physical-midpoint reconstruction is kept
            # as a second observable.
            "response": cb["displacement_raw"],
            "reference": REFERENCES["M3"],
            "relative_error": cb["relative_error_raw"],
            "response_midpoint": cb["displacement_midpoint"],
            "relative_error_midpoint": cb["relative_error_midpoint"],
            "raw_scalar_states": cb["raw_scalar_states"],
            "active_states": cb["active_states"],
            "status": "SUCCESS",
            "notes": "main = raw equator vertex; midpoint reconstruction in parentheses",
        })
        rows.append({
            "benchmark": "M3", "n": int(n), "method": "MidedgeTan",
            "response": mt.get("observable"),
            "reference": REFERENCES["M3"],
            "relative_error": relative_error(mt["observable"], REFERENCES["M3"])
            if mt.get("status") == "SUCCESS" else None,
            "raw_scalar_states": mt.get("raw_scalar_states"),
            "active_states": mt.get("active_states"),
            "tangent_nnz": mt.get("tangent_nnz"),
            "status": mt.get("status"),
            "notes": "official LibShell MidedgeAngleTanFormulation via uniform adapter; raw M3 vertex observable",
        })
    return rows


def markdown_table(rows: list[dict]) -> str:
    lines = [
        "# Unified CMAME comparison (Round 4)",
        "",
        "Checkerboard plate rows use the complete 3V production solve. "
        "MidedgeTan rows are rerun through one adapter on the same uniform "
        "logical grids. Roof errors are given against both references.",
        "",
        "| benchmark | n | Checkerboard response (error) | MidedgeTan response (error) | CB raw/active | MT raw/active |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for benchmark in ("plate", "roof", "M3"):
        for n in sorted({row["n"] for row in rows if row["benchmark"] == benchmark}):
            cb = next((row for row in rows if row["benchmark"] == benchmark
                       and row["n"] == n and row["method"].startswith("Checkerboard")
                       and row["method"] != "Checkerboard-scalar"), None)
            mt = next((row for row in rows if row["benchmark"] == benchmark
                       and row["n"] == n and row["method"] == "MidedgeTan"), None)
            if cb is None or mt is None:
                continue
            cb_error = 100.0*float(cb["relative_error"])
            mt_error = 100.0*float(mt["relative_error"])
            if benchmark == "roof":
                cb_error_text = (f"{cb_error:.3f}% / {100.0*float(cb['relative_error_modern']):.3f}%")
                mt_error_text = (f"{mt_error:.3f}% / {100.0*float(mt['relative_error_modern']):.3f}%")
            elif benchmark == "M3":
                cb_error_text = (f"{cb_error:.4f}% / "
                                 f"{100.0*float(cb['relative_error_midpoint']):.4f}%")
                mt_error_text = f"{mt_error:.4f}%"
            else:
                cb_error_text = f"{cb_error:.4f}%"
                mt_error_text = f"{mt_error:.4f}%"
            lines.append(
                f"| {benchmark} | {n} | {cb['response']:.8g} ({cb_error_text}) | "
                f"{mt['response']:.8g} ({mt_error_text}) | "
                f"{cb['raw_scalar_states']}/{cb['active_states']} | "
                f"{mt['raw_scalar_states']}/{mt['active_states']} |")
    if any(row["benchmark"] == "roof" for row in rows):
        lines.extend([
            "",
            "For roof, the two errors are relative to the historical `0.3024` "
            "and the modern `0.3005924566` reference, respectively.",
        ])
    if any(row["benchmark"] == "M3" for row in rows):
        lines.extend([
            "",
            "For M3, the two errors are the raw equator vertex (main, the same "
            "discrete observable as MidedgeTan) and the physical-midpoint "
            "moment reconstruction (in parentheses), respectively.",
        ])
    return "\n".join(lines) + "\n"


def make_figure(rows: list[dict]) -> None:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    FIGURES.mkdir(parents=True, exist_ok=True)
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))
    for axis, benchmark in zip(axes, ("plate", "roof", "M3")):
        checkerboard_method = (
            "Checkerboard-full-3V" if benchmark == "plate" else "Checkerboard")
        for method, marker, style, label in (
            (checkerboard_method, "o", "-", "Checkerboard"),
            ("MidedgeTan", "s", "--", "MidedgeTan"),
        ):
            selected = [row for row in rows if row["benchmark"] == benchmark
                        and row["method"] == method and row.get("relative_error") is not None]
            if not selected:
                continue
            x = [row["raw_scalar_states"] for row in selected]
            y = [100.0*float(row["relative_error"]) for row in selected]
            axis.plot(x, y, marker=marker, linestyle=style, color="0.15" if method == "Checkerboard" else "0.55",
                      label=label, markersize=5)
        if benchmark == "roof":
            selected = [row for row in rows if row["benchmark"] == "roof"
                        and row["method"] == "Checkerboard"
                        and row.get("relative_error_modern") is not None]
            axis.plot([row["raw_scalar_states"] for row in selected],
                      [100.0*float(row["relative_error_modern"]) for row in selected],
                      marker="^", linestyle=":", color="0.35",
                      label="Checkerboard vs modern", markersize=4)
        if benchmark == "M3":
            selected = [row for row in rows if row["benchmark"] == "M3"
                        and row["method"] == "Checkerboard"
                        and row.get("relative_error_midpoint") is not None]
            axis.plot([row["raw_scalar_states"] for row in selected],
                      [100.0*float(row["relative_error_midpoint"]) for row in selected],
                      marker="^", linestyle=":", color="0.35",
                      label="Checkerboard midpoint", markersize=4)
        axis.set_xscale("log")
        axis.set_yscale("log")
        axis.set_title(benchmark)
        axis.set_xlabel("raw scalar states")
        axis.grid(True, which="both", alpha=0.25)
        axis.legend(fontsize=7)
    axes[0].set_ylabel("continuum-reference error [%]")
    figure.tight_layout()
    figure.savefig(FIGURES/"Figure_R4_unified_comparison.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R4_unified_comparison.pdf")
    plt.close(figure)


def run(*, plate_meshes=PLATE_MESHES) -> dict:
    DATA.mkdir(parents=True, exist_ok=True)
    LOGS.mkdir(parents=True, exist_ok=True)
    rows = build_rows(plate_meshes=plate_meshes)
    write_csv(DATA/"unified_comparison.csv", rows)
    write_json(DATA/"unified_comparison.json", rows)
    (OUT/"UNIFIED_COMPARISON_TABLE.md").write_text(
        markdown_table(rows), encoding="utf-8")
    make_figure(rows)
    summary = {
        "rows": len(rows),
        "plate_meshes": list(plate_meshes),
        "roof_meshes": list(ROOF_MESHES),
        "M3_meshes": list(M3_MESHES),
        "midedgetan_binary": str(MIDEDGETAN_BIN),
        "failed_rows": sum(row.get("status") not in ("SUCCESS", "CONVERGED")
                           for row in rows),
    }
    write_json(LOGS/"unified_comparison_summary.json", summary)
    return summary


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--plate-meshes", type=int, nargs="+",
                        default=list(PLATE_MESHES))
    arguments = parser.parse_args()
    print(json.dumps(run(plate_meshes=tuple(arguments.plate_meshes)), indent=2))
