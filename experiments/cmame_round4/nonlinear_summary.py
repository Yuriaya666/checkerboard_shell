"""Collect nonlinear benchmark summaries into one table."""

from __future__ import annotations

import json
from pathlib import Path

from experiments.cmame_round4.common import DATA, OUT, write_csv, write_json


def _size(path: Path) -> int:
    numbers = [int(part) for part in path.stem.split("_")[-1].split("x")]
    product = 1
    for number in numbers:
        product *= number
    return product


def collect() -> list[dict]:
    rows: list[dict] = []
    for path in sorted((DATA/"cantilever").glob("mesh_*x*.json"), key=_size):
        payload = json.loads(path.read_text())
        rows.append({
            "case": "cantilever",
            "mesh": "x".join(str(value) for value in payload["mesh_cells"]),
            "classification": payload["classification"],
            "steps": len(payload.get("steps", [])),
            "rejected_steps": len(payload.get("rejected_steps", [])),
            "final_error_1": payload.get("path_RMS_normalized_error"),
            "final_error_2": payload.get("path_max_normalized_error"),
            "final_error_3": payload.get("final_tip_error"),
            "metric_1": "RMS path error",
            "metric_2": "max path error",
            "metric_3": "final tip error",
            "maximum_turning_angle": payload.get("maximum_turning_angle"),
        })
    for path in sorted((DATA/"hemisphere").glob("mesh_*.json"), key=_size):
        payload = json.loads(path.read_text())
        rows.append({
            "case": "hemisphere",
            "mesh": str(payload["mesh_cells"]),
            "classification": payload["classification"],
            "steps": len(payload.get("steps", [])),
            "rejected_steps": len(payload.get("rejected_steps", [])),
            "final_error_1": payload.get("final_A_error_vs_digitized"),
            "final_error_2": payload.get("final_B_error"),
            "final_error_3": payload.get("path_RMS_normalized_error"),
            "metric_1": "final A error",
            "metric_2": "final B error",
            "metric_3": "joint RMS path error",
        })
    for path in sorted((DATA/"open_cylinder"/"open_cylinder_pullout").glob("mesh_*x*.json"), key=_size):
        payload = json.loads(path.read_text())
        metrics = payload.get("path_metrics") or {}
        rows.append({
            "case": "open_cylinder",
            "mesh": f'{payload["mesh_cells_axial"]}x{payload["mesh_cells_circumferential"]}',
            "classification": payload["classification"],
            "steps": len(payload.get("steps", [])),
            "rejected_steps": len(payload.get("rejected_steps", [])),
            "final_error_1": (metrics.get("w_A") or {}).get("final_relative_error"),
            "final_error_2": (metrics.get("minus_u_B") or {}).get("final_relative_error"),
            "final_error_3": (metrics.get("minus_u_C") or {}).get("final_relative_error"),
            "metric_1": "final w_A error",
            "metric_2": "final -u_B error",
            "metric_3": "final -u_C error",
            "combined_RMS": (metrics.get("combined") or {}).get(
                "RMS_normalized_path_error"),
        })
    return rows


def markdown(rows: list[dict]) -> str:
    lines = [
        "# Nonlinear benchmark summary (Round 4)",
        "",
        "| case | mesh | classification | steps | metric 1 | metric 2 | metric 3 |",
        "|---|---:|---|---:|---:|---:|---:|",
    ]
    for row in rows:
        values = []
        for key in ("final_error_1", "final_error_2", "final_error_3"):
            value = row.get(key)
            values.append("—" if value is None else f"{100.0*float(value):.3f}%")
        lines.append(
            f"| {row['case']} | {row['mesh']} | {row['classification']} | "
            f"{row['steps']} | {values[0]} | {values[1]} | {values[2]} |")
    return "\n".join(lines) + "\n"


def run() -> dict:
    rows = collect()
    write_csv(DATA/"nonlinear_summary.csv", rows)
    write_json(DATA/"nonlinear_summary.json", rows)
    (OUT/"NONLINEAR_SUMMARY.md").write_text(markdown(rows), encoding="utf-8")
    return {"rows": len(rows)}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
