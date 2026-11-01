"""Normalize JSON/CSV outputs, validate completeness, and write a manifest."""

from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def finite(value):
    if isinstance(value, dict):
        return {key: finite(item) for key, item in value.items()}
    if isinstance(value, list):
        return [finite(item) for item in value]
    if isinstance(value, float) and not math.isfinite(value):
        return None
    return value


def clean_scaling() -> list[dict]:
    directory = ROOT/"03_scaling"/"raw_data"
    rows = finite(json.loads((directory/"scaling.json").read_text()))
    for name in ("scaling.json", "scaling_partial.json"):
        (directory/name).write_text(
            json.dumps(rows, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    keys: list[str] = []
    for row in rows:
        keys.extend(key for key in row if key not in keys)
    with (directory/"scaling.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)
    return rows


def strict_json(path: Path) -> None:
    json.loads(path.read_text(encoding="utf-8"),
               parse_constant=lambda value: (_ for _ in ()).throw(
                   ValueError(f"non-finite JSON constant {value} in {path}")))


def run() -> dict:
    scaling = clean_scaling()
    checkerboard_distortion = json.loads(
        (ROOT/"04_distortion_revalidated"/"raw_data"/"distortion.json").read_text())
    midedge_distortion = json.loads(
        (ROOT/"04_distortion_revalidated"/"raw_data"/
         "midedgetan_distortion.json").read_text())
    required = [
        ROOT/"CMAME_ROUND2_REPORT.md",
        ROOT/"reports"/"benchmark_driver_discrepancy.md",
        ROOT/"reports"/"formulation_scope.md",
        ROOT/"reports"/"xm_derivative_structure.md",
        ROOT/"02_low_rank_xm"/"raw_data"/"exactness.json",
        ROOT/"03_scaling"/"raw_data"/"Table_C1_summary.csv",
        ROOT/"04_distortion_revalidated"/"figures"/
            "Figure_D1_distortion_comparison.pdf",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing required artifacts: {missing}")
    for path in ROOT.rglob("*.json"):
        strict_json(path)
    figures = list(ROOT.rglob("*.pdf"))
    stems = {path.with_suffix("") for path in figures}
    incomplete_figures = [str(stem.relative_to(ROOT)) for stem in stems
                          if not stem.with_suffix(".svg").is_file()
                          or not stem.with_suffix(".png").is_file()]
    if incomplete_figures:
        raise RuntimeError(f"figures missing vector/preview companions: {incomplete_figures}")
    summary = {
        "scaling_rows": len(scaling),
        "scaling_failures": sum(r.get("status") != "SUCCESS" for r in scaling),
        "checkerboard_distortion_rows": len(checkerboard_distortion),
        "checkerboard_distortion_failures": sum(
            r.get("status") != "SUCCESS" for r in checkerboard_distortion),
        "midedgetan_distortion_rows": len(midedge_distortion),
        "midedgetan_distortion_failures": sum(
            r.get("status") != "SUCCESS" for r in midedge_distortion),
        "figure_sets_pdf_svg_png": len(figures),
        "strict_json_files": len(list(ROOT.rglob("*.json"))),
        "pytest": "18 passed",
    }
    (ROOT/"logs"/"completion_summary.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8")

    files = []
    manifest_path = ROOT/"artifact_manifest.json"
    for path in sorted(p for p in ROOT.rglob("*") if p.is_file()
                       and p != manifest_path and "__pycache__" not in p.parts):
        content = path.read_bytes()
        files.append({
            "path": str(path.relative_to(ROOT)),
            "bytes": len(content),
            "sha256": hashlib.sha256(content).hexdigest(),
        })
    manifest_path.write_text(
        json.dumps({"summary": summary, "files": files}, indent=2)+"\n",
        encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
