"""Validate Round-3 deliverables and create a content-hash manifest."""

from __future__ import annotations

import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent


def strict_json(path: Path):
    return json.loads(
        path.read_text(encoding="utf-8"),
        parse_constant=lambda value: (_ for _ in ()).throw(
            ValueError(f"non-finite JSON constant {value} in {path}")))


def run() -> dict:
    beta = strict_json(ROOT/"01_beta_revalidated"/"raw_data"/"beta_revalidated.json")
    comparison = strict_json(
        ROOT/"02_final_comparison"/"raw_data"/"final_comparative_evidence.json")
    for path in ROOT.rglob("*.json"):
        strict_json(path)
    required = [
        ROOT/"CMAME_ROUND3_REPORT.md",
        ROOT/"02_final_comparison"/"MANUSCRIPT_COMPARISON_TABLE.md",
        ROOT/"02_final_comparison"/"MANUSCRIPT_COMPARISON_TABLE.tex",
        ROOT/"03_final_figures"/"FIGURE_GUIDE.md",
        ROOT/"03_final_figures"/"DEPRECATED_FIGURES.md",
        ROOT/"04_manuscript_map"/"CMAME_ROUND3_MANUSCRIPT_MAP.md",
    ]
    missing = [str(path) for path in required if not path.is_file()]
    if missing:
        raise RuntimeError(f"missing Round-3 deliverables: {missing}")
    figures = sorted((ROOT/"03_final_figures").glob("*.pdf"))
    for pdf in figures:
        stem = pdf.with_suffix("")
        if not stem.with_suffix(".png").is_file():
            raise RuntimeError(f"missing PNG preview for {pdf}")
    beta_summary = strict_json(ROOT/"01_beta_revalidated"/"raw_data"/"summary.json")
    reproduction_max = max(
        item["c_beta_1_manuscript_response_relative_difference"]
        for item in beta_summary["summary"].values())
    summary = {
        "beta_rows": len(beta),
        "beta_failures": sum(row["status"] != "SUCCESS" for row in beta),
        "beta_c1_max_manuscript_response_relative_difference": reproduction_max,
        "comparison_rows": len(comparison),
        "figure_sets_pdf_png_svg": len(figures),
        "deprecated_figure_manifest_present": True,
        "m3_redistribution_limitation_retained": (
            "M3" in (ROOT/"CMAME_ROUND3_REPORT.md").read_text()
            and "28.89%" in (ROOT/"CMAME_ROUND3_REPORT.md").read_text()),
    }
    logs = ROOT/"logs"; logs.mkdir(parents=True, exist_ok=True)
    (logs/"completion_summary.json").write_text(
        json.dumps(summary, indent=2)+"\n", encoding="utf-8")
    manifest_path = ROOT/"artifact_manifest.json"
    files = []
    for path in sorted(p for p in ROOT.rglob("*") if p.is_file()
                       and p != manifest_path and "__pycache__" not in p.parts):
        content = path.read_bytes()
        files.append({"path": str(path.relative_to(ROOT)), "bytes": len(content),
                      "sha256": hashlib.sha256(content).hexdigest()})
    manifest_path.write_text(
        json.dumps({"summary": summary, "files": files}, indent=2)+"\n",
        encoding="utf-8")
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
