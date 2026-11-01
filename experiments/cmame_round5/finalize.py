"""Record environment, summary, and artifact hashes for Round 5."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy
import torch

from experiments.cmame_round5.common import DATA, HERE, LOGS, ROOT, write_json


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _git(*args: str) -> str:
    completed = subprocess.run(
        ("git", *args), cwd=ROOT, text=True, capture_output=True)
    return completed.stdout.strip()


def _load(name: str) -> dict:
    path = DATA/name
    return json.loads(path.read_text()) if path.exists() else {}


def collect() -> dict:
    adjoint = _load("adjoint_check.json")
    actuation = _load("actuation_design.json")
    shape = _load("shape_design.json")
    rotation = _load("large_rotation_design.json")
    mitc4 = _load("mitc4_actuation_baseline.json")
    comparison = _load("design_loop_comparison.json")
    fair = _load("fair_fd_adjoint.json")
    summary = {
        "adjoint_max_relative_error": adjoint.get("maximum_relative_error"),
        "actuation": {
            "displacement_relative_error": actuation.get(
                "displacement_relative_error"),
            "force_relative_error": actuation.get("force_relative_error"),
            "iterations": actuation.get("iterations"),
            "design_variables": actuation.get("design_variables"),
            "free_dofs": actuation.get("free_dofs"),
        },
        "shape": {
            "displacement_relative_error": shape.get(
                "displacement_relative_error"),
            "height_relative_error": shape.get("height_relative_error"),
            "iterations": shape.get("iterations"),
            "design_variables": shape.get("design_variables"),
            "free_dofs": shape.get("free_dofs"),
        },
        "large_rotation": {
            "tip_error": rotation.get("tip_error"),
            "iterations": rotation.get("iterations"),
            "design_variables": rotation.get("design_variables"),
            "free_dofs": rotation.get("free_dofs"),
        },
        "mitc4_design_loop": {
            "target_error": mitc4.get("target_displacement_relative_error"),
            "iterations": mitc4.get("iterations"),
            "analyses": mitc4.get("analyses"),
            "runtime_seconds": mitc4.get("runtime_seconds"),
            "free_dofs": mitc4.get("free_dofs"),
        },
        "design_loop_comparison": comparison,
        "fair_fd_adjoint": fair,
    }
    write_json(DATA/"design_summary.json", summary)
    environment = {
        "generated_at": __import__("datetime").datetime.now().astimezone().isoformat(),
        "workspace": str(ROOT),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_status_short": _git("status", "--short"),
        "python": sys.version.split()[0],
        "torch": torch.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "platform": platform.platform(),
        "commands": [
            "python -m experiments.cmame_round5.adjoint_check",
            "python -m experiments.cmame_round5.actuation_design",
            "python -m experiments.cmame_round5.shape_design",
            "python -m experiments.cmame_round5.large_rotation_design",
            "python -m experiments.cmame_round5.fair_fd_adjoint",
            "python -m experiments.cmame_round5.mitc4_actuation_baseline",
            "python -m experiments.cmame_round5.design_loop_comparison",
        ],
        "summary": summary,
    }
    write_json(LOGS/"environment.json", environment)
    files = []
    for path in sorted(HERE.rglob("*")):
        if not path.is_file() or "__pycache__" in path.parts:
            continue
        if path.name == "artifact_manifest.json":
            continue
        files.append({
            "path": str(path.relative_to(HERE)),
            "bytes": path.stat().st_size,
            "sha256": _sha(path),
        })
    write_json(HERE/"artifact_manifest.json", {
        "summary": {"files": len(files),
                    "git_commit": environment["git_commit"]},
        "files": files,
    })
    return {"summary": summary, "files": len(files)}


if __name__ == "__main__":
    print(json.dumps(collect(), indent=2))
