"""Record environment, provenance, and artifact hashes for Round 4."""

from __future__ import annotations

import json
import platform
import subprocess
import sys
from pathlib import Path

import numpy as np
import scipy
import torch

from experiments.cmame_round4.common import (
    DATA, HERE, LOGS, MIDEDGETAN_BIN, ROOT, sha256, write_json,
)


LIBSHELL_COMMIT = "6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf"
EIGEN_VERSION = "3.4.0"
LIBSHELL_SOURCE = "https://github.com/evouga/libshell"
EIGEN_SOURCE = "https://gitlab.com/libeigen/eigen/-/archive/3.4.0/eigen-3.4.0.tar.gz"


def _openseespy_version() -> str | None:
    try:
        import openseespy.opensees as ops
    except ImportError:
        import sys
        sys.path.insert(0, "/tmp/openseespy_pkg")
        try:
            import openseespy.opensees as ops
        except ImportError:
            return None
    return str(ops.version())


def _git(*arguments: str) -> str:
    completed = subprocess.run(
        ("git", *arguments), cwd=ROOT, text=True, capture_output=True)
    return completed.stdout.strip()


def collect() -> dict:
    LOGS.mkdir(parents=True, exist_ok=True)
    environment = {
        "generated_at": __import__("datetime").datetime.now().astimezone().isoformat(),
        "workspace": str(ROOT),
        "git_commit": _git("rev-parse", "HEAD"),
        "git_branch": _git("rev-parse", "--abbrev-ref", "HEAD"),
        "git_status_short": _git("status", "--short"),
        "python": sys.version.split()[0],
        "executable": sys.executable,
        "platform": platform.platform(),
        "torch": torch.__version__,
        "numpy": np.__version__,
        "scipy": scipy.__version__,
        "midedgetan_binary": str(MIDEDGETAN_BIN),
        "midedgetan_binary_sha256": (
            sha256(MIDEDGETAN_BIN) if MIDEDGETAN_BIN.exists() else None),
        "openseespy": _openseespy_version(),
        "openseespy_path": "/tmp/openseespy_pkg",
        "mitc4_element": "OpenSees ShellMITC4",
        "mitc4_section": "ElasticMembranePlateSection",
        "libshell_source": LIBSHELL_SOURCE,
        "libshell_commit": LIBSHELL_COMMIT,
        "eigen_source": EIGEN_SOURCE,
        "eigen_version": EIGEN_VERSION,
        "commands": [
            "python -m experiments.cmame_round4.unified_comparison",
            "python -m experiments.cmame_round4.cantilever --meshes 16 24 32",
            "python -m experiments.cmame_round4.hemisphere --meshes 6 8 12 16",
            "python -m experiments.cmame_round4.open_cylinder --meshes 8x8 12x12 16x16",
            "python -m experiments.cmame_round4.plate_full",
            "python -m experiments.cmame_round4.nonlinear_summary",
            "python -m experiments.cmame_round4.nonlinear_figures",
            "python -m experiments.cmame_round4.tangent_final_states",
            "python -m experiments.cmame_round4.mitc4",
            "python -m experiments.cmame_round4.combined_comparison",
        ],
    }
    write_json(LOGS/"environment.json", environment)

    files = []
    for path in sorted(HERE.rglob("*")):
        if not path.is_file():
            continue
        if path.name in {"artifact_manifest.json"}:
            continue
        if "__pycache__" in path.parts:
            continue
        files.append({
            "path": str(path.relative_to(HERE)),
            "bytes": path.stat().st_size,
            "sha256": sha256(path),
        })
    manifest = {"summary": {
        "files": len(files),
        "generated_at": environment["generated_at"],
        "git_commit": environment["git_commit"],
    }, "files": files}
    write_json(HERE/"artifact_manifest.json", manifest)
    return {"environment": environment, "manifest": {"files": len(files)}}


if __name__ == "__main__":
    print(json.dumps(collect(), indent=2))
