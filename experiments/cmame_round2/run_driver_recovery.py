"""Rerun the recovered manuscript sequences and verified LibShell baselines."""

from __future__ import annotations

import csv
import hashlib
import json
from pathlib import Path
import platform
import subprocess
import time

import numpy as np
import scipy
import torch

from experiments.cmame_round2.recovered_drivers import (
    solve_m3,
    solve_plate,
    solve_roof,
)


ROOT = Path(__file__).resolve().parent
RAW = ROOT / "00_driver_recovery" / "raw_data"


def _sha(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def _write_csv(path: Path, rows: list[dict]) -> None:
    keys = []
    for row in rows:
        keys.extend(key for key in row if key not in keys)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)


def _external(command: list[str]) -> dict:
    completed = subprocess.run(command, check=True, text=True, capture_output=True)
    return json.loads(completed.stdout.strip().splitlines()[-1])


def run() -> dict:
    RAW.mkdir(parents=True, exist_ok=True)
    started = time.time()
    rows = []
    for n in (8, 16, 32, 64):
        rows.append(solve_plate(n))
    for n in (12, 16, 24, 32):
        rows.append(solve_roof(n))
    for n in (6, 8, 12, 16):
        rows.append(solve_m3(n))
    _write_csv(RAW / "checkerboard_manuscript_rerun.csv", rows)

    executables = {
        "plate": Path("/tmp/platepinn_checkerboard_benchmark_build/plate_linear_driver"),
        "roof": Path("/tmp/platepinn_checkerboard_paper_build/scordelis_libshell_driver"),
        "m3": Path("/tmp/platepinn_checkerboard_thin_limit_build/libshell_linear_thickness_driver"),
    }
    baselines = []
    if all(path.is_file() for path in executables.values()):
        for n in (8, 16, 32, 64):
            baselines.append(_external([str(executables["plate"]), "MidedgeTan", str(n)]))
        for n in (12, 16, 24, 32):
            baselines.append(_external([str(executables["roof"]), "MidedgeTan", str(n)]))
        for n in (6, 8, 12, 16):
            baselines.append(_external([str(executables["m3"]), str(n), "0.04"]))
        _write_csv(RAW / "midedgetan_verified_rerun.csv", baselines)

    historical = Path("/home/song/projects/CB-shell")
    sources = [
        historical / "experiments/checkerboard_paper_benchmark/run_mechanics.py",
        historical / "experiments/checkerboard_paper_benchmark/formulation.py",
        historical / "experiments/checkerboard_benchmark/mechanics/linear_plate.py",
        historical / "experiments/checkerboard_benchmark/mechanics/hemisphere.py",
        historical / "experiments/intrinsic_m_general_locking/benchmark_crosscheck.py",
    ]
    payload = {
        "generated_at": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "elapsed_seconds": time.time() - started,
        "environment": {
            "python": platform.python_version(), "numpy": np.__version__,
            "scipy": scipy.__version__, "torch": torch.__version__,
        },
        "historical_repository": {
            "path": str(historical),
            "head": subprocess.check_output(
                ["git", "-C", str(historical), "rev-parse", "HEAD"], text=True).strip(),
            "tag": "v1.0-arxiv",
            "clean": not bool(subprocess.check_output(
                ["git", "-C", str(historical), "status", "--porcelain"], text=True).strip()),
            "source_sha256": {str(path.relative_to(historical)): _sha(path) for path in sources},
        },
        "baseline": {
            "implementation": "official LibShell MidedgeAngleTanFormulation",
            "libshell_commit": "6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf",
            "executable_sha256": {key: _sha(path) for key, path in executables.items()
                                  if path.is_file()},
            "rerun_available": bool(baselines),
        },
        "checkerboard_rows": rows,
        "midedgetan_rows": baselines,
    }
    (RAW / "driver_recovery.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return payload


if __name__ == "__main__":
    result = run()
    print(json.dumps({
        "checkerboard_rows": len(result["checkerboard_rows"]),
        "midedgetan_rows": len(result["midedgetan_rows"]),
        "elapsed_seconds": result["elapsed_seconds"],
    }, indent=2))

