"""Run isolated production-scaling workers and preserve every result."""

from __future__ import annotations

import csv
import json
import os
from pathlib import Path
import subprocess
import sys


ROOT = Path(__file__).resolve().parent
OUT = ROOT/"03_scaling"/"raw_data"
LOGS = ROOT/"03_scaling"/"logs"


def run() -> list[dict]:
    OUT.mkdir(parents=True, exist_ok=True); LOGS.mkdir(parents=True, exist_ok=True)
    rows = []
    for benchmark in ("plate", "M3"):
        for n in (16, 32, 64, 128):
            for method in ("CB-0", "CB-M", "CB-Full", "MidedgeTan"):
                command = [sys.executable, "-m", "experiments.cmame_round2.scaling_worker",
                           benchmark, str(n), method]
                environment = os.environ.copy()
                environment.update({
                    "CB_SCALING_THREADS": "8",
                    "OMP_NUM_THREADS": "8",
                    "MKL_NUM_THREADS": "8",
                    "OPENBLAS_NUM_THREADS": "8",
                })
                completed = subprocess.run(
                    command, text=True, capture_output=True, env=environment)
                (LOGS/f"{benchmark}_{n}_{method}.log").write_text(
                    completed.stdout+completed.stderr, encoding="utf-8")
                if completed.returncode == 0:
                    row = json.loads(completed.stdout.strip().splitlines()[-1])
                else:
                    row = {"benchmark": benchmark, "n": n, "method": method,
                           "status": "FAILED", "reason": completed.stderr[-4000:]}
                rows.append(row)
                (OUT/"scaling_partial.json").write_text(
                    json.dumps(rows, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    keys = []
    for row in rows:
        keys.extend(key for key in row if key not in keys)
    with (OUT/"scaling.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys)
        writer.writeheader(); writer.writerows(rows)
    (OUT/"scaling.json").write_text(
        json.dumps(rows, indent=2, allow_nan=False)+"\n", encoding="utf-8")
    (ROOT/"03_scaling"/"config.json").write_text(json.dumps({
        "mesh_sizes": [16, 32, 64, 128],
        "benchmarks": ["plate", "M3"],
        "methods": ["CB-0", "CB-M", "CB-Full", "MidedgeTan"],
        "formal_repeats": 3,
        "excluded_warmups": 1,
        "torch_threads": 8,
        "worker_isolation": "one fresh process per benchmark/n/method row",
        "timing_note": "MidedgeTan uses a different C++/LibShell backend; its CPU time is auxiliary only.",
    }, indent=2)+"\n", encoding="utf-8")
    return rows


if __name__ == "__main__":
    results = run()
    print(json.dumps({"rows": len(results),
                      "failed": sum(r.get("status") != "SUCCESS" for r in results)},
                     indent=2))
