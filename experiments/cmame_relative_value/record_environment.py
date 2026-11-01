"""Record software and artifact provenance for the measured study."""
from __future__ import annotations

import hashlib
import json
import os
import platform
import subprocess
import sys
from pathlib import Path

import numpy
import scipy
import torch

from .run_scan import BIN, DATA, HERE


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    try:
        import openseespy.opensees as ops
        opensees_version = ops.version()
    except Exception as error:
        opensees_version = f"unavailable: {error}"
    payload = {
        "python": sys.version,
        "platform": platform.platform(),
        "processor": next((line.split(":", 1)[1].strip() for line in
                           open("/proc/cpuinfo") if line.startswith("model name")),
                          platform.processor()),
        "cpu_count_reported": os.cpu_count(),
        "numpy": numpy.__version__, "scipy": scipy.__version__,
        "torch": torch.__version__, "opensees": opensees_version,
        "thread_settings": {name: os.environ.get(name) for name in
                            ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS")},
        "midedgetan_binary_sha256": sha(BIN) if BIN.exists() else None,
        "libshell_source_archive_sha256": sha(Path(os.environ["LIBSHELL_SOURCE_ARCHIVE"]))
            if "LIBSHELL_SOURCE_ARCHIVE" in os.environ else None,
        "source_sha256": {name: sha(HERE/name) for name in
                          ("run_scan.py", "run_completion.py", "analyze.py")},
        "raw_sha256": {name: sha(DATA/name) for name in
                       ("thickness_scan.csv", "reference_convergence.csv",
                        "completion_ablation.csv", "cylinder_isometry.csv")},
    }
    DATA.mkdir(exist_ok=True)
    (DATA/"environment.json").write_text(json.dumps(payload, indent=2)+"\n")


if __name__ == "__main__": main()
