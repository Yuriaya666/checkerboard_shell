"""Write hardware and numerical software versions next to timing records."""
from __future__ import annotations

import json
import os
import platform
from importlib.metadata import version,PackageNotFoundError
from pathlib import Path

from experiments.engineering_inverse.aircraft_skin import OUT


def run():
    cpu="unknown"
    for line in Path("/proc/cpuinfo").read_text().splitlines():
        if line.startswith("model name"):
            cpu=line.split(":",1)[1].strip();break
    def package_version(name):
        try:return version(name)
        except PackageNotFoundError:return "not available as package metadata"
    data=dict(cpu=cpu,platform=platform.platform(),python=platform.python_version(),
        packages={name:package_version(name) for name in
            ("numpy","scipy","torch","openseespy","openseespylinux","matplotlib")},
        environment={name:os.environ.get(name) for name in
            ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS")})
    try:
        from experiments.cmame_round4.mitc4 import _ops
        data["opensees_core_version"]=_ops().version()
    except Exception as error:
        data["opensees_core_version_error"]=str(error)
    (OUT/"efficiency_environment.json").write_text(json.dumps(data,indent=2)+"\n")
    return data


if __name__=="__main__":print(json.dumps(run(),indent=2))
