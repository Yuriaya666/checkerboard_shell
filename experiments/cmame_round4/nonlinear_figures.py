"""Compact nonlinear verification figure from the Round-4 reruns."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from experiments.cmame_round4.common import DATA, FIGURES
from experiments.cmame_round4.hemisphere import (
    REFERENCE_A_DIGITIZED, REFERENCE_B_DIGITIZED, REFERENCE_LOADS,
)
from experiments.cmame_round4.open_cylinder import (
    REFERENCE_A as CYL_A, REFERENCE_B as CYL_B, REFERENCE_C as CYL_C,
    REFERENCE_LOADS as CYL_LOADS,
)


def _size(path: Path) -> int:
    numbers = [int(part) for part in path.stem.split("_")[-1].split("x")]
    return int(np.prod(numbers))


def _steps(pattern: str) -> list[dict]:
    paths = sorted(Path(DATA).glob(pattern), key=_size)
    if not paths:
        return []
    return json.loads(paths[-1].read_text()).get("steps", [])


def make() -> dict:
    figure, axes = plt.subplots(1, 3, figsize=(13.2, 3.9))

    cantilever = _steps("cantilever/mesh_*x*.json")
    if cantilever:
        axis = axes[0]
        axis.plot([row["tip_x"] for row in cantilever],
                  [row["tip_z"] for row in cantilever],
                  color="0.1", linewidth=1.4, label="Checkerboard")
        axis.plot([row["analytical_tip_x"] for row in cantilever],
                  [row["analytical_tip_z"] for row in cantilever],
                  color="0.55", linestyle="--", linewidth=1.2,
                  label="analytical")
        axis.set_title("cantilever end-moment path")
        axis.set_xlabel("tip $x$")
        axis.set_ylabel("tip $z$")
        axis.legend(fontsize=8)

    hemisphere = _steps("hemisphere/mesh_*.json")
    if hemisphere:
        axis = axes[1]
        axis.plot([row["load_factor"] for row in hemisphere],
                  [row["u_A"] for row in hemisphere],
                  color="0.1", label="Checkerboard $u_A$")
        axis.plot([row["load_factor"] for row in hemisphere],
                  [row["minus_u_B"] for row in hemisphere],
                  color="0.35", linestyle="-.", label="Checkerboard $-u_B$")
        axis.plot(REFERENCE_LOADS, REFERENCE_A_DIGITIZED,
                  color="0.55", linestyle=":", label="digitized $u_A$")
        axis.plot(REFERENCE_LOADS, REFERENCE_B_DIGITIZED,
                  color="0.7", linestyle=":", label="digitized $-u_B$")
        axis.set_title("hemisphere alternating forces")
        axis.set_xlabel("$P/P_{\\max}$")
        axis.set_ylabel("displacement")
        axis.legend(fontsize=7)

    cylinder = _steps("open_cylinder/open_cylinder_pullout/mesh_*x*.json")
    if cylinder:
        axis = axes[2]
        axis.plot([row["load_factor"] for row in cylinder],
                  [row["w_A"] for row in cylinder],
                  color="0.1", label="Checkerboard $w_A$")
        axis.plot([row["load_factor"] for row in cylinder],
                  [row["minus_u_B"] for row in cylinder],
                  color="0.35", linestyle="-.", label="Checkerboard $-u_B$")
        axis.plot([row["load_factor"] for row in cylinder],
                  [row["minus_u_C"] for row in cylinder],
                  color="0.55", linestyle=":", label="Checkerboard $-u_C$")
        axis.plot(CYL_LOADS, CYL_A, color="0.2", linewidth=0.8, alpha=0.5)
        axis.plot(CYL_LOADS, CYL_B, color="0.4", linewidth=0.8, alpha=0.5)
        axis.plot(CYL_LOADS, CYL_C, color="0.6", linewidth=0.8, alpha=0.5)
        axis.set_title("open cylinder pull-out")
        axis.set_xlabel("$P/P_{\\max}$")
        axis.set_ylabel("displacement")
        axis.legend(fontsize=7)

    for axis in axes:
        axis.grid(True, alpha=0.25)
    figure.tight_layout()
    FIGURES.mkdir(parents=True, exist_ok=True)
    figure.savefig(FIGURES/"Figure_R4_nonlinear_paths.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R4_nonlinear_paths.pdf")
    plt.close(figure)
    return {"cantilever_steps": len(cantilever),
            "hemisphere_steps": len(hemisphere),
            "cylinder_steps": len(cylinder)}


if __name__ == "__main__":
    print(json.dumps(make(), indent=2))
