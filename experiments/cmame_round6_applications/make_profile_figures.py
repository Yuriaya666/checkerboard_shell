"""2D centerline profiles for the two application studies."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from experiments.cmame_round6_applications import aircraft_skin, automotive_panel
from experiments.cmame_round6_applications.make_deformation_figures import (
    EXAGGERATION, _positions,
)

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
FIGURES = HERE/"figures"


def _data(name):
    return json.loads((DATA/f"{name}.json").read_text())


def _profile(case_module, name, *, axis, index, xlabel):
    reference, target, optimized, _, _ = _positions(case_module, _data(name))
    if axis == "x":
        line = reference[index, :, :]
        target_line = target[index, :, :]
        optimized_line = optimized[index, :, :]
        coordinate = line[:, 0]
    else:
        line = reference[:, index, :]
        target_line = target[:, index, :]
        optimized_line = optimized[:, index, :]
        coordinate = line[:, 1]
    figure, axes = plt.subplots(1, 2, figsize=(10.6, 3.8))
    axes[0].plot(coordinate, line[:, 2], color="0.25", lw=2.0,
                 label="initial")
    axes[0].plot(coordinate, target_line[:, 2], color="0.65", lw=2.0,
                 ls="--", label="target")
    axes[0].plot(coordinate, optimized_line[:, 2], color="0.05", lw=1.6,
                 label="optimized")
    axes[0].set_xlabel(xlabel); axes[0].set_ylabel("$z$ (×10)")
    axes[0].set_title("shape profile (deformation ×10)")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=8)
    difference = (target_line[:, 2] - optimized_line[:, 2]) / EXAGGERATION
    axes[1].plot(coordinate, 1000.0 * difference, color="0.1", lw=1.8)
    axes[1].set_xlabel(xlabel)
    axes[1].set_ylabel("target − optimized (×1000)")
    axes[1].set_title("profile error")
    axes[1].grid(True, alpha=0.25)
    figure.tight_layout()
    figure.savefig(FIGURES/f"Figure_R6_{name}_profile.png", dpi=220)
    figure.savefig(FIGURES/f"Figure_R6_{name}_profile.pdf")
    plt.close(figure)


def main():
    _profile(aircraft_skin, "aircraft_skin", axis="x", index=4,
             xlabel="$x$")
    _profile(automotive_panel, "automotive_panel", axis="y", index=12,
             xlabel="$y$")


if __name__ == "__main__":
    main()
