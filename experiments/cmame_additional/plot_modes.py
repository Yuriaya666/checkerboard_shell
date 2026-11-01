#!/usr/bin/env python3
"""Replot stored eigenvectors without recomputing any eigensystem."""

from __future__ import annotations

from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parents[1]/"src"))
sys.path.insert(0, str(HERE))

import numpy as np
from common import ROOT, plot_style, save_figure


def draw_cases(experiment, selected, output):
    plot_style()
    import matplotlib.pyplot as plt

    fig = plt.figure(figsize=(7.1, 3.2))
    for panel, (tag, title) in enumerate(selected, 1):
        path = ROOT/experiment/"raw_data"/"mode_shapes"/f"modes_{tag}.npz"
        data = np.load(path)
        reference = data["reference"]
        mode = data["modes"][0]
        scale = 0.28*np.linalg.norm(np.ptp(reference.reshape(-1, 3), axis=0))
        scale /= max(float(np.linalg.norm(mode, axis=-1).max()), 1e-30)
        deformed = reference+scale*mode
        ax = fig.add_subplot(1, len(selected), panel, projection="3d")
        for field, color, alpha in ((reference, "0.55", 0.45),
                                    (deformed, "C0", 0.9)):
            for line in field:
                ax.plot(*line.T, c=color, alpha=alpha, lw=0.55)
            for line in field.transpose(1, 0, 2):
                ax.plot(*line.T, c=color, alpha=alpha, lw=0.55)
        ax.set_title(f"{title}\n$\\lambda_1$={data['eigenvalues'][0]:.2e}")
        ax.set_xlabel("x"); ax.set_ylabel("y"); ax.set_zlabel("z")
        ax.tick_params(labelsize=8, pad=0)
    fig.text(0.5, 0.01, "Gray: reference; blue: mode displaced (independently scaled for visibility)",
             ha="center", fontsize=8)
    save_figure(fig, ROOT/experiment/"figures"/output)
    plt.close(fig)


def main():
    draw_cases("05_curved_spectrum", [
        ("sphere_16_full_tangent", "positive K, n=16"),
        ("saddle_16_full_tangent", "negative K, n=16"),
    ], "Figure_K3")
    draw_cases("06_boundary_modes", [
        ("flat_16_full_tangent", "free flat strip, n=16"),
        ("cylinder_16_full_tangent", "free cylindrical strip, n=16"),
    ], "Figure_E3")


if __name__ == "__main__":
    main()
