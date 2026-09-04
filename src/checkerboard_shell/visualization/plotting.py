"""Geometry-only visualization; mechanics and solvers do not import it."""

from __future__ import annotations

from pathlib import Path

import numpy as np


def _draw_quad_mesh(ax, quads, vertices, *, color, linewidth, alpha, label=None):
    for index, quad in enumerate(quads):
        cycle = np.r_[quad, quad[0]]
        ax.plot(*vertices[cycle].T, color=color, linewidth=linewidth,
                alpha=alpha, label=label if index == 0 else None)


def plot_deformation(result, *, save_path=None, show=False,
                     displacement_scale=1.0):
    import matplotlib.pyplot as plt

    reference = result.reference_vertices.detach().cpu().numpy().reshape(-1, 3)
    displacement = result.displacement.detach().cpu().numpy().reshape(-1, 3)
    displayed = reference + float(displacement_scale) * displacement
    figure = plt.figure(figsize=(8, 5))
    axis = figure.add_subplot(111, projection="3d")
    _draw_quad_mesh(axis, result.mesh.quads, reference, color="0.65",
                    linewidth=0.7, alpha=0.7, label="reference")
    _draw_quad_mesh(axis, result.mesh.quads, displayed, color="C0",
                    linewidth=1.0, alpha=1.0, label="deformed")
    all_points = np.vstack((reference, displayed))
    center = 0.5 * (all_points.min(0) + all_points.max(0))
    radius = max(float(np.ptp(all_points, axis=0).max()) / 2.0, 1e-12)
    axis.set_xlim(center[0] - radius, center[0] + radius)
    axis.set_ylim(center[1] - radius, center[1] + radius)
    axis.set_zlim(center[2] - radius, center[2] + radius)
    axis.set_box_aspect((1, 1, 1))
    axis.set_xlabel("x")
    axis.set_ylabel("y")
    axis.set_zlabel("z")
    axis.legend()
    axis.set_title(f"Checkerboard shell: {result.status}")
    figure.tight_layout()
    if save_path is not None:
        target = Path(save_path)
        target.parent.mkdir(parents=True, exist_ok=True)
        figure.savefig(target, dpi=180)
    if show:
        plt.show()
    return figure

