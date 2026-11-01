"""3D deformation/morphology figures for the two application studies."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from experiments.cmame_round6_applications import aircraft_skin, automotive_panel
from experiments.cmame_round5.common import linear_solve

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
FIGURES = HERE/"figures"
FIGURES.mkdir(parents=True, exist_ok=True)

EXAGGERATION = 10.0


def _load(name):
    return json.loads((DATA/name).read_text())


def _positions(case_module, data):
    base, heights, fixed_nodes, force_nodes, obs_nodes, weights, problem = (
        case_module.build_case())
    mesh = problem["mesh"]
    reference = mesh.reference_vertices.detach().cpu().numpy()
    q_true, _ = linear_solve(problem["tangent"], problem["P"] @ weights)
    f_opt = np.asarray(data["checkerboard"]["f_opt"])
    q_opt, _ = linear_solve(problem["tangent"], problem["P"] @ f_opt)
    # Use the boundary-condition scatter to recover the full displacement.
    target = reference + EXAGGERATION * _full_displacement(
        problem["model"], q_true)
    optimized = reference + EXAGGERATION * _full_displacement(
        problem["model"], q_opt)
    return reference, target, optimized, q_true, q_opt


def _full_displacement(model, free_values):
    reference = model.mesh.reference_vertices
    free_values = reference.new_tensor(np.asarray(free_values, dtype=float))
    positions = model.boundary_conditions.positions(reference, free_values)
    return (positions - reference).detach().cpu().numpy()


def _plot_surface(axis, positions, values, *, title, vmin, vmax, zlim):
    x = positions[..., 0]
    y = positions[..., 1]
    z = positions[..., 2]
    norm = plt.Normalize(vmin=zlim[0], vmax=zlim[1])
    axis.plot_surface(
        x, y, z, cmap="viridis", norm=norm,
        rstride=1, cstride=1, linewidth=0.2, edgecolor="0.35",
        antialiased=True, shade=True)
    axis.set_xlabel("$x$", labelpad=2)
    axis.set_ylabel("$y$", labelpad=2)
    axis.set_zlabel("$z$", labelpad=2)
    axis.set_title(title, pad=4, fontsize=10)
    axis.set_zlim(*zlim)
    axis.view_init(elev=32, azim=-55)
    axis.set_proj_type("ortho")
    axis.set_box_aspect((1.0, 0.6, 0.8))
    axis.tick_params(labelsize=8)
    return norm


def _deformation_figure(case_module, data, name):
    reference, target, optimized, q_true, q_opt = _positions(case_module, data)
    # Actual (unscaled) displacement magnitude used for the color.
    target_u = target[..., 2] - reference[..., 2]
    optimized_u = optimized[..., 2] - reference[..., 2]
    target_u /= EXAGGERATION
    optimized_u /= EXAGGERATION
    vmin = min(target_u.min(), optimized_u.min(), 0.0)
    vmax = max(target_u.max(), optimized_u.max(), 1e-12)
    z_values = np.concatenate([
        reference[..., 2].ravel(), target[..., 2].ravel(),
        optimized[..., 2].ravel()])
    padding = 0.08 * (z_values.max() - z_values.min() + 1e-9)
    zlim = (z_values.min() - padding, z_values.max() + padding)
    figure = plt.figure(figsize=(13.2, 4.4))
    axes = [figure.add_subplot(1, 3, index, projection="3d")
            for index in range(1, 4)]
    norm = _plot_surface(axes[0], reference, np.zeros_like(reference[..., 2]),
                         title="initial", vmin=vmin, vmax=vmax, zlim=zlim)
    _plot_surface(axes[1], target, target_u,
                  title="target (deformation ×10)",
                  vmin=vmin, vmax=vmax, zlim=zlim)
    _plot_surface(axes[2], optimized, optimized_u,
                  title="optimized (deformation ×10)",
                  vmin=vmin, vmax=vmax, zlim=zlim)
    figure.suptitle(
        f"{name}: shape and deformation (×10; color = $z$)",
        fontsize=12)
    figure.subplots_adjust(left=0.02, right=0.99, bottom=0.04, top=0.84,
                           wspace=0.02)
    figure.savefig(FIGURES/f"Figure_R6_{name.lower().replace(' ', '_')}_deformation.png",
                   dpi=220, bbox_inches="tight")
    figure.savefig(FIGURES/f"Figure_R6_{name.lower().replace(' ', '_')}_deformation.pdf",
                   bbox_inches="tight")
    plt.close(figure)


def main():
    _deformation_figure(
        aircraft_skin, _load("aircraft_skin.json"), "Aircraft skin")
    _deformation_figure(
        automotive_panel, _load("automotive_panel.json"), "Automotive panel")
    print("deformation figures written")


if __name__ == "__main__":
    main()
