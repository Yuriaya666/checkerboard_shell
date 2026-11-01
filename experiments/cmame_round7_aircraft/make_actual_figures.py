"""3D and profile figures for the actual-units aircraft skin case."""

from __future__ import annotations

import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from experiments.cmame_round7_aircraft import aircraft_actual as case

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
FIGURES = HERE/"figures"
FIGURES.mkdir(parents=True, exist_ok=True)
EXAGGERATION = 10.0


def _full_displacement(model, free_values):
    reference = model.reference_vertices
    values = reference.new_tensor(np.asarray(free_values, dtype=float))
    positions = model.boundary_conditions.positions(reference, values)
    return (positions - reference).detach().cpu().numpy()


def _data():
    return json.loads((DATA/"aircraft_actual.json").read_text())


def _positions(problem, data):
    mesh = problem["mesh"]
    reference = mesh.reference_vertices.detach().cpu().numpy()
    x = problem["x"]
    target_delta = case.TARGET_DEFORMATION * (
        (x - case.X0) / (case.X1 - case.X0))**2
    target = reference.copy()
    target[..., 2] += target_delta
    f_opt = np.asarray(data["checkerboard"]["f_opt"])
    q = case.linear_solve(
        problem["tangent"],
        problem["P"] @ f_opt - problem["residual"])[0]
    optimized = reference + _full_displacement(problem["model"], q)
    return reference, target, optimized


def _surface(axis, positions, *, title, zlim):
    x = positions[..., 0] * 1000.0
    y = positions[..., 1] * 1000.0
    z = positions[..., 2] * 1000.0
    axis.plot_surface(
        x, y, z, cmap="viridis", rstride=1, cstride=1,
        linewidth=0.2, edgecolor="0.35", antialiased=True, shade=True)
    axis.set_xlabel("$x$ [mm]", labelpad=2)
    axis.set_ylabel("$y$ [mm]", labelpad=2)
    axis.set_zlabel("$z$ [mm]", labelpad=2)
    axis.set_title(title, pad=4, fontsize=10)
    axis.set_zlim(*zlim)
    axis.view_init(elev=32, azim=-55)
    axis.set_proj_type("ortho")
    axis.set_box_aspect((1.0, 0.6, 0.7))
    axis.tick_params(labelsize=8)


def deformation_figure():
    problem = case.build_case()
    data = _data()
    reference, target, optimized = _positions(problem, data)
    target = reference + EXAGGERATION * (target - reference)
    optimized = reference + EXAGGERATION * (optimized - reference)
    z_values = np.concatenate([
        reference[..., 2].ravel(), target[..., 2].ravel(),
        optimized[..., 2].ravel()]) * 1000.0
    padding = 0.08 * (z_values.max() - z_values.min() + 1e-9)
    zlim = (z_values.min() - padding, z_values.max() + padding)
    figure = plt.figure(figsize=(13.2, 4.4))
    axes = [figure.add_subplot(1, 3, index, projection="3d")
            for index in range(1, 4)]
    _surface(axes[0], reference, title="initial", zlim=zlim)
    _surface(axes[1], target, title="target (×10)", zlim=zlim)
    _surface(axes[2], optimized, title="optimized (×10)", zlim=zlim)
    figure.suptitle(
        "NACA 2412 morphing skin: shape and deformation (×10; color = $z$)",
        fontsize=12)
    figure.subplots_adjust(left=0.02, right=0.99, bottom=0.04, top=0.84,
                           wspace=0.02)
    figure.savefig(FIGURES/"Figure_R7_aircraft_actual_deformation.png",
                   dpi=220)
    figure.savefig(FIGURES/"Figure_R7_aircraft_actual_deformation.pdf")
    plt.close(figure)


def profile_figure():
    problem = case.build_case()
    data = _data()
    reference, target, optimized = _positions(problem, data)
    j = case.NY // 2
    x = reference[j, :, 0] * 1000.0
    figure, axes = plt.subplots(1, 2, figsize=(10.6, 3.8))
    axes[0].plot(x, reference[j, :, 2] * 1000.0, color="0.25", lw=2.0,
                 label="initial")
    axes[0].plot(x, target[j, :, 2] * 1000.0, color="0.65", lw=2.0,
                 ls="--", label="target")
    axes[0].plot(x, optimized[j, :, 2] * 1000.0, color="0.05", lw=1.6,
                 label="optimized")
    axes[0].set_xlabel("$x$ [mm]"); axes[0].set_ylabel("$z$ [mm]")
    axes[0].set_title("mid-span shape profile")
    axes[0].grid(True, alpha=0.25)
    axes[0].legend(fontsize=8)
    difference = (target[j, :, 2] - optimized[j, :, 2]) * 1000.0
    axes[1].plot(x, difference, color="0.1", lw=1.8)
    axes[1].axhline(0.0, color="0.5", lw=0.8)
    axes[1].set_xlabel("$x$ [mm]")
    axes[1].set_ylabel("target − optimized [mm]")
    axes[1].set_title("profile error")
    axes[1].grid(True, alpha=0.25)
    figure.tight_layout()
    figure.savefig(FIGURES/"Figure_R7_aircraft_actual_profile.png", dpi=220)
    figure.savefig(FIGURES/"Figure_R7_aircraft_actual_profile.pdf")
    plt.close(figure)


if __name__ == "__main__":
    deformation_figure()
    profile_figure()
    print("actual aircraft figures written")
