"""Shared utilities for the four design demos (Round 5).

The demos use the frozen Checkerboard production energy.  Design variables are
either nodal forces (actuation) or reference heights (shape design); the state
is always the 3V free displacement of the production model.
"""

from __future__ import annotations

import csv
import json
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

import numpy as np
import torch

from checkerboard_shell import (
    BoundaryConditions, CheckerboardMesh, CheckerboardShell,
    IsotropicMaterial, NodalLoad, NewtonSettings,
)
from checkerboard_shell.mechanics.energy import make_targets


HERE = Path(__file__).resolve().parent
OUT = HERE
ROOT = HERE.parents[1]
DATA = HERE/"raw_data"
FIGURES = HERE/"figures"
LOGS = HERE/"logs"


# ---------------------------------------------------------------------------
# Model construction
# ---------------------------------------------------------------------------

def cantilever_mesh(nx=12, ny=8, width=1.0, height=0.6,
                    heights: np.ndarray | None = None):
    mesh = CheckerboardMesh.rectangle(
        width=width, height=height, nx=nx, ny=ny)
    if heights is None:
        return mesh
    reference = mesh.reference_vertices.clone()
    reference[..., 2] = torch.as_tensor(heights, dtype=reference.dtype)
    return CheckerboardMesh(
        reference_vertices=reference, hx=mesh.hx, hy=mesh.hy,
        quads=mesh.quads, boundaries=mesh.boundaries)


def clamped_model(mesh, material, load=None, settings=None):
    boundary = BoundaryConditions.clamped_boundary(mesh, "left")
    load = NodalLoad.zero(mesh) if load is None else load
    settings = (NewtonSettings(residual_tolerance=1.0e-12,
                               maximum_newton_iterations=8)
                if settings is None else settings)
    return CheckerboardShell(mesh, material, boundary, load, settings)


def free_layout(model):
    mask = model.boundary_conditions.free_mask.detach().cpu().numpy()
    flat = np.flatnonzero(mask.reshape(-1))
    lookup = {int(global_index): local for local, global_index in enumerate(flat)}
    return flat, lookup


def node_id(mesh, i, j):
    return j * (mesh.nx + 1) + i


def z_dof(mesh, node):
    return 3 * node + 2


def selection_matrix(model, nodes, component=2):
    """Map free DOFs to selected nodal components."""
    flat, lookup = free_layout(model)
    matrix = np.zeros((len(nodes), flat.size))
    for row, node in enumerate(nodes):
        global_index = 3 * int(node) + component
        if global_index not in lookup:
            raise ValueError(f"node/component {global_index} is fixed")
        matrix[row, lookup[global_index]] = 1.0
    return matrix


def force_matrix(model, nodes, component=2):
    return selection_matrix(model, nodes, component).T


def tangent_at_reference(model):
    n_free = int(model.boundary_conditions.free_mask.sum())
    coordinates = torch.zeros(n_free, dtype=model.reference_vertices.dtype)
    tangent, residual, audit = model.low_rank_tangent(coordinates)
    return tangent, residual, audit


def linear_solve(tangent, rhs):
    solution, audit = tangent.solve(np.asarray(rhs, dtype=float))
    return solution, audit


# ---------------------------------------------------------------------------
# IO and optimisation
# ---------------------------------------------------------------------------

def write_json(path: Path, payload) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, allow_nan=True) + "\n",
                    encoding="utf-8")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        return
    keys = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def lbfgs(function: Callable[[np.ndarray], tuple[float, np.ndarray]],
          initial: np.ndarray, *, max_iterations=40, callback=None):
    from scipy.optimize import minimize
    history: list[dict] = []

    def wrapped(values):
        value, gradient = function(np.asarray(values, dtype=float))
        return float(value), np.asarray(gradient, dtype=float)

    def record(values):
        value, gradient = wrapped(values)
        history.append({
            "iteration": len(history),
            "objective": float(value),
            "gradient_norm": float(np.linalg.norm(gradient)),
            "design": np.asarray(values, dtype=float).tolist(),
        })
        if callback is not None:
            extra = callback(values, value, gradient)
            if isinstance(extra, dict):
                history[-1].update(extra)

    record(initial)
    result = minimize(
        wrapped, np.asarray(initial, dtype=float), jac=True,
        method="L-BFGS-B",
        callback=record,
        options={"maxiter": int(max_iterations), "ftol": 1e-14,
                 "gtol": 1e-10, "maxls": 40})
    record(result.x)
    return result, history


# ---------------------------------------------------------------------------
# Targets
# ---------------------------------------------------------------------------

def gaussian_patch(mesh, x0, y0, sigma_x, sigma_y, total, x_range):
    """Gaussian nodal force distribution on a cantilever plate."""
    vertices = mesh.reference_vertices.detach().cpu().numpy()
    x = vertices[..., 0]
    y = vertices[..., 1]
    weights = np.exp(-((x - x0) / sigma_x) ** 2
                     - ((y - y0) / sigma_y) ** 2)
    mask = (x >= x_range[0]) & (x <= x_range[1])
    weights = weights * mask
    weights = weights / weights.sum()
    return total * weights


def control_grid(mesh, nx_control=4, ny_control=4, width=1.0, height=0.6):
    x = np.linspace(0.0, width, nx_control)
    y = np.linspace(0.0, height, ny_control)
    return x, y


def bilinear_controls(mesh, values, x_control, y_control):
    """Bilinear interpolation of control heights to all mesh vertices."""
    values = np.asarray(values, dtype=float).reshape(
        len(y_control), len(x_control))
    x = mesh.reference_vertices[..., 0].detach().cpu().numpy()
    y = mesh.reference_vertices[..., 1].detach().cpu().numpy()
    ix = np.clip(np.searchsorted(x_control, x) - 1, 0, len(x_control) - 2)
    iy = np.clip(np.searchsorted(y_control, y) - 1, 0, len(y_control) - 2)
    tx = (x - x_control[ix]) / (x_control[ix + 1] - x_control[ix])
    ty = (y - y_control[iy]) / (y_control[iy + 1] - y_control[iy])
    z = ((1 - tx) * (1 - ty) * values[iy, ix]
         + tx * (1 - ty) * values[iy, ix + 1]
         + (1 - tx) * ty * values[iy + 1, ix]
         + tx * ty * values[iy + 1, ix + 1])
    return z


def design_indices(nx_control=4, ny_control=4):
    """Interior control indices; the clamped x=0 column is excluded."""
    indices = []
    for j in range(ny_control):
        for i in range(1, nx_control):
            indices.append((j, i))
    return indices
