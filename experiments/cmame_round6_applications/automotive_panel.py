"""Automotive panel study: fixture-force design for gap/flush control."""

from __future__ import annotations

import json
from pathlib import Path

import numpy as np

from checkerboard_shell import IsotropicMaterial
from experiments.cmame_round5.common import cantilever_mesh, node_id
from experiments.cmame_round6_applications.common import (
    checkerboard_design, checkerboard_problem, gaussian_weights,
    mitc4_design, mitc4_response, write_json,
)

HERE = Path(__file__).resolve().parent
DATA = HERE/"raw_data"
DATA.mkdir(parents=True, exist_ok=True)

NX, NY = 12, 8
WIDTH, HEIGHT = 1.0, 0.6
THICKNESS = 0.02
YOUNG, POISSON = 2.0e5, 0.3
CURVATURE = 0.04

# Fixture forces: a 3x3 patch near the free top/right corner.
FORCE_IJ = [(i, j) for j in (4, 6, 8) for i in (8, 10, 11)]
# Gap/flush measurements: right edge and the rear half of the top edge.
OBS_IJ = ([(NX, j) for j in range(2, 7)]
          + [(i, NY) for i in range(6, NX)])


def build_case():
    base = cantilever_mesh(NX, NY, WIDTH, HEIGHT)
    vertices = base.reference_vertices.detach().cpu().numpy()
    x = vertices[..., 0]
    y = vertices[..., 1]
    heights = CURVATURE * np.sin(np.pi * y / HEIGHT)
    force_nodes = [node_id(base, i, j) for i, j in FORCE_IJ]
    obs_nodes = [node_id(base, i, j) for i, j in OBS_IJ]
    # A rear-edge clamp: representative of a hinged/attached panel.
    fixed_nodes = [node_id(base, i, j) for j in (0, 1)
                   for i in range(NX + 1)]
    weights = gaussian_weights(
        base, force_nodes, x0=0.90, y0=0.45, sx=0.12, sy=0.12,
        total=0.01, x_range=(0.60, 1.0))
    material = IsotropicMaterial(YOUNG, POISSON, THICKNESS)
    problem = checkerboard_problem(
        nx=NX, ny=NY, width=WIDTH, height=HEIGHT, heights=heights,
        material=material, fixed_nodes=fixed_nodes,
        force_nodes=force_nodes, obs_nodes=obs_nodes, f_true=weights)
    return base, heights, fixed_nodes, force_nodes, obs_nodes, weights, problem


def run():
    base, heights, fixed_nodes, force_nodes, obs_nodes, weights, problem = build_case()
    cb = checkerboard_design(problem)
    mt_target = mitc4_response(
        nx=NX, ny=NY, width=WIDTH, height=HEIGHT, heights=heights,
        thickness=THICKNESS, young=YOUNG, poisson=POISSON,
        fixed_nodes=fixed_nodes, force_nodes=force_nodes,
        forces=weights, obs_nodes=obs_nodes)
    mt = mitc4_design(
        nx=NX, ny=NY, width=WIDTH, height=HEIGHT, heights=heights,
        thickness=THICKNESS, young=YOUNG, poisson=POISSON,
        fixed_nodes=fixed_nodes, force_nodes=force_nodes,
        obs_nodes=obs_nodes, target=mt_target, f_true=weights)
    mt["free_dofs"] = 6 * (NX + 1) * (NY + 1) - 6 * len(fixed_nodes)
    payload = {
        "case": "automotive_panel_fixture",
        "geometry": "rear-edge-clamped cylindrical panel, 12x8",
        "width": WIDTH, "height": HEIGHT, "thickness": THICKNESS,
        "young": YOUNG, "poisson": POISSON, "curvature": CURVATURE,
        "force_nodes": force_nodes, "observation_nodes": obs_nodes,
        "force_xy": [[float(base.reference_vertices[j, i, 0]),
                      float(base.reference_vertices[j, i, 1])]
                     for i, j in FORCE_IJ],
        "observation_xy": [[float(base.reference_vertices[j, i, 0]),
                            float(base.reference_vertices[j, i, 1])]
                           for i, j in OBS_IJ],
        "design_variables": len(force_nodes),
        "free_dofs": int(problem["model"].boundary_conditions.free_mask.sum()),
        "observations": len(obs_nodes),
        "checkerboard": cb,
        "mitc4": mt,
        "target": problem["target"].tolist(),
        "mitc4_target": mt_target.tolist(),
        "target_policy": "each model uses its own target from the same f_true",
        "f_true": weights.tolist(),
    }
    write_json(DATA/"automotive_panel.json", payload)
    return payload


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
