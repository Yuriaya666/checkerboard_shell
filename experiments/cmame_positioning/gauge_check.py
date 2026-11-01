"""Check that a checkerboard-gauge shift changes neither midpoints nor energy."""

from __future__ import annotations

import json

import numpy as np
import torch

from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components

from .diagnostic import RESULTS, grid, midpoint_map


def main():
    n = 8
    q = grid(n, 0.0)
    P, _, _, _ = midpoint_map(n)
    jj, ii = torch.meshgrid(torch.arange(n + 1), torch.arange(n + 1),
                            indexing="ij")
    parity = torch.where((ii + jj) % 2 == 0, 1., -1.).to(q)
    gauge = parity[..., None] * torch.tensor([.007, -.005, .003], dtype=q.dtype)
    deformed = q.clone()
    deformed[..., 2] += .002 * torch.sin(2 * torch.pi * q[..., 0]) * torch.sin(
        torch.pi * q[..., 1])
    material = IsotropicMaterial(2e6, .3, .008)

    original = paper_components(deformed, q, 1 / n, 1 / n,
                                make_targets(q, 1 / n, 1 / n), material)["total"]
    shifted = paper_components(deformed + gauge, q + gauge, 1 / n, 1 / n,
                               make_targets(q + gauge, 1 / n, 1 / n),
                               material)["total"]
    result = {
        "n": n,
        "gauge_shift_m": [.007, -.005, .003],
        "physical_midpoint_change_inf_m": float(np.max(np.abs(
            P @ gauge.numpy().reshape(-1)))),
        "energy_original_J": float(original),
        "energy_gauge_shift_J": float(shifted),
        "energy_abs_difference_J": float(torch.abs(original - shifted)),
    }
    RESULTS.mkdir(exist_ok=True)
    (RESULTS / "gauge_check.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
