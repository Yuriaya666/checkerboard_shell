"""Forward-simulation result object."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import numpy as np
import torch


@dataclass
class SimulationResult:
    mesh: object
    status: str
    reference_vertices: torch.Tensor
    deformed_vertices: torch.Tensor
    coordinates: torch.Tensor
    history: list[dict]
    energies: dict[str, float]
    iterations: int
    residual_inf: float
    runtime_seconds: float

    @property
    def displacement(self):
        return self.deformed_vertices - self.reference_vertices

    @property
    def max_displacement(self) -> float:
        return float(torch.linalg.vector_norm(self.displacement, dim=-1).max())

    def save(self, path):
        """Save reference, deformation and displacement in a portable NPZ file."""
        target = Path(path)
        target.parent.mkdir(parents=True, exist_ok=True)
        np.savez(
            target,
            reference=self.reference_vertices.detach().cpu().numpy(),
            deformed=self.deformed_vertices.detach().cpu().numpy(),
            displacement=self.displacement.detach().cpu().numpy(),
            status=np.asarray(self.status),
            residual_inf=np.asarray(self.residual_inf),
        )
        return target

    def plot(self, *, save_path=None, show=False, displacement_scale=1.0):
        """Plot reference and deformed quad meshes without affecting the solver."""
        from checkerboard_shell.visualization.plotting import plot_deformation

        return plot_deformation(
            self, save_path=save_path, show=show,
            displacement_scale=displacement_scale)

