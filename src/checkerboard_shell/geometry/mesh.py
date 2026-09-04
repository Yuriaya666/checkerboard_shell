"""Structured geometry used by the Checkerboard shell operators."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


@dataclass(frozen=True)
class CheckerboardMesh:
    """A structured ``(ny + 1) x (nx + 1)`` surface in material order."""

    reference_vertices: torch.Tensor
    hx: float
    hy: float
    quads: np.ndarray
    boundaries: dict[str, np.ndarray]

    @property
    def nx(self) -> int:
        return self.reference_vertices.shape[1] - 1

    @property
    def ny(self) -> int:
        return self.reference_vertices.shape[0] - 1

    @property
    def shape(self) -> torch.Size:
        return self.reference_vertices.shape

    @property
    def num_vertices(self) -> int:
        return self.reference_vertices.shape[0] * self.reference_vertices.shape[1]

    @classmethod
    def rectangle(
        cls,
        width: float = 1.0,
        height: float = 1.0,
        nx: int = 10,
        ny: int = 10,
        origin: tuple[float, float, float] = (0.0, 0.0, 0.0),
        *,
        dtype: torch.dtype = torch.float64,
        device: torch.device | str | None = None,
    ) -> "CheckerboardMesh":
        """Create the regular rectangular grid used by the paper operators."""
        if nx < 3 or ny < 3:
            raise ValueError("Checkerboard boundary closure requires nx >= 3 and ny >= 3")
        if width <= 0.0 or height <= 0.0:
            raise ValueError("width and height must be positive")

        x0, y0, z0 = origin
        x = torch.linspace(x0, x0 + width, nx + 1, dtype=dtype, device=device)
        y = torch.linspace(y0, y0 + height, ny + 1, dtype=dtype, device=device)
        yy, xx = torch.meshgrid(y, x, indexing="ij")
        zz = torch.full_like(xx, float(z0))
        reference = torch.stack((xx, yy, zz), dim=-1)

        quads = []
        for j in range(ny):
            for i in range(nx):
                a = j * (nx + 1) + i
                quads.append((a, a + 1, a + nx + 2, a + nx + 1))
        boundaries = {
            "bottom": np.arange(nx + 1, dtype=int),
            "top": ny * (nx + 1) + np.arange(nx + 1, dtype=int),
            "left": np.arange(ny + 1, dtype=int) * (nx + 1),
            "right": np.arange(ny + 1, dtype=int) * (nx + 1) + nx,
        }
        return cls(
            reference_vertices=reference,
            hx=float(width / nx),
            hy=float(height / ny),
            quads=np.asarray(quads, dtype=int),
            boundaries=boundaries,
        )

