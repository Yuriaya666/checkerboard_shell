"""Essential boundary conditions expressed independently of shell mechanics."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


@dataclass(frozen=True)
class BoundaryConditions:
    """Boolean fixed-DOF mask over a structured ``[..., 3]`` vertex field."""

    fixed_mask: torch.Tensor
    prescribed_displacement: torch.Tensor | None = None

    def __post_init__(self):
        if self.prescribed_displacement is not None and self.prescribed_displacement.shape != self.fixed_mask.shape:
            raise ValueError("prescribed displacement must match the fixed mask shape")

    def with_prescribed_displacement(self, values):
        """Return constraints with displacement values at fixed components.

        Values at free components are ignored. Torch values retain their design
        derivative through the position scatter.
        """
        prescribed = torch.as_tensor(values, dtype=torch.float64, device=self.fixed_mask.device)
        return type(self)(self.fixed_mask, prescribed)

    @classmethod
    def fixed_nodes(cls, mesh, node_ids, components=(True, True, True)):
        mask = torch.zeros(mesh.num_vertices, 3, dtype=torch.bool,
                           device=mesh.reference_vertices.device)
        ids = torch.as_tensor(np.asarray(node_ids), dtype=torch.long,
                              device=mask.device)
        selected = torch.as_tensor(components, dtype=torch.bool, device=mask.device)
        mask[ids] = selected
        return cls(mask.reshape(mesh.shape))

    @classmethod
    def fixed_boundary(cls, mesh, name="left", components=(True, True, True)):
        """Fix one named boundary line (a displacement-only constraint)."""
        return cls.fixed_nodes(mesh, mesh.boundaries[name], components)

    @classmethod
    def clamped_boundary(cls, mesh, name="left"):
        """Fix a boundary and its adjacent line to constrain discrete tangent."""
        grid = np.arange(mesh.num_vertices).reshape(mesh.ny + 1, mesh.nx + 1)
        if name == "left":
            node_ids = grid[:, :2].reshape(-1)
        elif name == "right":
            node_ids = grid[:, -2:].reshape(-1)
        elif name == "bottom":
            node_ids = grid[:2, :].reshape(-1)
        elif name == "top":
            node_ids = grid[-2:, :].reshape(-1)
        else:
            raise ValueError(f"unknown structured boundary {name!r}")
        return cls.fixed_nodes(mesh, node_ids)

    @property
    def free_mask(self):
        return ~self.fixed_mask

    def initial_coordinates(self, reference):
        return reference.new_zeros(int(self.free_mask.sum()))

    def positions(self, reference, free_displacement):
        if free_displacement.numel() != int(self.free_mask.sum()):
            raise ValueError("free displacement size does not match boundary conditions")
        displacement = torch.zeros_like(reference)
        displacement[self.free_mask] = free_displacement
        if self.prescribed_displacement is not None:
            displacement[self.fixed_mask] = self.prescribed_displacement.to(
                dtype=reference.dtype, device=reference.device)[self.fixed_mask]
        return reference + displacement

    def coordinates_from_positions(self, reference, positions):
        return (positions - reference)[self.free_mask]
