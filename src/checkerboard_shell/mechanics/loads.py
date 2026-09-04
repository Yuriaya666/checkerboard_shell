"""Parameter-driven dead loads for forward simulation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
import torch


@dataclass(frozen=True)
class NodalLoad:
    """Dead nodal forces with the same structured shape as the mesh."""

    values: torch.Tensor
    name: str = "nodal_load"

    @classmethod
    def zero(cls, mesh):
        return cls(torch.zeros_like(mesh.reference_vertices), name="zero")

    @classmethod
    def total_force(cls, mesh, force, node_ids=None, name="total_force"):
        """Distribute a prescribed resultant force equally over selected nodes."""
        if node_ids is None:
            ids = torch.arange(mesh.num_vertices, device=mesh.reference_vertices.device)
        else:
            ids = torch.as_tensor(np.asarray(node_ids), dtype=torch.long,
                                  device=mesh.reference_vertices.device)
        if ids.numel() == 0:
            raise ValueError("at least one loaded node is required")
        values = torch.zeros(mesh.num_vertices, 3,
                             dtype=mesh.reference_vertices.dtype,
                             device=mesh.reference_vertices.device)
        resultant = torch.as_tensor(force, dtype=values.dtype, device=values.device)
        values[ids] = resultant / ids.numel()
        return cls(values.reshape(mesh.shape), name=name)

    def scaled(self, factor):
        return NodalLoad(float(factor) * self.values, name=self.name)

