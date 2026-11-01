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

    @classmethod
    def quad_patch_force(cls, mesh, cell_i, cell_j, local_x, local_y,
                         force, name="patch_force"):
        """Transfer one resultant through a quad attachment point.

        Bilinear weights preserve the resultant and moment about any origin.
        ``local_x`` and ``local_y`` are in [0,1] within the cell.
        """
        if not (0 <= cell_i < mesh.nx and 0 <= cell_j < mesh.ny):
            raise ValueError("attachment cell outside mesh")
        if not (0 <= local_x <= 1 and 0 <= local_y <= 1):
            raise ValueError("attachment coordinate outside cell")
        values = torch.zeros_like(mesh.reference_vertices)
        fx = torch.as_tensor(force, dtype=values.dtype, device=values.device)
        for di, dj, weight in (
            (0, 0, (1-local_x)*(1-local_y)),
            (1, 0, local_x*(1-local_y)),
            (1, 1, local_x*local_y),
            (0, 1, (1-local_x)*local_y),
        ):
            values[cell_j+dj, cell_i+di] += weight*fx
        return cls(values, name=name)

    @classmethod
    def dead_pressure(cls, mesh, pressure, name="dead_pressure"):
        """Consistent corner forces for reference-normal cell pressure.

        ``pressure`` may be a scalar or one value per cell. The pressure is
        fixed in the reference configuration; it is not a follower load.
        """
        field = torch.as_tensor(pressure, dtype=mesh.reference_vertices.dtype,
                                device=mesh.reference_vertices.device)
        if field.ndim == 0:
            field = field.expand(mesh.ny, mesh.nx)
        if field.shape != (mesh.ny, mesh.nx):
            raise ValueError("pressure must have one value per cell")
        vertices = mesh.reference_vertices
        values = torch.zeros_like(vertices)
        for j in range(mesh.ny):
            for i in range(mesh.nx):
                p00, p10 = vertices[j, i], vertices[j, i+1]
                p01, p11 = vertices[j+1, i], vertices[j+1, i+1]
                area_vector = 0.5*(torch.linalg.cross(p10-p00, p01-p00)
                                   + torch.linalg.cross(p11-p01, p11-p10))
                corner_force = field[j, i]*area_vector/4.0
                for di, dj in ((0,0), (1,0), (1,1), (0,1)):
                    values[j+dj, i+di] += corner_force
        return cls(values, name=name)

