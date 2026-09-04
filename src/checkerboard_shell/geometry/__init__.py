"""Mesh construction and B/W geometric operators."""

from checkerboard_shell.geometry.boundary import (
    BoundaryGeometry, BoundaryTargets, fast_boundary_geometry,
    koiter_targets, targets_from_reference,
)
from checkerboard_shell.geometry.core import EdgeMidpoints, TorchGeometry, edge_midpoints, torch_geometry
from checkerboard_shell.geometry.intrinsic import intrinsic_m_coordinate
from checkerboard_shell.geometry.mesh import CheckerboardMesh

__all__ = [
    "BoundaryGeometry", "BoundaryTargets", "CheckerboardMesh",
    "EdgeMidpoints", "TorchGeometry", "edge_midpoints",
    "fast_boundary_geometry", "intrinsic_m_coordinate", "koiter_targets",
    "targets_from_reference", "torch_geometry",
]

