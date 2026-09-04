"""Minimal, parameterized Checkerboard Shell forward simulation framework."""

from checkerboard_shell.geometry.mesh import CheckerboardMesh
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.loads import NodalLoad
from checkerboard_shell.mechanics.model import CheckerboardShell
from checkerboard_shell.solver.boundary_conditions import BoundaryConditions
from checkerboard_shell.solver.newton import NewtonSettings

__all__ = [
    "BoundaryConditions", "CheckerboardMesh", "CheckerboardShell",
    "IsotropicMaterial", "NodalLoad", "NewtonSettings",
]

