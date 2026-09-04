"""Public forward-model assembly for a Checkerboard shell."""

from __future__ import annotations

import torch

from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.loads import NodalLoad
from checkerboard_shell.solver.boundary_conditions import BoundaryConditions
from checkerboard_shell.solver.newton import NewtonSettings, NewtonSolver
from checkerboard_shell.solver.result import SimulationResult


class CheckerboardShell:
    """Assemble geometry, material, constraints, load, and the Newton solver."""

    def __init__(self, mesh, material, boundary_conditions=None, load=None,
                 solver_settings=None):
        self.mesh = mesh
        self.material = material
        self.boundary_conditions = (
            BoundaryConditions.clamped_boundary(mesh)
            if boundary_conditions is None else boundary_conditions)
        self.load = NodalLoad.zero(mesh) if load is None else load
        if self.load.values.shape != mesh.shape:
            raise ValueError("load values must have the same shape as reference vertices")
        self.targets = make_targets(
            mesh.reference_vertices, mesh.hx, mesh.hy)
        self.solver = NewtonSolver(
            NewtonSettings() if solver_settings is None else solver_settings)

    @property
    def reference_vertices(self):
        return self.mesh.reference_vertices

    def positions(self, coordinates):
        return self.boundary_conditions.positions(
            self.reference_vertices, coordinates)

    def energy_components(self, coordinates):
        return paper_components(
            self.positions(coordinates), self.reference_vertices,
            self.mesh.hx, self.mesh.hy, self.targets, self.material)

    def potential(self, coordinates, load_factor=1.0):
        positions = self.positions(coordinates)
        internal = paper_components(
            positions, self.reference_vertices, self.mesh.hx, self.mesh.hy,
            self.targets, self.material)["total"]
        return internal - float(load_factor) * (self.load.values * positions).sum()

    def residual(self, coordinates, load_factor=1.0, *, create_graph=False):
        """Exact reduced residual: the gradient of total potential."""
        variable = coordinates if coordinates.requires_grad else (
            coordinates.detach().clone().requires_grad_(True))
        residual, = torch.autograd.grad(
            self.potential(variable, load_factor), variable,
            create_graph=create_graph)
        return residual

    def tangent(self, coordinates, load_factor=1.0):
        """Exact consistent tangent of :meth:`residual`."""
        variable = coordinates.detach().clone().requires_grad_(True)
        return torch.autograd.functional.jacobian(
            lambda value: self.residual(
                value, load_factor, create_graph=True),
            variable, vectorize=True)

    def solve(self, initial_displacement=None, load_factor=1.0):
        """Run the extracted Newton solver and return deformation results."""
        if initial_displacement is None:
            initial = self.boundary_conditions.initial_coordinates(
                self.reference_vertices)
        elif initial_displacement.shape == self.reference_vertices.shape:
            initial = initial_displacement[self.boundary_conditions.free_mask]
        else:
            initial = initial_displacement
        raw = self.solver.solve(
            lambda coordinates: self.potential(coordinates, load_factor), initial)
        coordinates = raw["coordinates"]
        positions = self.positions(coordinates).detach()
        parts = self.energy_components(coordinates)
        energies = {
            name: float(value.detach()) for name, value in parts.items()
            if isinstance(value, torch.Tensor) and value.ndim == 0
        }
        return SimulationResult(
            mesh=self.mesh,
            status=raw["status"],
            reference_vertices=self.reference_vertices.detach(),
            deformed_vertices=positions,
            coordinates=coordinates.detach(),
            history=raw["history"],
            energies=energies,
            iterations=raw["iterations"],
            residual_inf=raw["history"][-1]["residual_inf"],
            runtime_seconds=raw["runtime_seconds"],
        )

