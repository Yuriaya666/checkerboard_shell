"""Public forward-model assembly for a Checkerboard shell."""

from __future__ import annotations

import math
import time

import numpy as np
import scipy.sparse as sp
import torch

from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.loads import NodalLoad
from checkerboard_shell.mechanics.low_rank import (
    LowRankTangent, assemble_paper_low_rank_tangent,
)
from checkerboard_shell.solver.boundary_conditions import BoundaryConditions
from checkerboard_shell.solver.newton import NewtonSettings, NewtonSolver
from checkerboard_shell.solver.result import SimulationResult


class CheckerboardShell:
    """Assemble geometry, material, constraints, load, and the Newton solver."""

    def __init__(self, mesh, material, boundary_conditions=None, load=None,
                 solver_settings=None, extra_stiffness=None,
                 hessian_assembly_strategy="cached"):
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
        self.extra_stiffness = None
        self.hessian_assembly_strategy = hessian_assembly_strategy
        self._extra_torch = None
        if extra_stiffness is not None:
            matrix = sp.csr_matrix(extra_stiffness)
            if matrix.shape != (3*mesh.num_vertices, 3*mesh.num_vertices):
                raise ValueError("extra stiffness must act on all 3V positions")
            self.extra_stiffness = matrix
            self._extra_torch = torch.as_tensor(matrix.toarray(),
                dtype=mesh.reference_vertices.dtype,
                device=mesh.reference_vertices.device)
        self.solver = NewtonSolver(
            NewtonSettings() if solver_settings is None else solver_settings)

    @property
    def reference_vertices(self):
        return self.mesh.reference_vertices

    def positions(self, coordinates):
        return self.boundary_conditions.positions(
            self.reference_vertices, coordinates)

    def energy_components(self, coordinates):
        parts = paper_components(
            self.positions(coordinates), self.reference_vertices,
            self.mesh.hx, self.mesh.hy, self.targets, self.material)
        if self._extra_torch is not None:
            displacement = (self.positions(coordinates)-self.reference_vertices).reshape(-1)
            extra = .5*displacement @ self._extra_torch @ displacement
            parts["reinforcement"] = extra
            parts["total"] = parts["total"]+extra
        return parts

    def potential(self, coordinates, load_factor=1.0):
        positions = self.positions(coordinates)
        internal = paper_components(
            positions, self.reference_vertices, self.mesh.hx, self.mesh.hy,
            self.targets, self.material)["total"]
        if self._extra_torch is not None:
            displacement = (positions-self.reference_vertices).reshape(-1)
            internal = internal + .5*displacement @ self._extra_torch @ displacement
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
        """Dense exact tangent retained as a small-problem reference API."""
        variable = coordinates.detach().clone().requires_grad_(True)
        return torch.autograd.functional.jacobian(
            lambda value: self.residual(
                value, load_factor, create_graph=True),
            variable, vectorize=True)

    def low_rank_tangent(self, coordinates, load_factor=1.0):
        """Production sparse/local plus exact rank-one tangent and residual."""
        displacement = self.positions(coordinates) - self.reference_vertices
        tangent, residual, audit = assemble_paper_low_rank_tangent(
            displacement, self.reference_vertices, self.mesh.hx, self.mesh.hy,
            self.targets, self.material, variant="CB-Full",
            colour_batch=7 if self.hessian_assembly_strategy=="legacy" else 21,
            assembly_strategy=self.hessian_assembly_strategy)
        if self.extra_stiffness is not None:
            tangent = LowRankTangent(
                (tangent.local+self.extra_stiffness).tocsr(),
                tangent.vector,tangent.coefficient)
            residual = residual+self.extra_stiffness @ displacement.detach().cpu().numpy().reshape(-1)
        active = np.flatnonzero(
            self.boundary_conditions.free_mask.reshape(-1).detach().cpu().numpy())
        reduced = LowRankTangent(
            tangent.local[active][:, active].tocsr(),
            tangent.vector[active], tangent.coefficient)
        external = self.load.values.reshape(-1).detach().cpu().numpy()[active]
        return reduced, residual[active]-float(load_factor)*external, audit

    def full_residual_and_reactions(self, coordinates, load_factor=1.0):
        """Return full nodal residual and reaction on each constrained DOF."""
        positions = self.positions(coordinates).detach().clone().requires_grad_(True)
        internal = paper_components(
            positions, self.reference_vertices, self.mesh.hx, self.mesh.hy,
            self.targets, self.material)["total"]
        forces, = torch.autograd.grad(internal, positions)
        if self._extra_torch is not None:
            displacement = (positions-self.reference_vertices).reshape(-1)
            forces = forces+(self._extra_torch @ displacement).reshape_as(forces)
        residual = forces - float(load_factor) * self.load.values
        reaction = torch.zeros_like(residual)
        reaction[self.boundary_conditions.fixed_mask] = residual.detach()[
            self.boundary_conditions.fixed_mask]
        return residual.detach(), reaction

    def solve(self, initial_displacement=None, load_factor=1.0,
              *, backend="sparse_low_rank"):
        """Solve with the production sparse/rank-one or dense reference path."""
        if initial_displacement is None:
            initial = self.boundary_conditions.initial_coordinates(
                self.reference_vertices)
        elif initial_displacement.shape == self.reference_vertices.shape:
            initial = initial_displacement[self.boundary_conditions.free_mask]
        else:
            initial = initial_displacement
        if backend == "dense_reference":
            raw = self.solver.solve(
                lambda coordinates: self.potential(coordinates, load_factor), initial)
        elif backend == "sparse_low_rank":
            raw = self._solve_sparse_low_rank(initial, load_factor)
        else:
            raise ValueError(f"unknown solver backend {backend!r}")
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

    def _solve_sparse_low_rank(self, initial, load_factor):
        settings = self.solver.settings
        coordinates = initial.detach().clone()
        history = []
        started = time.perf_counter()
        status = "MAXIMUM_ITERATIONS"
        for iteration in range(settings.maximum_newton_iterations + 1):
            tangent, residual, audit = self.low_rank_tangent(
                coordinates, load_factor)
            residual_inf = float(np.max(np.abs(residual))) if residual.size else 0.0
            residual_l2 = float(np.linalg.norm(residual))
            value = float(self.potential(coordinates, load_factor).detach())
            history.append({
                "iteration": iteration, "potential": value,
                "residual_inf": residual_inf, "residual_l2": residual_l2,
                "local_tangent_nnz": int(tangent.local.nnz),
                "rank_one_storage": tangent.rank_one_storage,
                "assembly_seconds": (audit["graph_seconds"]+audit["hvp_seconds"]),
            })
            if residual_inf <= settings.residual_tolerance:
                status = "CONVERGED"; break
            if iteration == settings.maximum_newton_iterations:
                break
            try:
                direction, solve_audit = tangent.solve(-residual)
            except Exception as error:
                status = f"LINEAR_SOLVE_FAILED:{type(error).__name__}"; break
            directional = float(residual @ direction)
            if not np.all(np.isfinite(direction)) or directional >= 0.0:
                status = "NON_DESCENT_DIRECTION"; break
            alpha = 1.0
            accepted = False
            while alpha >= settings.minimum_alpha:
                trial = coordinates + torch.as_tensor(
                    alpha*direction, dtype=coordinates.dtype,
                    device=coordinates.device)
                try:
                    trial_value = float(self.potential(trial, load_factor).detach())
                except (RuntimeError, ValueError):
                    trial_value = math.inf
                if (math.isfinite(trial_value)
                        and trial_value <= value
                        + settings.armijo*alpha*directional):
                    coordinates = trial.detach(); accepted = True
                    history[-1].update(solve_audit)
                    history[-1]["line_search_alpha"] = alpha
                    break
                alpha *= 0.5
            if not accepted:
                status = "LINE_SEARCH_FAILED"; break
        return {
            "status": status, "coordinates": coordinates,
            "iterations": history[-1]["iteration"], "history": history,
            "runtime_seconds": time.perf_counter()-started,
        }
