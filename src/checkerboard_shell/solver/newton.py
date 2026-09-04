"""The energy-consistent Newton/line-search solver from the paper code."""

from __future__ import annotations

from dataclasses import dataclass
import math
import time

import torch


@dataclass
class SolverSettings:
    initial_step: float = 0.02
    minimum_step: float = 0.00125
    maximum_step: float = 0.05
    residual_tolerance: float = 1.0e-7
    maximum_newton_iterations: int = 24
    minimum_alpha: float = 2.0**-20
    armijo: float = 1.0e-4
    easy_iteration_count: int = 4
    growth: float = 1.4
    shrink: float = 0.5
    eigenvalue_relative_floor: float = 1.0e-10


# Public name used by the minimal forward API.
NewtonSettings = SolverSettings


def value_gradient_hessian(function, value):
    variable = value.detach().clone().requires_grad_(True)
    scalar = function(variable)
    gradient, = torch.autograd.grad(scalar, variable, create_graph=True)
    hessian = torch.autograd.functional.jacobian(
        lambda x: torch.autograd.grad(function(x), x, create_graph=True)[0],
        variable, vectorize=True)
    hessian = 0.5 * (hessian + hessian.T)
    return scalar.detach(), gradient.detach(), hessian.detach()


def _descent_direction(hessian, gradient, relative_floor):
    values, vectors = torch.linalg.eigh(hessian)
    scale = max(float(values.abs().max()), 1.0)
    floor = relative_floor * scale
    retained = values.abs() > floor
    direction = -(vectors[:, retained] @ (
        (vectors[:, retained].T @ gradient) / values[retained]))
    if (not torch.isfinite(direction).all()
            or float(gradient @ direction) >= 0.0):
        shifted = torch.clamp(values, min=floor)
        direction = -(vectors @ ((vectors.T @ gradient) / shifted))
    return direction, values


def solve_load_level(potential, initial: torch.Tensor, settings: SolverSettings):
    """Solve one dead-load level while preserving the supplied potential."""
    coordinates = initial.detach().clone()
    history = []
    minimum_alpha = 1.0
    last_values = None
    started = time.perf_counter()
    for iteration in range(settings.maximum_newton_iterations + 1):
        scalar, gradient, hessian = value_gradient_hessian(potential, coordinates)
        residual_inf = float(torch.linalg.vector_norm(gradient, ord=float("inf")))
        residual_l2 = float(torch.linalg.vector_norm(gradient))
        history.append({"iteration": iteration, "potential": float(scalar),
                        "residual_inf": residual_inf, "residual_l2": residual_l2})
        if residual_inf <= settings.residual_tolerance:
            return {"status": "CONVERGED", "coordinates": coordinates,
                    "iterations": iteration, "history": history,
                    "minimum_alpha": minimum_alpha,
                    "tangent_eigenvalues": last_values,
                    "runtime_seconds": time.perf_counter()-started}
        if iteration == settings.maximum_newton_iterations:
            break
        direction, values = _descent_direction(
            hessian, gradient, settings.eigenvalue_relative_floor)
        last_values = values.detach()
        directional = float(gradient @ direction)
        alpha = 1.0
        accepted = False
        while alpha >= settings.minimum_alpha:
            trial = coordinates + alpha * direction
            try:
                trial_value = float(potential(trial).detach())
            except (RuntimeError, ValueError):
                trial_value = math.inf
            trial_variable = trial.detach().clone().requires_grad_(True)
            trial_scalar = potential(trial_variable)
            trial_gradient, = torch.autograd.grad(trial_scalar, trial_variable)
            trial_residual = float(torch.linalg.vector_norm(
                trial_gradient, ord=float("inf")))
            energy_descent = (directional < 0.0 and math.isfinite(trial_value)
                              and trial_value <= float(scalar)
                              + settings.armijo*alpha*directional)
            residual_descent = trial_residual < (1.0-1e-4*alpha)*residual_inf
            if math.isfinite(trial_value) and (energy_descent or residual_descent):
                accepted = True
                coordinates = trial.detach()
                minimum_alpha = min(minimum_alpha, alpha)
                break
            alpha *= 0.5
        if not accepted:
            return {"status": "LINE_SEARCH_FAILED", "coordinates": coordinates,
                    "iterations": iteration+1, "history": history,
                    "minimum_alpha": min(minimum_alpha, alpha),
                    "tangent_eigenvalues": last_values,
                    "runtime_seconds": time.perf_counter()-started}
    return {"status": "MAXIMUM_ITERATIONS", "coordinates": coordinates,
            "iterations": settings.maximum_newton_iterations,
            "history": history, "minimum_alpha": minimum_alpha,
            "tangent_eigenvalues": last_values,
            "runtime_seconds": time.perf_counter()-started}


class NewtonSolver:
    """Small object wrapper convenient for ``CheckerboardShell.solve()``."""

    def __init__(self, settings=None):
        self.settings = NewtonSettings() if settings is None else settings

    def solve(self, potential, initial):
        return solve_load_level(potential, initial, self.settings)
