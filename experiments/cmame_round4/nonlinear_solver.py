"""Reduced dense Newton / load continuation for the nonlinear benchmarks.

Ported from the archived manuscript driver
``CB-shell/experiments/checkerboard_nonlinear_benchmarks/solver.py`` at
commit ``093ad89`` (tag ``v1.0-arxiv``).  The only change is that the
unused pseudo-arc-length import is removed: none of the three benchmarks
requires arc-length continuation.

The solver never changes the shell energy.  A symmetric spectral shift is
used only to obtain a descent direction when the exact tangent is
indefinite; the unshifted potential and residual remain the acceptance and
convergence gates.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
import math
import time
from typing import Callable

import numpy as np
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


def value_gradient_hessian(function: Callable[[torch.Tensor], torch.Tensor],
                           value: torch.Tensor):
    variable = value.detach().clone().requires_grad_(True)
    scalar = function(variable)
    gradient, = torch.autograd.grad(scalar, variable, create_graph=True)
    hessian = torch.autograd.functional.jacobian(
        lambda x: torch.autograd.grad(function(x), x, create_graph=True)[0],
        variable, vectorize=True)
    hessian = 0.5 * (hessian + hessian.T)
    return scalar.detach(), gradient.detach(), hessian.detach()


def _descent_direction(hessian: torch.Tensor, gradient: torch.Tensor,
                       relative_floor: float):
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


def solve_load_level(potential, initial: torch.Tensor,
                     settings: SolverSettings):
    coordinates = initial.detach().clone()
    history = []
    minimum_alpha = 1.0
    last_values = None
    started = time.perf_counter()
    for iteration in range(settings.maximum_newton_iterations + 1):
        scalar, gradient, hessian = value_gradient_hessian(
            potential, coordinates)
        residual_inf = float(torch.linalg.vector_norm(
            gradient, ord=float("inf")))
        history.append({"iteration": iteration, "potential": float(scalar),
                        "residual_inf": residual_inf,
                        "residual_l2": float(torch.linalg.vector_norm(gradient))})
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
            trial_gradient, = torch.autograd.grad(
                trial_scalar, trial_variable)
            trial_residual = float(torch.linalg.vector_norm(
                trial_gradient, ord=float("inf")))
            energy_descent = (directional < 0.0 and math.isfinite(trial_value)
                              and trial_value <= float(scalar)
                              + settings.armijo*alpha*directional)
            residual_descent = trial_residual < (
                1.0-1e-4*alpha)*residual_inf
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


def adaptive_continuation(make_potential, initial: torch.Tensor,
                          on_accept, settings=SolverSettings(),
                          maximum_load=1.0,
                          mandatory_loads: tuple[float, ...] = ()) -> dict:
    """Load continuation with rollback and automatic step adaptation."""
    coordinates = initial.detach().clone()
    load = 0.0
    step = min(settings.initial_step, maximum_load)
    accepted = []
    rejected = []
    gates = sorted(set(float(value) for value in mandatory_loads
                       if 0.0 < value <= maximum_load))
    while load < maximum_load - 1e-14:
        target = min(maximum_load, load + step)
        future = [value for value in gates if value > load+1e-14]
        if future:
            target = min(target, future[0])
        potential = make_potential(target, coordinates.detach().clone())
        result = solve_load_level(potential, coordinates, settings)
        if result["status"] == "CONVERGED":
            coordinates = result["coordinates"].detach()
            load = target
            row = on_accept(load, coordinates, result, step)
            accepted.append(row)
            if result["iterations"] <= settings.easy_iteration_count:
                step = min(settings.maximum_step, settings.growth*step)
        else:
            rejected.append({"attempted_load": target, "accepted_load": load,
                             "attempted_step": step, "status": result["status"],
                             "iterations": result["iterations"]})
            step *= settings.shrink
            if step < settings.minimum_step - 1e-15:
                return {"status": "CONTINUATION_FAILED",
                        "coordinates": coordinates,
                        "last_converged_load": load, "accepted": accepted,
                        "rejected": rejected, "settings": asdict(settings),
                        "arc_length_recommended": True}
    return {"status": "CONVERGED", "coordinates": coordinates,
            "last_converged_load": load, "accepted": accepted,
            "rejected": rejected, "settings": asdict(settings),
            "arc_length_recommended": False}


def interpolate_curve(rows: list[dict], loads: np.ndarray,
                      keys: tuple[str, ...]):
    x = np.asarray([0.0] + [row["load_factor"] for row in rows])
    result = {}
    for key in keys:
        y = np.asarray([0.0] + [row[key] for row in rows])
        result[key] = np.interp(loads, x, y)
    return result
