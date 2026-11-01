"""Generic matrix-free nonlinear equilibrium solvers.

The algorithms in this module know nothing about Checkerboard geometry or FPC
loads.  A model exposes an energy, its gradient/residual, and a Hessian-vector
product through :class:`NonlinearProblem`.  The default HVP is obtained by
autodiff; models with a cheaper tangent action can override it.
"""

from __future__ import annotations

from dataclasses import dataclass, field, replace
import math
import time
from typing import Callable

import torch

from checkerboard_shell.solver.config import SolverConfig
from checkerboard_shell.solver.convergence import (
    convergence_metrics, is_converged, residual_norms,
)


class NonlinearProblem:
    """Minimal model interface consumed by the nonlinear solver.

    ``update_parameter`` is intentionally a no-op: a load continuation driver
    can use it for a model-specific scalar/vector parameter without requiring
    the core solver to know what that parameter means.
    """

    def energy(self, x: torch.Tensor) -> torch.Tensor:
        raise NotImplementedError

    def residual(self, x: torch.Tensor, *, create_graph: bool = False) -> torch.Tensor:
        variable = x if x.requires_grad else x.detach().clone().requires_grad_(True)
        return torch.autograd.grad(
            self.energy(variable), variable, create_graph=create_graph)[0]

    def gradient(self, x: torch.Tensor, *, create_graph: bool = False) -> torch.Tensor:
        """Alias for residual, useful for models using gradient terminology."""
        return self.residual(x, create_graph=create_graph)

    def hvp(self, x: torch.Tensor, vector: torch.Tensor) -> torch.Tensor:
        """Default Hessian-vector product with no Hessian materialization."""
        variable = x.detach().clone().requires_grad_(True)
        residual = self.residual(variable, create_graph=True)
        action, = torch.autograd.grad(residual, variable, vector)
        return action.detach()

    def linearize(self, x: torch.Tensor):
        """Build one residual graph and return ``(energy, residual, action)``."""
        # Honour a model-supplied tangent action.  The base implementation
        # below is faster for the common autodiff case because it retains one
        # graph for all Krylov products at the current Newton state.
        if type(self).hvp is not NonlinearProblem.hvp:
            state = x.detach()
            scalar = self.energy(state).detach()
            residual = self.residual(state).detach()

            def custom_action(vector: torch.Tensor) -> torch.Tensor:
                return self.hvp(state, vector).detach()

            return scalar, residual, custom_action
        if (type(self).residual is not NonlinearProblem.residual
                or type(self).gradient is not NonlinearProblem.gradient):
            variable = x.detach().clone().requires_grad_(True)
            scalar = self.energy(variable)
            residual = (self.residual(variable, create_graph=True)
                        if type(self).residual is not NonlinearProblem.residual
                        else self.gradient(variable, create_graph=True))

            def residual_action(vector: torch.Tensor) -> torch.Tensor:
                value, = torch.autograd.grad(
                    residual, variable, vector, retain_graph=True)
                return value.detach()

            return scalar.detach(), residual.detach(), residual_action
        variable = x.detach().clone().requires_grad_(True)
        scalar = self.energy(variable)
        residual = torch.autograd.grad(scalar, variable, create_graph=True)[0]

        def action(vector: torch.Tensor) -> torch.Tensor:
            value, = torch.autograd.grad(
                residual, variable, vector, retain_graph=True)
            return value.detach()

        return scalar.detach(), residual.detach(), action

    def update_parameter(self, parameter) -> None:
        """Update the continuation parameter; default models need no action."""
        del parameter

    def residual_scale(self, parameter=None) -> float:
        """Scale for an optional relative residual criterion."""
        del parameter
        return 1.0


@dataclass
class NonlinearResult:
    status: str
    coordinates: torch.Tensor
    iterations: int
    gmres_iterations: int
    hvp_count: int
    gradient_evaluations: int
    line_search_rejections: int
    residual_inf: float
    residual_l2: float
    history: list[dict] = field(default_factory=list)
    energy: float = math.inf
    runtime_seconds: float = 0.0
    minimum_alpha: float = 1.0
    linear_history: list[dict] = field(default_factory=list)

    @property
    def krylov_iterations(self) -> int:
        """Model-neutral alias for the retained prototype counter name."""
        return self.gmres_iterations


def _safe_float(value) -> float:
    try:
        result = float(value.detach() if isinstance(value, torch.Tensor) else value)
    except (TypeError, ValueError, RuntimeError):
        return math.inf
    return result if math.isfinite(result) else math.inf


def _gmres(operator: Callable[[torch.Tensor], torch.Tensor], rhs: torch.Tensor,
           maximum: int, tolerance: float,
           preconditioner: Callable[[torch.Tensor], torch.Tensor] | None = None):
    """Deterministic unrestarted left-preconditioned GMRES."""
    apply = (lambda value: value) if preconditioner is None else preconditioner
    transformed_rhs = apply(rhs)
    beta = torch.linalg.vector_norm(transformed_rhs)
    if not torch.isfinite(beta) or float(beta) == 0.0:
        return torch.zeros_like(rhs), 0, 0.0
    basis = [transformed_rhs / beta]
    hessenberg = torch.zeros(
        (maximum + 1, maximum), dtype=rhs.dtype, device=rhs.device)
    target = torch.zeros(maximum + 1, dtype=rhs.dtype, device=rhs.device)
    target[0] = beta
    solution = torch.zeros_like(rhs)
    relative = math.inf
    for column in range(maximum):
        work = apply(operator(basis[column]))
        for row in range(column + 1):
            coefficient = torch.dot(basis[row], work)
            hessenberg[row, column] = coefficient
            work = work - coefficient * basis[row]
        next_norm = torch.linalg.vector_norm(work)
        hessenberg[column + 1, column] = next_norm
        matrix = hessenberg[:column + 2, :column + 1]
        coefficients = torch.linalg.lstsq(
            matrix, target[:column + 2, None]).solution[:, 0]
        solution = sum(
            (coefficients[index] * basis[index]
             for index in range(column + 1)), torch.zeros_like(rhs))
        relative = float(torch.linalg.vector_norm(
            target[:column + 2] - matrix @ coefficients)) / max(float(beta), 1.0e-30)
        if relative <= tolerance or float(next_norm) <= 1.0e-14:
            return solution, column + 1, relative
        basis.append(work / next_norm)
    return solution, maximum, relative


def _cg(operator, rhs, maximum, tolerance, preconditioner=None):
    apply = (lambda value: value) if preconditioner is None else preconditioner
    solution = torch.zeros_like(rhs)
    residual = rhs.clone()
    preconditioned = apply(residual)
    direction = preconditioned.clone()
    rr = torch.dot(residual, preconditioned)
    scale = max(float(torch.linalg.vector_norm(rhs)), 1.0e-30)
    relative = 1.0
    for iteration in range(maximum):
        action = operator(direction)
        denominator = torch.dot(direction, action)
        threshold = 1.0e-14 * max(float(
            torch.linalg.vector_norm(direction)
            * torch.linalg.vector_norm(action)), 1.0e-30)
        if not torch.isfinite(denominator) or float(denominator) <= threshold:
            # On first-step negative curvature the right-hand side is a valid
            # descent direction; later, keep the truncated Krylov solution.
            return (preconditioned if iteration == 0 else solution), iteration + 1, relative
        alpha = rr / denominator
        solution = solution + alpha * direction
        residual = residual - alpha * action
        relative = float(torch.linalg.vector_norm(residual)) / scale
        if relative <= tolerance:
            return solution, iteration + 1, relative
        next_preconditioned = apply(residual)
        next_rr = torch.dot(residual, next_preconditioned)
        if not torch.isfinite(next_rr) or float(next_rr) <= 0.0:
            return solution, iteration + 1, relative
        direction = next_preconditioned + (next_rr / rr) * direction
        preconditioned = next_preconditioned
        rr = next_rr
    return solution, maximum, relative


def _line_search(problem, x, direction, scalar, residual, config,
                 minimum_alpha):
    directional = float(torch.dot(residual, direction))
    current = float(residual.abs().max())
    alpha, rejections = 1.0, 0
    for _ in range(config.max_backtracks + 1):
        trial = x + alpha * direction
        try:
            trial_energy = _safe_float(problem.energy(trial))
            trial_residual = problem.residual(trial)
            trial_norm = float(trial_residual.abs().max())
            finite = math.isfinite(trial_energy) and math.isfinite(trial_norm)
        except (RuntimeError, ValueError, torch.linalg.LinAlgError):
            trial_energy, trial_norm, finite = math.inf, math.inf, False
        energy_ok = (directional < 0.0 and finite
                     and trial_energy <= _safe_float(scalar)
                     + config.armijo * alpha * directional)
        residual_ok = finite and trial_norm < current * (1.0 - 1.0e-4 * alpha)
        accept = ((config.line_search_merit in ("energy", "mixed") and energy_ok)
                  or (config.line_search_merit in ("residual", "mixed")
                      and residual_ok))
        if accept:
            return trial.detach(), alpha, rejections, min(minimum_alpha, alpha), True
        alpha *= 0.5
        rejections += 1
    return x, alpha, rejections, min(minimum_alpha, alpha), False


def solve_level(problem: NonlinearProblem, initial: torch.Tensor,
                config: SolverConfig | None = None, *, parameter=None,
                preconditioner=None) -> NonlinearResult:
    """Solve one parameterized equilibrium level with HVP Newton--Krylov."""
    config = SolverConfig() if config is None else config
    if config.trust_region:
        raise NotImplementedError(
            "the trust-region switch is reserved; use the validated line-search core")
    if parameter is not None:
        problem.update_parameter(parameter)
    started = time.perf_counter()
    x = initial.detach().clone()
    scale = problem.residual_scale(parameter)
    history, linear_history = [], []
    hvp_count = gmres_total = gradient_evaluations = rejections = 0
    minimum_alpha, status, energy = 1.0, "MAXIMUM_ITERATIONS", math.inf
    for iteration in range(config.max_iterations + 1):
        scalar, residual_graph, action = problem.linearize(x)
        gradient_evaluations += 1
        residual_inf, residual_l2 = residual_norms(residual_graph)
        metrics = convergence_metrics(residual_graph, scale)
        energy = _safe_float(scalar)
        record = {"iteration": iteration, "energy": energy, **metrics}
        history.append(record)
        if is_converged(metrics, config.tolerance, config.relative_tolerance):
            status = "CONVERGED"
            break
        if iteration == config.max_iterations:
            break

        def hvp(vector):
            nonlocal hvp_count
            hvp_count += 1
            return action(vector)

        def operator(vector):
            return hvp(vector) + config.damping * vector

        forcing = config.gmres_tolerance
        if len(history) >= 2:
            ratio = residual_inf / max(float(history[-2]["residual_inf"]), 1.0e-30)
            forcing = min(forcing, max(config.forcing_floor, math.sqrt(abs(ratio))))
        forcing = min(forcing, max(
            config.forcing_floor,
            0.1 * math.sqrt(residual_inf / max(config.tolerance, 1.0e-30))))
        rhs = -residual_graph
        if config.method == "newton-gmres":
            direction, used, linear_relative = _gmres(
                operator, rhs, config.max_krylov_iterations, forcing,
                preconditioner)
        else:
            direction, used, linear_relative = _cg(
                operator, rhs, config.max_krylov_iterations, forcing,
                preconditioner)
        gmres_total += used
        linear_history.append({"iteration": iteration,
                               "krylov_iterations": used,
                               "relative_residual": linear_relative,
                               "forcing_tolerance": forcing})
        if (not torch.isfinite(direction).all()
                or float(torch.dot(residual_graph, direction)) >= 0.0):
            direction = -residual_graph
        if config.line_search:
            x, alpha, rejected, minimum_alpha, accepted = _line_search(
                problem, x, direction, scalar, residual_graph, config,
                minimum_alpha)
            rejections += rejected
            if not accepted:
                status = "LINE_SEARCH_FAILED"
                break
        else:
            x = (x + direction).detach()
    final_residual = problem.residual(x)
    final_inf, final_l2 = residual_norms(final_residual)
    final_energy = _safe_float(problem.energy(x))
    return NonlinearResult(
        status=status, coordinates=x.detach(), iterations=len(history) - 1,
        gmres_iterations=gmres_total, hvp_count=hvp_count,
        gradient_evaluations=gradient_evaluations,
        line_search_rejections=rejections, residual_inf=final_inf,
        residual_l2=final_l2, history=history, energy=final_energy,
        runtime_seconds=time.perf_counter() - started,
        minimum_alpha=minimum_alpha, linear_history=linear_history)


def solve(problem: NonlinearProblem, initial: torch.Tensor,
          config: SolverConfig | None = None, *, parameter=None,
          preconditioner=None):
    """Dispatch to one corrector or configured adaptive continuation."""
    config = SolverConfig() if config is None else config
    if config.continuation:
        if parameter is not None:
            raise ValueError(
                "parameter is set by adaptive continuation; do not pass it explicitly")
        if preconditioner is not None:
            raise ValueError(
                "continuation rebuilds correctors; attach preconditioning in the problem")
        from checkerboard_shell.solver.continuation import adaptive_continuation
        return adaptive_continuation(problem, initial, config)
    return solve_level(
        problem, initial, config, parameter=parameter,
        preconditioner=preconditioner)


def solve_newton_gmres(problem, initial, config=None, **kwargs):
    config = SolverConfig(method="newton-gmres") if config is None else config
    config = replace(config, method="newton-gmres")
    return solve_level(problem, initial, config, **kwargs)


def solve_newton_cg(problem, initial, config=None, **kwargs):
    config = SolverConfig(method="newton-cg") if config is None else config
    config = replace(config, method="newton-cg")
    return solve_level(problem, initial, config, **kwargs)


# Short internal names are intentionally private; continuation imports the
# public ``solve`` entry point so applications can use one stable API.
__all__ = ["NonlinearProblem", "NonlinearResult", "solve", "solve_level",
           "solve_newton_gmres", "solve_newton_cg"]
