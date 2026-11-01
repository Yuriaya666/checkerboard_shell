"""Adaptive rollback/grow continuation for generic nonlinear problems."""

from __future__ import annotations

from dataclasses import dataclass, field

import torch

from checkerboard_shell.solver.config import SolverConfig
from checkerboard_shell.solver.nonlinear import (
    NonlinearProblem, NonlinearResult, solve_level,
)


@dataclass
class ContinuationResult:
    status: str
    coordinates: torch.Tensor
    accepted_parameter: float
    attempts: list[dict] = field(default_factory=list)
    inner_history: list[dict] = field(default_factory=list)
    rollback_count: int = 0
    corrector_results: list[NonlinearResult] = field(default_factory=list)


def adaptive_continuation(problem: NonlinearProblem, initial: torch.Tensor,
                          config: SolverConfig | None = None) -> ContinuationResult:
    """Track parameter ``0 -> 1`` with rollback and adaptive step sizes."""
    config = SolverConfig(continuation=True) if config is None else config
    q = initial.detach().clone()
    parameter, step = 0.0, float(config.initial_load_step)
    problem.update_parameter(parameter)
    easy_count = rollback_count = 0
    attempts, inner, correctors = [], [], []
    while parameter < 1.0 - 1.0e-12 and len(attempts) < config.max_continuation_attempts:
        start = parameter
        target = min(1.0, parameter + step)
        result: NonlinearResult = solve_level(
            problem, q, config, parameter=target)
        correctors.append(result)
        accepted = result.status == "CONVERGED"
        if accepted:
            q = result.coordinates
            parameter = target
            easy = (result.iterations <= config.easy_newton_iterations
                    and result.minimum_alpha >= config.easy_alpha)
            easy_count = easy_count + 1 if easy else 0
            if easy_count >= 2:
                step = min(config.max_load_step, step * config.step_growth)
                easy_count = 0
        else:
            rollback_count += 1
            step *= config.step_shrink
            easy_count = 0
        row = {
            "attempt": len(attempts) + 1,
            "parameter_start": start, "parameter_target": target,
            "parameter_accepted": parameter, "step_size_next": step,
            "accepted": accepted, "status": result.status,
            "iterations": result.iterations,
            "krylov_iterations": result.gmres_iterations,
            "hvp_count": result.hvp_count,
            "residual_inf": result.residual_inf,
            "residual_l2": result.residual_l2,
            "energy": result.energy,
            "runtime_seconds": result.runtime_seconds,
        }
        attempts.append(row)
        inner.extend({"attempt": len(attempts), "parameter": target, **item}
                     for item in result.history)
        if step < config.min_load_step:
            break
    status = "CONVERGED" if parameter >= 1.0 - 1.0e-12 else "CONTINUATION_FAILED"
    # A failed attempt may have left an application model at its rejected
    # target.  Restore the last accepted parameter together with its state.
    problem.update_parameter(parameter)
    return ContinuationResult(
        status=status, coordinates=q, accepted_parameter=parameter,
        attempts=attempts, inner_history=inner,
        rollback_count=rollback_count, corrector_results=correctors)


__all__ = ["ContinuationResult", "adaptive_continuation"]
