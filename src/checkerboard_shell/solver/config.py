"""Configuration objects for the generic nonlinear equilibrium solvers."""

from __future__ import annotations

from dataclasses import dataclass
from typing import Literal


@dataclass
class SolverConfig:
    """Explicit options shared by HVP Newton--Krylov and continuation.

    The configuration contains algorithmic controls only.  A physical model
    supplies the energy, residual, tangent action, and (optionally) a residual
    scale through :class:`~checkerboard_shell.solver.nonlinear.NonlinearProblem`.
    """

    method: Literal["newton-gmres", "newton-cg"] = "newton-gmres"
    max_iterations: int = 12
    max_krylov_iterations: int = 32
    tolerance: float = 1.0e-7
    relative_tolerance: float | None = None
    gmres_tolerance: float = 0.15
    forcing_floor: float = 0.02
    damping: float = 0.0
    line_search: bool = True
    line_search_merit: Literal["energy", "residual", "mixed"] = "residual"
    max_backtracks: int = 22
    armijo: float = 1.0e-4
    # Reserved extension switch.  The extracted core intentionally does not
    # yet claim the experimental geometry-aware trust-region implementation.
    trust_region: bool = False

    # Adaptive continuation controls.  They are ignored by a one-level solve.
    continuation: bool = False
    initial_load_step: float = 0.1
    min_load_step: float = 0.01
    max_load_step: float = 0.2
    step_growth: float = 1.5
    step_shrink: float = 0.5
    easy_newton_iterations: int = 5
    easy_alpha: float = 0.5
    max_continuation_attempts: int = 80

    def __post_init__(self):
        if self.method not in ("newton-gmres", "newton-cg"):
            raise ValueError("method must be 'newton-gmres' or 'newton-cg'")
        if self.max_iterations < 0 or self.max_krylov_iterations < 1:
            raise ValueError("iteration limits must be non-negative/positive")
        if self.tolerance < 0.0 or (self.relative_tolerance is not None
                                    and self.relative_tolerance < 0.0):
            raise ValueError("tolerances must be non-negative")
        if self.gmres_tolerance <= 0.0 or self.forcing_floor <= 0.0:
            raise ValueError("Krylov tolerances must be positive")
        if self.forcing_floor > self.gmres_tolerance:
            raise ValueError("forcing_floor cannot exceed gmres_tolerance")
        if self.damping < 0.0:
            raise ValueError("damping must be non-negative")
        if self.max_backtracks < 0 or not 0.0 < self.armijo < 1.0:
            raise ValueError("invalid line-search controls")
        if self.line_search_merit not in ("energy", "residual", "mixed"):
            raise ValueError("line_search_merit must be energy, residual, or mixed")
        if self.min_load_step <= 0.0 or self.max_load_step < self.min_load_step:
            raise ValueError("invalid continuation step bounds")
        if self.initial_load_step <= 0.0:
            raise ValueError("initial_load_step must be positive")
        if self.max_continuation_attempts < 1 or self.easy_newton_iterations < 0:
            raise ValueError("invalid continuation iteration limits")
        if self.easy_alpha < 0.0 or self.easy_alpha > 1.0:
            raise ValueError("easy_alpha must lie in [0, 1]")
        if not 0.0 < self.step_shrink < 1.0:
            raise ValueError("step_shrink must lie strictly between zero and one")
        if self.step_growth < 1.0:
            raise ValueError("step_growth must be at least one")
