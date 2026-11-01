"""Dense reference and matrix-free nonlinear solver interfaces."""

from checkerboard_shell.solver.boundary_conditions import BoundaryConditions
from checkerboard_shell.solver.config import SolverConfig
from checkerboard_shell.solver.continuation import (
    ContinuationResult, adaptive_continuation,
)
from checkerboard_shell.solver.convergence import (
    ConvergenceMonitor, convergence_metrics, is_converged, residual_norms,
)
from checkerboard_shell.solver.newton import NewtonSettings, NewtonSolver, SolverSettings, solve_load_level
from checkerboard_shell.solver.nonlinear import (
    NonlinearProblem, NonlinearResult, solve as solve_nonlinear,
    solve_level, solve_newton_cg, solve_newton_gmres,
)
from checkerboard_shell.solver.result import SimulationResult

__all__ = [
    "BoundaryConditions", "ContinuationResult", "ConvergenceMonitor",
    "NewtonSettings", "NewtonSolver", "NonlinearProblem", "NonlinearResult",
    "SimulationResult", "SolverConfig", "SolverSettings",
    "adaptive_continuation", "convergence_metrics", "is_converged",
    "residual_norms", "solve_load_level", "solve_newton_cg",
    "solve_level", "solve_newton_gmres", "solve_nonlinear",
]
