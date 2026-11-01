import pytest
import torch

from checkerboard_shell.solver import (
    NonlinearProblem,
    SolverConfig,
    adaptive_continuation,
    convergence_metrics,
    is_converged,
    solve_newton_cg,
    solve_newton_gmres,
)


class QuadraticProblem(NonlinearProblem):
    def __init__(self, matrix, right_hand_side):
        self.matrix = matrix
        self.right_hand_side = right_hand_side
        self.parameter = 1.0

    def energy(self, x):
        rhs = self.parameter * self.right_hand_side
        return 0.5 * torch.dot(x, self.matrix @ x) - torch.dot(rhs, x)

    def update_parameter(self, parameter):
        self.parameter = float(parameter)

    def residual_scale(self, parameter=None):
        factor = self.parameter if parameter is None else float(parameter)
        return float((factor * self.right_hand_side).abs().max())


class QuadraticWithCustomHVP(QuadraticProblem):
    def __init__(self, matrix, right_hand_side):
        super().__init__(matrix, right_hand_side)
        self.hvp_calls = 0

    def hvp(self, x, vector):
        del x
        self.hvp_calls += 1
        return self.matrix @ vector


def _problem(kind=QuadraticProblem):
    matrix = torch.tensor(
        [[4.0, 1.0, 0.0], [1.0, 3.0, 0.5], [0.0, 0.5, 2.0]],
        dtype=torch.float64,
    )
    rhs = torch.tensor([1.0, -2.0, 0.5], dtype=torch.float64)
    return kind(matrix, rhs)


def _tight_config(method="newton-gmres"):
    return SolverConfig(
        method=method,
        max_iterations=4,
        max_krylov_iterations=6,
        tolerance=1.0e-11,
        gmres_tolerance=1.0e-12,
        forcing_floor=1.0e-12,
        line_search_merit="mixed",
    )


def test_default_hvp_matches_explicit_quadratic_action():
    problem = _problem()
    x = torch.tensor([0.2, -0.4, 0.7], dtype=torch.float64)
    vector = torch.tensor([-0.3, 0.1, 0.8], dtype=torch.float64)
    assert torch.allclose(
        problem.hvp(x, vector), problem.matrix @ vector,
        atol=1.0e-13, rtol=1.0e-13,
    )


def test_newton_gmres_uses_custom_hvp_and_converges():
    problem = _problem(QuadraticWithCustomHVP)
    initial = torch.zeros(3, dtype=torch.float64)
    result = solve_newton_gmres(problem, initial, _tight_config())
    expected = torch.linalg.solve(problem.matrix, problem.right_hand_side)
    assert result.status == "CONVERGED"
    assert result.hvp_count == problem.hvp_calls > 0
    assert torch.allclose(result.coordinates, expected, atol=1.0e-11, rtol=1.0e-11)


def test_newton_cg_converges_on_spd_problem_without_mutating_config():
    problem = _problem()
    initial = torch.zeros(3, dtype=torch.float64)
    config = _tight_config(method="newton-gmres")
    result = solve_newton_cg(problem, initial, config)
    expected = torch.linalg.solve(problem.matrix, problem.right_hand_side)
    assert result.status == "CONVERGED"
    assert config.method == "newton-gmres"
    assert torch.allclose(result.coordinates, expected, atol=1.0e-11, rtol=1.0e-11)


def test_relative_convergence_criterion():
    residual = torch.tensor([2.0e-5, -1.0e-5], dtype=torch.float64)
    metrics = convergence_metrics(residual, scale=1.0e4)
    assert not is_converged(metrics, absolute_tolerance=1.0e-6)
    assert is_converged(
        metrics, absolute_tolerance=1.0e-6, relative_tolerance=1.0e-8)


def test_adaptive_continuation_reaches_full_parameter_and_grows_step():
    problem = _problem()
    initial = torch.zeros(3, dtype=torch.float64)
    config = _tight_config()
    config.continuation = True
    config.initial_load_step = 0.2
    config.min_load_step = 0.025
    config.max_load_step = 0.45
    config.step_growth = 1.5
    config.easy_newton_iterations = 1
    result = adaptive_continuation(problem, initial, config)
    assert result.status == "CONVERGED"
    assert result.accepted_parameter == 1.0
    assert problem.parameter == 1.0
    assert any(row["step_size_next"] > config.initial_load_step
               for row in result.attempts)
    expected = torch.linalg.solve(problem.matrix, problem.right_hand_side)
    assert torch.allclose(result.coordinates, expected, atol=1.0e-10, rtol=1.0e-10)


def test_linearize_honours_a_custom_residual_graph():
    class CustomResidualProblem(QuadraticProblem):
        def residual(self, x, *, create_graph=False):
            value = 2.0 * (self.matrix @ x - self.right_hand_side)
            return value if create_graph else value.detach()

    problem = _problem(CustomResidualProblem)
    x = torch.tensor([0.2, -0.4, 0.7], dtype=torch.float64)
    vector = torch.tensor([-0.3, 0.1, 0.8], dtype=torch.float64)
    _energy, residual, action = problem.linearize(x)
    assert torch.allclose(
        residual, 2.0 * (problem.matrix @ x - problem.right_hand_side))
    assert torch.allclose(action(vector), 2.0 * problem.matrix @ vector)


def test_newton_cg_accepts_a_generic_preconditioner():
    problem = _problem()
    initial = torch.zeros(3, dtype=torch.float64)
    calls = 0

    def diagonal_inverse(value):
        nonlocal calls
        calls += 1
        return value / torch.diag(problem.matrix)

    result = solve_newton_cg(
        problem, initial, _tight_config("newton-cg"),
        preconditioner=diagonal_inverse)
    expected = torch.linalg.solve(problem.matrix, problem.right_hand_side)
    assert result.status == "CONVERGED"
    assert calls > 0
    assert torch.allclose(result.coordinates, expected, atol=1.0e-11, rtol=1.0e-11)


@pytest.mark.parametrize("values", [
    {"gmres_tolerance": 0.0},
    {"gmres_tolerance": 0.1, "forcing_floor": 0.2},
    {"damping": -1.0},
    {"armijo": 1.0},
    {"initial_load_step": 0.0},
])
def test_solver_config_rejects_invalid_algorithmic_controls(values):
    with pytest.raises(ValueError):
        SolverConfig(**values)
