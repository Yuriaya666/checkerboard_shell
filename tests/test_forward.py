import torch

from checkerboard_shell import (
    BoundaryConditions, CheckerboardMesh, CheckerboardShell,
    IsotropicMaterial, NodalLoad, NewtonSettings,
)


def test_rectangle_and_stress_free_energy():
    mesh = CheckerboardMesh.rectangle(nx=4, ny=3)
    model = CheckerboardShell(
        mesh, IsotropicMaterial(young=1e5, thickness=0.05),
        BoundaryConditions.clamped_boundary(mesh), NodalLoad.zero(mesh))
    coordinates = model.boundary_conditions.initial_coordinates(mesh.reference_vertices)
    parts = model.energy_components(coordinates)
    assert torch.allclose(parts["total"], torch.zeros((), dtype=mesh.reference_vertices.dtype), atol=1e-12)


def test_zero_load_forward_solver_converges():
    mesh = CheckerboardMesh.rectangle(nx=4, ny=3)
    model = CheckerboardShell(
        mesh, IsotropicMaterial(young=1e5, thickness=0.05),
        BoundaryConditions.clamped_boundary(mesh), NodalLoad.zero(mesh),
        NewtonSettings(maximum_newton_iterations=3))
    result = model.solve()
    assert result.status == "CONVERGED"
    assert torch.allclose(result.deformed_vertices, mesh.reference_vertices)


def test_nonzero_load_forward_solver_converges_and_deforms():
    mesh = CheckerboardMesh.rectangle(width=1.0, height=0.6, nx=4, ny=3)
    load = NodalLoad.total_force(
        mesh, (0.0, 0.0, -1.0e-3), mesh.boundaries["right"])
    model = CheckerboardShell(
        mesh, IsotropicMaterial(young=2e5, thickness=0.04),
        BoundaryConditions.clamped_boundary(mesh), load,
        NewtonSettings(residual_tolerance=1e-9, maximum_newton_iterations=10))
    result = model.solve()
    assert result.status == "CONVERGED"
    assert result.residual_inf < 1e-9
    assert result.max_displacement > 0.0


def test_residual_and_tangent_are_finite_and_consistent_shapes():
    mesh = CheckerboardMesh.rectangle(nx=4, ny=3)
    model = CheckerboardShell(
        mesh, IsotropicMaterial(young=1e5, thickness=0.05),
        BoundaryConditions.clamped_boundary(mesh), NodalLoad.zero(mesh))
    coordinates = model.boundary_conditions.initial_coordinates(mesh.reference_vertices)
    residual = model.residual(coordinates)
    tangent = model.tangent(coordinates)
    assert residual.shape == coordinates.shape
    assert tangent.shape == (coordinates.numel(), coordinates.numel())
    assert torch.isfinite(tangent).all()
    assert torch.allclose(tangent, tangent.T, atol=1e-10)
