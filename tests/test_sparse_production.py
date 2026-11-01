"""Regression gates for local boundary storage and exact low-rank tangents."""

from __future__ import annotations

import numpy as np
import torch

from checkerboard_shell.geometry.boundary import (
    boundary_geometry,
    fast_boundary_geometry,
    koiter_targets,
)
from checkerboard_shell.geometry.mesh import CheckerboardMesh
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent


def test_local_boundary_stencil_matches_validated_slow_closure():
    mesh = CheckerboardMesh.rectangle(nx=5, ny=4, width=1.3, height=0.7)
    reference = mesh.reference_vertices.clone()
    reference[..., 2] = 0.11 * reference[..., 0] ** 2 - 0.07 * reference[..., 1] ** 2
    perturbation = 2.0e-3 * torch.sin(
        torch.arange(reference.numel(), dtype=reference.dtype)).reshape_as(reference)
    current = reference + perturbation
    targets = koiter_targets(reference, mesh.hx, mesh.hy)
    slow = boundary_geometry(current, mesh.hx, mesh.hy)
    fast = fast_boundary_geometry(current, mesh.hx, mesh.hy, targets)
    for name in (
        "t1_w", "t2_w", "normal_x", "normal_y", "d1_normal",
        "d2_normal", "first_form_w", "second_form", "x_alternating",
    ):
        torch.testing.assert_close(
            getattr(fast, name), getattr(slow, name), rtol=2.0e-13, atol=2.0e-13)
    sites = (mesh.nx + 1) * (mesh.ny + 1)
    assert targets.reconstruction.indices.shape == (sites, 9)
    dense_old_bytes = 4 * sites * mesh.nx * mesh.ny * reference.element_size()
    assert targets.reconstruction.storage_bytes < dense_old_bytes


def test_exact_nonlinear_local_plus_rank_one_matches_dense_hessian():
    mesh = CheckerboardMesh.rectangle(nx=3, ny=3)
    reference = mesh.reference_vertices.clone()
    reference[..., 2] = 0.08 * (
        reference[..., 0] ** 2 - 0.4 * reference[..., 1] ** 2)
    material = IsotropicMaterial(young=17.0, poisson=0.27, thickness=0.12)
    targets = make_targets(reference, mesh.hx, mesh.hy)
    phase = torch.arange(reference.numel(), dtype=reference.dtype).reshape_as(reference)
    displacement = 3.0e-3 * torch.sin(0.37 * phase)
    tangent, _, _ = assemble_paper_low_rank_tangent(
        displacement, reference, mesh.hx, mesh.hy, targets, material,
        variant="CB-Full", colour_batch=11,
    )

    def energy(flat):
        return paper_components(
            reference + flat.reshape_as(reference), reference,
            mesh.hx, mesh.hy, targets, material,
        )["total"]

    dense = torch.autograd.functional.hessian(
        energy, displacement.reshape(-1), vectorize=True).detach().numpy()
    reconstructed = tangent.local.toarray() + tangent.coefficient * np.outer(
        tangent.vector, tangent.vector)
    np.testing.assert_allclose(reconstructed, dense, rtol=3.0e-10, atol=3.0e-10)
    direction = np.cos(np.arange(dense.shape[0], dtype=float))
    np.testing.assert_allclose(
        tangent.action(direction), dense @ direction,
        rtol=3.0e-10, atol=3.0e-10,
    )


def test_xm_sensitivity_multiplier_only_scales_the_reference_rank_one_term():
    mesh = CheckerboardMesh.rectangle(nx=4, ny=4)
    reference = mesh.reference_vertices
    material = IsotropicMaterial(young=13.0, poisson=0.25, thickness=0.09)
    targets = make_targets(reference, mesh.hx, mesh.hy)
    zero = torch.zeros_like(reference)
    base, _, _ = assemble_paper_low_rank_tangent(
        zero, reference, mesh.hx, mesh.hy, targets, material,
        xm_multiplier=0.0)
    nominal, _, _ = assemble_paper_low_rank_tangent(
        zero, reference, mesh.hx, mesh.hy, targets, material,
        xm_multiplier=1.0)
    doubled, _, _ = assemble_paper_low_rank_tangent(
        zero, reference, mesh.hx, mesh.hy, targets, material,
        xm_multiplier=2.0)
    np.testing.assert_allclose(base.local.toarray(), nominal.local.toarray())
    np.testing.assert_allclose(nominal.local.toarray(), doubled.local.toarray())
    assert base.coefficient == 0.0
    assert doubled.coefficient == 2.0*nominal.coefficient
    np.testing.assert_allclose(nominal.vector, doubled.vector)


def test_cached_coloured_hessian_matches_legacy_at_deformed_states():
    mesh = CheckerboardMesh.rectangle(nx=4, ny=4)
    reference = mesh.reference_vertices.clone()
    reference[..., 2] = 0.08 * reference[..., 0] ** 2
    material = IsotropicMaterial(young=17.0, poisson=0.27, thickness=0.12)
    targets = make_targets(reference, mesh.hx, mesh.hy)
    phase = torch.arange(reference.numel(), dtype=reference.dtype).reshape_as(reference)
    rotation = torch.tensor([[0.,0.,1.],[0.,1.,0.],[-1.,0.,0.]],
                            dtype=reference.dtype)
    rotated = reference @ rotation.T - reference
    rotated_energy = paper_components(reference+rotated, reference,
                                      mesh.hx, mesh.hy, targets, material)["total"]
    assert abs(float(rotated_energy)) < 1e-10
    states = [scale * torch.sin(0.37 * phase)
              for scale in (0.0, 0.003, 0.02)]
    states.append(rotated + 0.003*torch.sin(0.37*phase))
    for displacement in states:
        original, original_force, _ = assemble_paper_low_rank_tangent(
            displacement, reference, mesh.hx, mesh.hy, targets, material,
            colour_batch=7, assembly_strategy="legacy")
        accelerated, accelerated_force, _ = assemble_paper_low_rank_tangent(
            displacement, reference, mesh.hx, mesh.hy, targets, material,
            colour_batch=21, assembly_strategy="cached")
        np.testing.assert_allclose(accelerated_force, original_force,
                                   rtol=1e-12, atol=1e-12)
        np.testing.assert_allclose(accelerated.local.toarray(),
                                   original.local.toarray(), rtol=1e-11, atol=1e-11)
        np.testing.assert_allclose(accelerated.vector, original.vector,
                                   rtol=1e-12, atol=1e-12)
        assert accelerated.coefficient == original.coefficient
        direction = np.cos(np.arange(reference.numel(), dtype=float))
        np.testing.assert_allclose(accelerated.action(direction),
                                   original.action(direction), rtol=1e-11, atol=1e-11)
        skew = accelerated.local-accelerated.local.T
        assert not skew.nnz or np.max(abs(skew.data)) < 1e-12
