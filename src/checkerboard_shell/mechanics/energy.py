"""Frozen Checkerboard shell energy, including ``X_M`` and ``X_W_rel``."""

from __future__ import annotations

from dataclasses import dataclass

import torch

from checkerboard_shell.geometry.boundary import (
    BoundaryTargets, fast_boundary_geometry, koiter_targets,
    transported_boundary_x,
)
from checkerboard_shell.geometry.intrinsic import intrinsic_m_coordinate
from checkerboard_shell.material.isotropic import lame_surface_coefficients, stvk_density


@dataclass
class EnergyTargets:
    """Cached reference quantities; detached so they remain stress-free."""

    checkerboard: object


def make_targets(reference, hx, hy) -> EnergyTargets:
    return EnergyTargets(checkerboard=koiter_targets(reference, hx, hy))


def paper_components(vertices, reference, hx, hy, targets, material,
                     *, include_xw=True):
    """Evaluate the accepted paper energy without altering its coefficients."""
    alpha, mu = lame_surface_coefficients(material.young, material.poisson)
    boundary = targets.checkerboard.boundary
    geometry = fast_boundary_geometry(vertices, hx, hy, targets.checkerboard)
    delta_a_b = geometry.bulk.first_form - boundary.first_form_b
    membrane_density_b = stvk_density(
        delta_a_b, boundary.first_form_b, alpha, mu)
    membrane_b = (float(material.thickness) / 4.0) * (
        boundary.area_b * membrane_density_b).sum()

    xm = intrinsic_m_coordinate(
        vertices, reference, hx, hy, boundary.area_b)
    membrane_xm = (float(material.thickness) / 4.0) * (8.0 * float(mu)) * xm.square()

    delta_b_w = geometry.second_form - boundary.second_form
    bending_density = stvk_density(
        delta_b_w, boundary.first_form_w, alpha, mu)
    bending = (float(material.thickness) ** 3 / 12.0) * (
        boundary.dual_fraction * boundary.area_w * bending_density).sum()

    if bool(torch.all(boundary.x_alternating == 0.0)):
        relative_xw = geometry.x_alternating
        rotation = deformation = None
    else:
        relative_xw, rotation, deformation = transported_boundary_x(
            geometry, boundary)
    if include_xw:
        bending_compatibility = (
            (float(material.thickness) ** 3 / 12.0) * float(mu)
            * (boundary.gamma_w[1:-1, 1:-1]
               * relative_xw[1:-1, 1:-1].square().sum(-1)).sum())
    else:
        bending_compatibility = vertices.new_zeros(())
    total = membrane_b + membrane_xm + bending + bending_compatibility
    return {
        "membrane_b": membrane_b,
        "membrane_xm": membrane_xm,
        "membrane": membrane_b + membrane_xm,
        "bending": bending,
        "bending_compatibility": bending_compatibility,
        "total": total,
        "X_M": xm,
        "X_W_rel": relative_xw,
        "geometry": geometry,
        "rotation_w": rotation,
        "deformation_w": deformation,
    }


def potential_from_positions(positions, reference, hx, hy, targets, material,
                             nodal_load=None, include_xw=True):
    """Internal energy minus dead-load potential for a full position field."""
    value = paper_components(
        positions, reference, hx, hy, targets, material,
        include_xw=include_xw)["total"]
    if nodal_load is not None:
        value = value - (nodal_load * positions).sum()
    return value

