"""Intrinsic ``X_M`` completion of the Checkerboard membrane channel."""

from __future__ import annotations

import torch

from checkerboard_shell.geometry.core import torch_geometry


def _inverse_square_root(metric: torch.Tensor) -> torch.Tensor:
    determinant = torch.linalg.det(metric)
    if not bool(torch.all(determinant > 0.0)):
        raise ValueError("nondegenerate positively oriented B frames required")
    root_det = torch.sqrt(determinant)
    trace = torch.diagonal(metric, dim1=-2, dim2=-1).sum(-1)
    scale = torch.sqrt(trace + 2.0 * root_det)
    identity = torch.eye(2, dtype=metric.dtype, device=metric.device)
    return scale[..., None, None] * torch.linalg.inv(
        metric + root_det[..., None, None] * identity)


def polar_b_triads(vertices: torch.Tensor, hx: float, hy: float) -> torch.Tensor:
    geometry = torch_geometry(vertices, hx, hy)
    tangents = torch.stack((geometry.t1_b, geometry.t2_b), dim=-1)
    normalized = tangents @ _inverse_square_root(geometry.first_form)
    normal = torch.linalg.cross(
        normalized[..., :, 0], normalized[..., :, 1], dim=-1)
    normal = normal / torch.linalg.vector_norm(normal, dim=-1, keepdim=True)
    return torch.cat((normalized, normal[..., None]), dim=-1)


def relative_b_rotations(vertices, reference, hx, hy):
    current_frame = polar_b_triads(vertices, hx, hy)
    reference_frame = polar_b_triads(reference, hx, hy).detach()
    rotations = current_frame @ reference_frame.transpose(-1, -2)
    return rotations, reference_frame


def _axial_cayley(connection: torch.Tensor) -> torch.Tensor:
    identity = torch.eye(3, dtype=connection.dtype, device=connection.device)
    cayley = torch.linalg.solve(
        (connection + identity).transpose(-1, -2),
        (connection - identity).transpose(-1, -2)).transpose(-1, -2)
    return torch.stack((
        cayley[..., 2, 1] - cayley[..., 1, 2],
        cayley[..., 0, 2] - cayley[..., 2, 0],
        cayley[..., 1, 0] - cayley[..., 0, 1],
    ), dim=-1)


def _edge_normal(first, second):
    value = first + second
    return value / torch.linalg.vector_norm(value, dim=-1, keepdim=True)


def intrinsic_connection_field(vertices, reference, hx, hy):
    rotations, reference_frame = relative_b_rotations(
        vertices, reference, hx, hy)
    normals = reference_frame[..., :, 2]
    left, right = rotations[:, :-1], rotations[:, 1:]
    connection_x = left.transpose(-1, -2) @ right
    normal_x = _edge_normal(normals[:, :-1], normals[:, 1:])
    theta_x = (_axial_cayley(connection_x) * normal_x).sum(dim=-1)
    lower, upper = rotations[:-1], rotations[1:]
    connection_y = lower.transpose(-1, -2) @ upper
    normal_y = _edge_normal(normals[:-1], normals[1:])
    theta_y = (_axial_cayley(connection_y) * normal_y).sum(dim=-1)
    mixed_from_x = theta_x[1:] - theta_x[:-1]
    mixed_from_y = theta_y[:, 1:] - theta_y[:, :-1]
    mixed = 0.5 * (mixed_from_x + mixed_from_y)
    return {"rotation_b": rotations, "theta_x": theta_x, "theta_y": theta_y,
            "mixed_from_x": mixed_from_x, "mixed_from_y": mixed_from_y,
            "mixed": mixed}


def intrinsic_m_coordinate(vertices, reference, hx, hy, area_b=None):
    """Return the normalized M-character moment used by the paper energy."""
    field = intrinsic_connection_field(vertices, reference, hx, hy)["mixed"]
    ny, nx = field.shape
    jj, ii = torch.meshgrid(
        torch.arange(ny, device=field.device),
        torch.arange(nx, device=field.device), indexing="ij")
    parity = torch.where((ii + jj).remainder(2) == 0, 1.0, -1.0).to(field)
    if area_b is None:
        weights = torch.full_like(field, float(hx * hy))
    else:
        weights = 0.25 * (area_b[:-1, :-1] + area_b[:-1, 1:]
                          + area_b[1:, :-1] + area_b[1:, 1:])
    return (weights * parity * field).sum() / (4.0 * torch.sqrt(weights.sum()))

