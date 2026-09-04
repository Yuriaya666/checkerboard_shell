"""Moment-exact finite boundary closure for the frozen B/W geometry."""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache

import numpy as np
import torch
import torch.nn.functional as torch_functional

from checkerboard_shell.geometry.core import torch_geometry


@dataclass
class BoundaryGeometry:
    bulk: object
    t1_w: torch.Tensor
    t2_w: torch.Tensor
    normal_x: torch.Tensor
    normal_y: torch.Tensor
    d1_normal: torch.Tensor
    d2_normal: torch.Tensor
    first_form_w: torch.Tensor
    b_tilde: torch.Tensor
    second_form: torch.Tensor
    skew_diagnostic: torch.Tensor
    x_alternating: torch.Tensor
    frame_w: torch.Tensor


@dataclass
class BoundaryTargets:
    geometry: BoundaryGeometry
    first_form_b: torch.Tensor
    first_form_w: torch.Tensor
    second_form: torch.Tensor
    x_alternating: torch.Tensor
    frame_w: torch.Tensor
    area_b: torch.Tensor
    area_w: torch.Tensor
    gamma_w: torch.Tensor
    dual_fraction: torch.Tensor


@dataclass
class KoiterTargets:
    boundary: BoundaryTargets
    value_operator: torch.Tensor
    derivative_x_operator: torch.Tensor
    derivative_y_operator: torch.Tensor
    mixed_operator: torch.Tensor


def _moment_weights(nodes, target, derivative):
    """Polynomial moment weights exact through degree ``len(nodes)-1``."""
    nodes = np.asarray(nodes, dtype=float)
    degree = len(nodes) - 1
    vandermonde = np.stack([nodes**power for power in range(degree + 1)])
    rhs = np.zeros(degree + 1)
    for power in range(degree + 1):
        if derivative == 0:
            rhs[power] = target**power
        elif derivative == 1 and power:
            rhs[power] = power * target**(power - 1)
    return np.linalg.solve(vandermonde, rhs)


def _axis_stencil(index, count, spacing, route="moment"):
    """Cell-center samples used at raw vertex ``index``."""
    if count < 3:
        raise ValueError("boundary closure requires at least three raw cells")
    if index == 0:
        indices = np.array([0, 1, 2])
        nodes = np.array([0.5, 1.5, 2.5])
        target = 0.0
        direct_value = np.array([15.0/8.0, -5.0/4.0, 3.0/8.0])
        direct_derivative = np.array([-2.0, 3.0, -1.0]) / spacing
    elif index == count:
        indices = np.array([count - 3, count - 2, count - 1])
        nodes = np.array([count - 2.5, count - 1.5, count - 0.5])
        target = float(count)
        direct_value = np.array([3.0/8.0, -5.0/4.0, 15.0/8.0])
        direct_derivative = np.array([1.0, -3.0, 2.0]) / spacing
    else:
        indices = np.array([index - 1, index])
        nodes = np.array([index - 0.5, index + 0.5])
        target = float(index)
        direct_value = np.array([0.5, 0.5])
        direct_derivative = np.array([-1.0, 1.0]) / spacing
    if route == "direct":
        return indices, direct_value, direct_derivative
    if route != "moment":
        raise ValueError(f"unknown route {route}")
    value = _moment_weights(nodes, target, 0)
    derivative = _moment_weights(nodes, target, 1) / spacing
    return indices, value, derivative


def _tensor_reconstruct(field, j, i, hx, hy, route):
    """Return value, x/y derivatives and mixed derivative at a raw vertex."""
    ny, nx = field.shape[:2]
    ix, vx, dx = _axis_stencil(i, nx, hx, route)
    iy, vy, dy = _axis_stencil(j, ny, hy, route)
    value = field.new_zeros(field.shape[2:])
    grad_x = field.new_zeros(field.shape[2:])
    grad_y = field.new_zeros(field.shape[2:])
    mixed = field.new_zeros(field.shape[2:])
    for ay, by in enumerate(iy):
        for ax, bx in enumerate(ix):
            sample = field[int(by), int(bx)]
            value = value + float(vy[ay] * vx[ax]) * sample
            grad_x = grad_x + float(vy[ay] * dx[ax]) * sample
            grad_y = grad_y + float(dy[ay] * vx[ax]) * sample
            mixed = mixed + float(dy[ay] * dx[ax]) * sample
    return value, grad_x, grad_y, mixed


def _normalize(vector):
    return vector / torch.linalg.vector_norm(vector, dim=-1, keepdim=True)


def boundary_geometry(vertices, hx, hy, route="moment"):
    """Evaluate frozen bulk plus edge/corner W closure on every raw vertex."""
    bulk = torch_geometry(vertices, hx, hy)
    ny, nx = bulk.normals_b.shape[:2]
    t1_rows, t2_rows = [], []
    nx_rows, ny_rows, mixed_rows = [], [], []
    for j in range(ny + 1):
        row_t1, row_t2 = [], []
        row_nx, row_ny, row_mixed = [], [], []
        for i in range(nx + 1):
            normal_value, normal_x, normal_y, normal_mixed = (
                _tensor_reconstruct(bulk.normals_b, j, i, hx, hy, route))
            del normal_value
            if 0 < j < ny and 0 < i < nx:
                t1 = bulk.t1_w[j - 1, i - 1]
                t2 = bulk.t2_w[j - 1, i - 1]
            else:
                t1 = _tensor_reconstruct(bulk.t1_b, j, i, hx, hy, route)[0]
                t2 = _tensor_reconstruct(bulk.t2_b, j, i, hx, hy, route)[0]
            row_t1.append(t1)
            row_t2.append(t2)
            row_nx.append(normal_x)
            row_ny.append(normal_y)
            row_mixed.append(float(hx * hy) * normal_mixed)
        t1_rows.append(torch.stack(row_t1))
        t2_rows.append(torch.stack(row_t2))
        nx_rows.append(torch.stack(row_nx))
        ny_rows.append(torch.stack(row_ny))
        mixed_rows.append(torch.stack(row_mixed))
    t1_w = torch.stack(t1_rows)
    t2_w = torch.stack(t2_rows)
    normal_x = torch.stack(nx_rows)
    normal_y = torch.stack(ny_rows)
    x_alternating = torch.stack(mixed_rows)
    ell = float(np.hypot(hx, hy))
    e1x, e1y = float(hx / ell), float(hy / ell)
    e2x, e2y = float(-hx / ell), float(hy / ell)
    d1 = e1x * normal_x + e1y * normal_y
    d2 = e2x * normal_x + e2y * normal_y
    b_tilde = -torch.stack((
        torch.stack(((d1*t1_w).sum(-1), (d1*t2_w).sum(-1)), dim=-1),
        torch.stack(((d2*t1_w).sum(-1), (d2*t2_w).sum(-1)), dim=-1),
    ), dim=-2)
    second = 0.5 * (b_tilde + b_tilde.transpose(-1, -2))
    skew = 0.5 * (b_tilde - b_tilde.transpose(-1, -2))
    first = torch.stack((
        torch.stack(((t1_w*t1_w).sum(-1), (t1_w*t2_w).sum(-1)), dim=-1),
        torch.stack(((t2_w*t1_w).sum(-1), (t2_w*t2_w).sum(-1)), dim=-1),
    ), dim=-2)
    normal_w = _normalize(torch.linalg.cross(t1_w, t2_w, dim=-1))
    frame = torch.stack((t1_w, t2_w, normal_w), dim=-1)
    return BoundaryGeometry(
        bulk=bulk, t1_w=t1_w, t2_w=t2_w,
        normal_x=normal_x, normal_y=normal_y,
        d1_normal=d1, d2_normal=d2,
        first_form_w=first, b_tilde=b_tilde, second_form=second,
        skew_diagnostic=skew, x_alternating=x_alternating, frame_w=frame)


def checkerboard_to_raw_matrix(hx, hy, dtype, device=None):
    ell = float(np.hypot(hx, hy))
    return torch.tensor(
        [[hx / ell, -hx / ell], [hy / ell, hy / ell]],
        dtype=dtype, device=device)


def raw_metric_from_checkerboard(metric, hx, hy):
    frame = checkerboard_to_raw_matrix(hx, hy, metric.dtype, metric.device)
    inverse = torch.linalg.inv(frame)
    return inverse.transpose(-1, -2) @ metric @ inverse


def reference_area_and_gamma(metric, hx, hy):
    """Physical reference area and auxiliary-Q1 coefficient."""
    raw = raw_metric_from_checkerboard(metric, hx, hy)
    determinant = torch.linalg.det(raw)
    if bool(torch.any(determinant <= 0.0)):
        raise ValueError("reference tangent metric must be positive definite")
    sqrtg = torch.sqrt(determinant)
    inverse = torch.linalg.inv(raw)
    area = float(hx * hy) * sqrtg
    gamma = sqrtg * (
        inverse[..., 0, 0] * float(hy / hx)
        + inverse[..., 1, 1] * float(hx / hy)) / 12.0
    return area, gamma


def dual_fraction(ny, nx, dtype, device=None):
    wx = torch.ones(nx + 1, dtype=dtype, device=device)
    wy = torch.ones(ny + 1, dtype=dtype, device=device)
    wx[[0, -1]] = 0.5
    wy[[0, -1]] = 0.5
    return wy[:, None] * wx[None, :]


def targets_from_reference(reference, hx, hy, route="moment"):
    geometry = boundary_geometry(reference, hx, hy, route)
    area_b, _ = reference_area_and_gamma(geometry.bulk.first_form, hx, hy)
    area_w, gamma_w = reference_area_and_gamma(geometry.first_form_w, hx, hy)
    fractions = dual_fraction(
        reference.shape[0] - 1, reference.shape[1] - 1,
        reference.dtype, reference.device)
    return BoundaryTargets(
        geometry=geometry,
        first_form_b=geometry.bulk.first_form.detach().clone(),
        first_form_w=geometry.first_form_w.detach().clone(),
        second_form=geometry.second_form.detach().clone(),
        x_alternating=geometry.x_alternating.detach().clone(),
        frame_w=geometry.frame_w.detach().clone(),
        area_b=area_b.detach().clone(), area_w=area_w.detach().clone(),
        gamma_w=gamma_w.detach().clone(), dual_fraction=fractions)


def polar_rotation_newton(deformation, iterations=12):
    """Right polar factor on GL+(3), evaluated without SVD eigenvector AD."""
    rotation = deformation
    for _ in range(iterations):
        rotation = 0.5 * (
            rotation + torch.linalg.inv(rotation).transpose(-1, -2))
    return rotation


def transported_boundary_x(geometry, targets, iterations=12):
    deformation = geometry.frame_w @ torch.linalg.inv(targets.frame_w)
    rotation = polar_rotation_newton(deformation, iterations)
    transported = (rotation @ targets.x_alternating[..., None]).squeeze(-1)
    return geometry.x_alternating - transported, rotation, deformation


@lru_cache(maxsize=32)
def _operator_arrays(ny, nx, hx, hy):
    shape = ((ny + 1)*(nx + 1), ny*nx)
    value = np.zeros(shape)
    derivative_x = np.zeros(shape)
    derivative_y = np.zeros(shape)
    mixed = np.zeros(shape)
    for j in range(ny + 1):
        iy, vy, dy = _axis_stencil(j, ny, hy, "moment")
        for i in range(nx + 1):
            ix, vx, dx = _axis_stencil(i, nx, hx, "moment")
            row = j*(nx + 1) + i
            for ay, by in enumerate(iy):
                for ax, bx in enumerate(ix):
                    col = int(by)*nx + int(bx)
                    value[row, col] += vy[ay]*vx[ax]
                    derivative_x[row, col] += vy[ay]*dx[ax]
                    derivative_y[row, col] += dy[ay]*vx[ax]
                    mixed[row, col] += hx*hy*dy[ay]*dx[ax]
    return value, derivative_x, derivative_y, mixed


def fast_boundary_geometry(vertices, hx, hy, targets):
    """Algebraically identical vectorization of the validated closure."""
    bulk = torch_geometry(vertices, hx, hy)
    ny, nx = bulk.normals_b.shape[:2]
    matrices = (targets.value_operator, targets.derivative_x_operator,
                targets.derivative_y_operator, targets.mixed_operator)
    flat_n = bulk.normals_b.reshape(ny*nx, 3)
    value, dx, dy, mixed = (matrix @ flat_n for matrix in matrices)
    del value
    reconstructed_t1 = targets.value_operator @ bulk.t1_b.reshape(ny*nx, 3)
    reconstructed_t2 = targets.value_operator @ bulk.t2_b.reshape(ny*nx, 3)
    shape = (ny + 1, nx + 1, 3)
    reconstructed_t1 = reconstructed_t1.reshape(shape)
    reconstructed_t2 = reconstructed_t2.reshape(shape)
    padded_t1 = torch_functional.pad(
        bulk.t1_w.permute(2, 0, 1), (1, 1, 1, 1)).permute(1, 2, 0)
    padded_t2 = torch_functional.pad(
        bulk.t2_w.permute(2, 0, 1), (1, 1, 1, 1)).permute(1, 2, 0)
    interior = torch.zeros((ny + 1, nx + 1, 1), dtype=vertices.dtype,
                           device=vertices.device)
    interior[1:-1, 1:-1] = 1.0
    t1 = interior*padded_t1 + (1.0-interior)*reconstructed_t1
    t2 = interior*padded_t2 + (1.0-interior)*reconstructed_t2
    normal_x = dx.reshape(shape)
    normal_y = dy.reshape(shape)
    x_alternating = mixed.reshape(shape)
    ell = float(np.hypot(hx, hy))
    d1 = float(hx/ell)*normal_x + float(hy/ell)*normal_y
    d2 = float(-hx/ell)*normal_x + float(hy/ell)*normal_y
    b_tilde = -torch.stack((
        torch.stack(((d1*t1).sum(-1), (d1*t2).sum(-1)), -1),
        torch.stack(((d2*t1).sum(-1), (d2*t2).sum(-1)), -1)), -2)
    second = 0.5*(b_tilde + b_tilde.transpose(-1, -2))
    skew = 0.5*(b_tilde - b_tilde.transpose(-1, -2))
    first = torch.stack((
        torch.stack(((t1*t1).sum(-1), (t1*t2).sum(-1)), -1),
        torch.stack(((t2*t1).sum(-1), (t2*t2).sum(-1)), -1)), -2)
    normal = torch.linalg.cross(t1, t2, dim=-1)
    normal = normal/torch.linalg.vector_norm(normal, dim=-1, keepdim=True)
    frame = torch.stack((t1, t2, normal), -1)
    return BoundaryGeometry(
        bulk=bulk, t1_w=t1, t2_w=t2, normal_x=normal_x,
        normal_y=normal_y, d1_normal=d1, d2_normal=d2,
        first_form_w=first, b_tilde=b_tilde, second_form=second,
        skew_diagnostic=skew, x_alternating=x_alternating, frame_w=frame)


def koiter_targets(reference, hx, hy):
    boundary = targets_from_reference(reference, hx, hy, route="moment")
    ny, nx = reference.shape[0] - 1, reference.shape[1] - 1
    arrays = _operator_arrays(ny, nx, float(hx), float(hy))
    matrices = [torch.as_tensor(value, dtype=reference.dtype,
                                device=reference.device) for value in arrays]
    return KoiterTargets(boundary, *matrices)

