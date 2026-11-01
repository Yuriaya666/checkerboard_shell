"""Recovered manuscript benchmark definitions and executable drivers.

The definitions in this file were recovered from the clean historical Git
repository ``/home/song/projects/CB-shell`` at tag ``v1.0-arxiv``.  They are
ported to the current package without changing geometry, loading, boundary
conditions, quadrature, observables, or the accepted energy.
"""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import scipy.linalg
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch

from checkerboard_shell.geometry.boundary import _axis_stencil
from checkerboard_shell.material.isotropic import IsotropicMaterial, lame_surface_coefficients
from checkerboard_shell.mechanics.energy import make_targets
from checkerboard_shell.mechanics.low_rank import (
    LowRankTangent,
    assemble_paper_low_rank_tangent,
)


PLATE_REFERENCE = 4.4360891e-4
ROOF_REFERENCE = 0.3024
ROOF_MODERN_REFERENCE = 0.3005924566
M3_REFERENCE = 0.09352155


def _cell_normal_operators(n: int, h: float) -> tuple[sp.csr_matrix, sp.csr_matrix]:
    rows_x, cols_x, data_x = [], [], []
    rows_y, cols_y, data_y = [], [], []
    vertex = lambda i, j: j * (n + 1) + i
    cell = lambda i, j: j * n + i
    for j in range(n):
        for i in range(n):
            row = cell(i, j)
            q00, q10 = vertex(i, j), vertex(i + 1, j)
            q01, q11 = vertex(i, j + 1), vertex(i + 1, j + 1)
            for col, value in ((q00, 0.5/h), (q01, 0.5/h),
                               (q10, -0.5/h), (q11, -0.5/h)):
                rows_x.append(row); cols_x.append(col); data_x.append(value)
            for col, value in ((q00, 0.5/h), (q10, 0.5/h),
                               (q01, -0.5/h), (q11, -0.5/h)):
                rows_y.append(row); cols_y.append(col); data_y.append(value)
    shape = (n*n, (n + 1)*(n + 1))
    return (
        sp.coo_matrix((data_x, (rows_x, cols_x)), shape=shape).tocsr(),
        sp.coo_matrix((data_y, (rows_y, cols_y)), shape=shape).tocsr(),
    )


def _plate_reconstruction_operators(n: int, h: float):
    entries = [([], [], []) for _ in range(3)]
    site = lambda i, j: j * (n + 1) + i
    cell = lambda i, j: j * n + i
    for j in range(n + 1):
        iy, vy, dy = _axis_stencil(j, n, h, "moment")
        for i in range(n + 1):
            ix, vx, dx = _axis_stencil(i, n, h, "moment")
            row = site(i, j)
            for ay, by in enumerate(iy):
                for ax, bx in enumerate(ix):
                    col = cell(int(bx), int(by))
                    for target, value in zip(
                        entries,
                        (vy[ay] * dx[ax], dy[ay] * vx[ax],
                         h*h * dy[ay] * dx[ax]),
                    ):
                        target[0].append(row); target[1].append(col)
                        target[2].append(float(value))
    shape = ((n + 1)*(n + 1), n*n)
    return tuple(sp.coo_matrix((data, (rows, cols)), shape=shape).tocsr()
                 for rows, cols, data in entries)


def plate_stiffness(n: int, *, include_xw: bool = True) -> sp.csr_matrix:
    """Recovered manuscript transverse Hessian for the simply supported plate."""
    h = 1.0 / n
    gx, gy = _cell_normal_operators(n, h)
    rx, ry, mixed = _plate_reconstruction_operators(n, h)
    nx_x, nx_y = rx @ gx, ry @ gx
    ny_x, ny_y = rx @ gy, ry @ gy
    root2 = math.sqrt(2.0)
    d1x, d1y = (nx_x + nx_y)/root2, (ny_x + ny_y)/root2
    d2x, d2y = (-nx_x + nx_y)/root2, (-ny_x + ny_y)/root2
    b11 = -(d1x + d1y)/root2
    b12 = -0.5*((-d1x + d1y + d2x + d2y)/root2)
    b22 = -(-d2x + d2y)/root2
    xx, xy = mixed @ gx, mixed @ gy
    alpha, mu = lame_surface_coefficients(1.0, 0.3)
    edge = np.ones(n + 1); edge[[0, -1]] = 0.5
    area = sp.diags(np.outer(edge, edge).reshape(-1) * h*h)
    d = 0.5*alpha + mu
    coefficient = 0.01**3 / 12.0
    stiffness = coefficient * (
        2.0*d*(b11.T @ area @ b11)
        + 4.0*mu*(b12.T @ area @ b12)
        + 2.0*d*(b22.T @ area @ b22)
        + alpha*(b11.T @ area @ b22 + b22.T @ area @ b11))
    if include_xw:
        interior = np.ones((n + 1, n + 1))
        interior[[0, -1], :] = 0.0; interior[:, [0, -1]] = 0.0
        xweight = sp.diags(interior.reshape(-1) / 6.0)
        stiffness += 2.0*coefficient*mu*(
            xx.T @ xweight @ xx + xy.T @ xweight @ xy)
    return stiffness.tocsr()


def solve_plate(n: int, *, variant: str = "CB-Full") -> dict:
    if variant not in ("CB-0", "CB-M", "CB-Full"):
        raise ValueError(variant)
    stiffness = plate_stiffness(n, include_xw=(variant == "CB-Full"))
    grid = np.arange((n + 1)**2).reshape(n + 1, n + 1)
    active = grid[1:-1, 1:-1].reshape(-1)
    reduced = stiffness[active][:, active].tocsc()
    force = np.full(active.size, 1.0e-8/(n*n))
    solution = spla.spsolve(reduced, force)
    center_raw = grid[n//2, n//2]
    center = float(solution[np.flatnonzero(active == center_raw)[0]])
    return {
        "benchmark": "simply_supported_plate",
        "n": n, "variant": variant,
        "displacement": center,
        "reference": PLATE_REFERENCE,
        "relative_error": abs(center-PLATE_REFERENCE)/PLATE_REFERENCE,
        "raw_scalar_states": 3*(n + 1)**2,
        "active_states": int(active.size),
        "local_tangent_nnz": int(reduced.nnz),
        "observable": "raw transverse center vertex",
    }


def roof_problem(n: int) -> dict:
    radius, half_length = 25.0, 25.0
    half_angle = math.radians(40.0)
    x = torch.linspace(0.0, half_length, n + 1, dtype=torch.float64)
    arc = torch.linspace(-radius*half_angle, radius*half_angle,
                         n + 1, dtype=torch.float64)
    ss, xx = torch.meshgrid(arc, x, indexing="ij")
    theta = ss/radius
    reference = torch.stack((xx, radius*torch.sin(theta),
                             radius*torch.cos(theta)), -1)
    hx, hy = half_length/n, 2.0*radius*half_angle/n
    fixed = set()
    columns = n + 1
    for j in range(n + 1):
        fixed.add(3*(j*columns) + 0)
        fixed.add(3*(j*columns + n) + 1)
        fixed.add(3*(j*columns + n) + 2)
    active = np.asarray([k for k in range(reference.numel()) if k not in fixed])
    wx = torch.ones(n + 1, dtype=reference.dtype)
    wy = torch.ones(n + 1, dtype=reference.dtype)
    wx[[0, -1]] = 0.5; wy[[0, -1]] = 0.5
    weights = hx*hy * wy[:, None] * wx[None, :]
    load = torch.zeros_like(reference); load[..., 2] = -90.0*weights
    output_raw = 3*(n*columns + 0) + 2
    return {
        "reference": reference, "hx": hx, "hy": hy,
        "targets": make_targets(reference, hx, hy),
        "material": IsotropicMaterial(4.32e8, 0.0, 0.25),
        "active": active, "force": load.reshape(-1).numpy()[active],
        "output_raw": output_raw,
        "output_active": int(np.flatnonzero(active == output_raw)[0]),
    }


def _vertical_midpoint_operator(n: int, i: int, j: int) -> np.ndarray:
    operator = np.zeros((3, 3*(n + 1)**2))
    for row in (j, j + 1):
        vertex = row*(n + 1) + i
        for component in range(3):
            operator[component, 3*vertex + component] += 0.5
    return operator


def _equator_point_operator(n: int, i: int) -> np.ndarray:
    coefficients = (15.0/8.0, -5.0/4.0, 3.0/8.0)
    return sum(coefficient*_vertical_midpoint_operator(n, i, n-1-offset)
               for offset, coefficient in enumerate(coefficients))


def m3_problem(n: int) -> dict:
    theta0 = math.radians(18.0)
    theta = torch.linspace(theta0, math.pi/2.0, n + 1, dtype=torch.float64)
    phi = torch.linspace(0.0, math.pi/2.0, n + 1, dtype=torch.float64)
    tt, pp = torch.meshgrid(theta, phi, indexing="ij")
    reference = 10.0*torch.stack((
        torch.sin(tt)*torch.cos(pp), torch.sin(tt)*torch.sin(pp), torch.cos(tt)), -1)
    hx, hy = (math.pi/2.0)/n, (math.pi/2.0-theta0)/n
    count = reference.numel()
    constraints = []
    for j in range(n):
        left = [_vertical_midpoint_operator(n, i, j) for i in (0, 1, 2)]
        right = [_vertical_midpoint_operator(n, i, j) for i in (n, n-1, n-2)]
        constraints.append(left[0][1]); constraints.append(right[0][0])
        for component in (0, 2):
            constraints.append((-3.0*left[0] + 4.0*left[1] - left[2])[component])
        for component in (1, 2):
            constraints.append((3.0*right[0] - 4.0*right[1] + right[2])[component])
    parity = np.fromfunction(
        lambda j, i: (-1.0)**(i+j), (n+1, n+1), dtype=int).reshape(-1)
    for component in range(3):
        row = np.zeros(count); row[component::3] = parity/parity.size
        constraints.append(row)
    mean_z = np.zeros(count)
    for j in range(n):
        for i in range(n + 1):
            mean_z += _vertical_midpoint_operator(n, i, j)[2]
    constraints.append(mean_z/(n*(n + 1)))
    constraint = sp.csr_matrix(np.stack(constraints))
    point_a, point_b = _equator_point_operator(n, 0), _equator_point_operator(n, n)
    force = point_a[0] - point_b[1]
    return {
        "reference": reference, "hx": hx, "hy": hy,
        "targets": make_targets(reference, hx, hy),
        "material": IsotropicMaterial(6.825e7, 0.3, 0.04),
        "constraints": constraint, "point_a": point_a, "point_b": point_b,
        "force": force,
    }


def solve_curved_sparse(problem: dict, *, variant: str = "CB-Full") -> tuple[np.ndarray, dict]:
    reference = problem["reference"]
    tangent, residual, audit = assemble_paper_low_rank_tangent(
        torch.zeros_like(reference), reference, problem["hx"], problem["hy"],
        problem["targets"], problem["material"], variant=variant,
    )
    if "active" in problem:
        active = problem["active"]
        reduced = LowRankTangent(
            tangent.local[active][:, active].tocsr(),
            tangent.vector[active], tangent.coefficient,
        )
        solution, solve_audit = reduced.solve(problem["force"])
        full = np.zeros(reference.numel()); full[active] = solution
    else:
        full, solve_audit = tangent.solve(
            problem["force"], constraints=problem["constraints"])
    audit.update(solve_audit)
    audit["reference_residual_inf"] = float(np.max(np.abs(residual)))
    return full, audit


def solve_roof(n: int, *, variant: str = "CB-Full") -> dict:
    problem = roof_problem(n)
    full, audit = solve_curved_sparse(problem, variant=variant)
    value = -float(full[problem["output_raw"]])
    return {
        "benchmark": "Scordelis_Lo_roof", "n": n, "variant": variant,
        "displacement": value, "reference": ROOF_REFERENCE,
        "modern_reference": ROOF_MODERN_REFERENCE,
        "relative_error": abs(value-ROOF_REFERENCE)/ROOF_REFERENCE,
        "raw_scalar_states": int(problem["reference"].numel()),
        "active_states": int(problem["active"].size),
        **audit,
    }


def solve_m3(n: int, *, variant: str = "CB-Full") -> dict:
    problem = m3_problem(n)
    full, audit = solve_curved_sparse(problem, variant=variant)
    midpoint = -float(problem["point_b"][1] @ full)
    # The raw equator vertex (i=n, j=n) is the same continuum point B, but
    # it is the discrete observable used by the LibShell MidedgeTan driver.
    raw_index = 3*(n*(n + 1) + n) + 1
    raw = -float(full[raw_index])
    return {
        "benchmark": "MacNeal_Harder_M3", "n": n, "variant": variant,
        "displacement": midpoint,
        "displacement_midpoint": midpoint,
        "displacement_raw": raw,
        "reference": M3_REFERENCE,
        "relative_error": abs(midpoint-M3_REFERENCE)/M3_REFERENCE,
        "relative_error_midpoint": abs(midpoint-M3_REFERENCE)/M3_REFERENCE,
        "relative_error_raw": abs(raw-M3_REFERENCE)/M3_REFERENCE,
        "observable_midpoint": "physical-midpoint moment reconstruction at B",
        "observable_raw": "raw equator vertex (i=n,j=n) y displacement",
        "raw_scalar_states": int(problem["reference"].numel()),
        "active_states": int(problem["reference"].numel()-problem["constraints"].shape[0]),
        "constraint_count": int(problem["constraints"].shape[0]),
        **audit,
    }
