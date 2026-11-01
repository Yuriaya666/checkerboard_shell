"""Exact sparse/local plus rank-one tangent for the Checkerboard shell.

The accepted membrane completion has the scalar form ``c * X_M(q)**2``.
Its consistent tangent is split exactly as

``K_XM = 2*c*X_M*D2X_M + 2*c*DX_M*DX_M.T``.

The first term is local because ``X_M`` is a weighted sum of local
plaquette observations.  The second term is the only global contribution and
is retained as one vector plus one scalar; it is never explicitly assembled
as a dense outer product in this module.
"""

from __future__ import annotations

from dataclasses import dataclass
from functools import lru_cache
import time
from typing import Callable

import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import torch

from checkerboard_shell.material.isotropic import lame_surface_coefficients
from checkerboard_shell.mechanics.energy import paper_components


TensorFunction = Callable[[torch.Tensor], torch.Tensor]


def xm_energy_coefficient(material) -> float:
    """Return ``c`` from the source-level energy ``E_XM = c X_M^2``."""
    _, mu = lame_surface_coefficients(material.young, material.poisson)
    return (float(material.thickness) / 4.0) * (8.0 * float(mu))


@dataclass(frozen=True)
class LowRankTangent:
    """A local sparse matrix plus one exact global rank-one correction."""

    local: sp.csr_matrix
    vector: np.ndarray
    coefficient: float

    @property
    def shape(self) -> tuple[int, int]:
        return self.local.shape

    @property
    def rank_one_storage(self) -> int:
        return int(self.vector.size + 1)

    def action(self, value: np.ndarray) -> np.ndarray:
        vector = np.asarray(value, dtype=float)
        return self.local @ vector + self.coefficient * self.vector * (
            self.vector @ vector)

    def augmented_matrix(self, constraints: sp.spmatrix | None = None) -> sp.csc_matrix:
        """Return the exact bordered system, optionally with linear constraints."""
        if self.coefficient == 0.0:
            if constraints is None:
                return self.local.tocsc()
            cmat = sp.csr_matrix(constraints)
            return sp.bmat(
                [[self.local, cmat.T],
                 [cmat, sp.csc_matrix((cmat.shape[0], cmat.shape[0]))]],
                format="csc",
            )
        g = sp.csc_matrix(self.vector.reshape(-1, 1))
        inverse = -1.0 / float(self.coefficient)
        if constraints is None:
            return sp.bmat(
                [[self.local, g], [g.T, sp.csc_matrix([[inverse]])]],
                format="csc",
            )
        cmat = sp.csr_matrix(constraints)
        zero_rank_constraint = sp.csc_matrix((1, cmat.shape[0]))
        zero_constraints = sp.csc_matrix((cmat.shape[0], cmat.shape[0]))
        return sp.bmat(
            [
                [self.local, g, cmat.T],
                [g.T, sp.csc_matrix([[inverse]]), zero_rank_constraint],
                [cmat, zero_rank_constraint.T, zero_constraints],
            ],
            format="csc",
        )

    def solve(
        self,
        rhs: np.ndarray,
        constraints: sp.spmatrix | None = None,
        constraint_rhs: np.ndarray | None = None,
    ) -> tuple[np.ndarray, dict[str, float]]:
        matrix = self.augmented_matrix(constraints)
        right = np.asarray(rhs, dtype=float)
        extra = ([] if self.coefficient == 0.0 else [0.0])
        if constraints is not None:
            count = constraints.shape[0]
            values = (np.zeros(count) if constraint_rhs is None
                      else np.asarray(constraint_rhs, dtype=float))
            extra.extend(values.tolist())
        bordered_rhs = np.concatenate((right, np.asarray(extra)))
        started = time.perf_counter()
        factor = spla.splu(matrix)
        factor_seconds = time.perf_counter() - started
        started = time.perf_counter()
        solution = factor.solve(bordered_rhs)
        solve_seconds = time.perf_counter() - started
        return solution[: self.local.shape[0]], {
            "factorization_seconds": factor_seconds,
            "linear_solve_seconds": solve_seconds,
            "augmented_nnz": int(matrix.nnz),
        }


def _grid_colours(ny: int, nx: int, radius: int) -> np.ndarray:
    period = 2 * int(radius) + 1
    colours = np.empty(3 * (ny + 1) * (nx + 1), dtype=np.int64)
    for j in range(ny + 1):
        for i in range(nx + 1):
            node = j * (nx + 1) + i
            for component in range(3):
                colours[3 * node + component] = (
                    component * period * period
                    + (j % period) * period + (i % period)
                )
    return colours


def _neighbour_rows(column: int, ny: int, nx: int, radius: int) -> np.ndarray:
    node, _ = divmod(int(column), 3)
    j, i = divmod(node, nx + 1)
    rows = []
    for jj in range(max(0, j-radius), min(ny, j+radius) + 1):
        for ii in range(max(0, i-radius), min(nx, i+radius) + 1):
            base = 3 * (jj * (nx + 1) + ii)
            rows.extend((base, base + 1, base + 2))
    return np.asarray(rows, dtype=np.int64)


@lru_cache(maxsize=12)
def _cached_colour_batches(ny: int, nx: int, radius: int, batch_size: int):
    """Topology-only seed and COO maps; no configuration-dependent values."""
    colours = _grid_colours(ny, nx, radius)
    unique = np.unique(colours)
    count = colours.size
    batches = []
    for first in range(0, unique.size, batch_size):
        colour_batch = unique[first:first + batch_size]
        seeds = np.zeros((len(colour_batch), count), dtype=np.float64)
        rows, columns, batch_rows = [], [], []
        for local_index, colour in enumerate(colour_batch):
            for column in np.flatnonzero(colours == colour):
                seeds[local_index, column] = 1.0
                local_rows = _neighbour_rows(int(column), ny, nx, radius)
                rows.append(local_rows)
                columns.append(np.full(local_rows.size, column, dtype=np.int64))
                batch_rows.append(np.full(local_rows.size, local_index, dtype=np.int64))
        batches.append((seeds, np.concatenate(rows), np.concatenate(columns),
                        np.concatenate(batch_rows)))
    return tuple(batches), int(unique.size)


def assemble_coloured_local_hessian(
    local_proxy: TensorFunction,
    coordinates: torch.Tensor,
    ny: int,
    nx: int,
    *,
    radius: int = 3,
    colour_batch: int = 21,
    numerical_drop_tolerance: float = 0.0,
    assembly_strategy: str = "cached",
) -> tuple[sp.csr_matrix, dict[str, float]]:
    """Assemble an AD Hessian using a fixed local dependency colouring.

    ``radius=3`` is the conservative raw-node radius of the accepted
    boundary moment closure (three B samples imply at most four raw nodes).
    Same-colour columns have disjoint possible row supports, so each batched
    Hessian-vector product recovers multiple columns without approximation.
    """
    expected = 3 * (ny + 1) * (nx + 1)
    if coordinates.ndim != 1 or coordinates.numel() != expected:
        raise ValueError("coordinates must contain all raw scalar displacements")
    variable = coordinates.detach().clone().requires_grad_(True)
    started = time.perf_counter()
    value = local_proxy(variable)
    gradient, = torch.autograd.grad(value, variable, create_graph=True)
    graph_seconds = time.perf_counter() - started
    row_parts: list[np.ndarray] = []
    col_parts: list[np.ndarray] = []
    data_parts: list[np.ndarray] = []
    hvp_started = time.perf_counter()
    if assembly_strategy == "cached":
        batches, colour_count = _cached_colour_batches(
            ny, nx, radius, int(colour_batch))
        for index, (seed_np, rows, columns, batch_rows) in enumerate(batches):
            seeds = torch.as_tensor(seed_np, dtype=variable.dtype,
                                    device=variable.device)
            hvp, = torch.autograd.grad(
                gradient, variable, grad_outputs=seeds,
                is_grads_batched=True, retain_graph=(index+1 < len(batches)))
            values = hvp.detach().cpu().numpy()[batch_rows, rows]
            if numerical_drop_tolerance > 0.0:
                keep = np.abs(values) > numerical_drop_tolerance
                rows, columns, values = rows[keep], columns[keep], values[keep]
            row_parts.append(rows)
            col_parts.append(columns)
            data_parts.append(values)
    elif assembly_strategy == "legacy":
        colours = _grid_colours(ny, nx, radius)
        unique = np.unique(colours)
        colour_count = int(unique.size)
        for first in range(0, unique.size, int(colour_batch)):
            batch = unique[first:first + int(colour_batch)]
            seeds = torch.zeros(
                (batch.size, expected), dtype=variable.dtype, device=variable.device)
            for local_index, colour in enumerate(batch):
                columns = np.flatnonzero(colours == colour)
                seeds[local_index, torch.as_tensor(columns, device=variable.device)] = 1.0
            hvp, = torch.autograd.grad(
                gradient, variable, grad_outputs=seeds,
                is_grads_batched=True, retain_graph=(first + batch.size < unique.size),
            )
            hvp_np = hvp.detach().cpu().numpy()
            for local_index, colour in enumerate(batch):
                for column in np.flatnonzero(colours == colour):
                    rows = _neighbour_rows(int(column), ny, nx, radius)
                    values = hvp_np[local_index, rows]
                    if numerical_drop_tolerance > 0.0:
                        keep = np.abs(values) > numerical_drop_tolerance
                        rows, values = rows[keep], values[keep]
                    row_parts.append(rows)
                    col_parts.append(np.full(rows.size, int(column), dtype=np.int64))
                    data_parts.append(values)
    else:
        raise ValueError(f"unknown Hessian assembly strategy {assembly_strategy!r}")
    hvp_seconds = time.perf_counter() - hvp_started
    matrix = sp.coo_matrix(
        (np.concatenate(data_parts), (np.concatenate(row_parts),
                                      np.concatenate(col_parts))),
        shape=(expected, expected),
    ).tocsr()
    # Both triangles were independently recovered.  Averaging removes only AD
    # roundoff asymmetry and preserves the exact consistent Hessian.
    matrix = (0.5 * (matrix + matrix.T)).tocsr()
    matrix.eliminate_zeros()
    return matrix, {
        "graph_seconds": graph_seconds,
        "hvp_seconds": hvp_seconds,
        "colours": colour_count,
        "assembly_strategy": assembly_strategy,
        "radius": int(radius),
        "structural_nnz": int(matrix.nnz),
    }


def assemble_paper_low_rank_tangent(
    displacement: torch.Tensor,
    reference: torch.Tensor,
    hx: float,
    hy: float,
    targets,
    material,
    *,
    variant: str = "CB-Full",
    xm_multiplier: float = 1.0,
    radius: int = 3,
    colour_batch: int = 21,
    assembly_strategy: str = "cached",
) -> tuple[LowRankTangent, np.ndarray, dict[str, float]]:
    """Assemble the exact nonlinear tangent without a dense outer product."""
    if variant not in ("CB-0", "CB-M", "CB-W", "CB-Full"):
        raise ValueError(f"unknown formulation variant {variant!r}")
    if float(xm_multiplier) < 0.0:
        raise ValueError("xm_multiplier must be nonnegative")
    ny, nx = reference.shape[0] - 1, reference.shape[1] - 1
    zero = displacement.detach().reshape(-1)
    coefficient = xm_energy_coefficient(material)
    alpha = 2.0 * coefficient * float(xm_multiplier)

    def parts(values: torch.Tensor):
        positions = reference + values.reshape_as(reference)
        return paper_components(
            positions, reference, hx, hy, targets, material,
            include_xw=(variant in ("CB-W", "CB-Full")),
        )

    evaluation = parts(zero)
    x_value = evaluation["X_M"].detach()
    at_reference = bool(torch.count_nonzero(zero).item() == 0)
    if at_reference:
        # X_M(reference; reference) is identically zero by construction.  Use
        # that exact value instead of propagating a ~1e-34 roundoff remainder
        # through the otherwise unnecessary X_M*D2X_M graph.
        x_value = zero.new_zeros(())

    def ordinary(values: torch.Tensor) -> torch.Tensor:
        values_parts = parts(values)
        result = values_parts["membrane_b"] + values_parts["bending"]
        if variant in ("CB-W", "CB-Full"):
            result = result + values_parts["bending_compatibility"]
        return result

    def full_energy(values: torch.Tensor) -> torch.Tensor:
        value = ordinary(values)
        if variant in ("CB-M", "CB-Full"):
            value = (
                value + float(xm_multiplier) * coefficient
                * parts(values)["X_M"].square())
        return value

    variable = zero.detach().clone().requires_grad_(True)
    residual_started = time.perf_counter()
    total = full_energy(variable)
    residual, = torch.autograd.grad(total, variable)
    residual_seconds = time.perf_counter() - residual_started
    if variant in ("CB-0", "CB-W") or alpha == 0.0:
        g = np.zeros(variable.numel())
        rank_coefficient = 0.0
        rank_vector_seconds = 0.0

        def proxy(values: torch.Tensor) -> torch.Tensor:
            return ordinary(values)
    else:
        rank_started = time.perf_counter()
        x = parts(variable)["X_M"]
        gradient_x, = torch.autograd.grad(x, variable)
        g = gradient_x.detach().cpu().numpy()
        rank_vector_seconds = time.perf_counter() - rank_started
        rank_coefficient = alpha

        if at_reference:
            def proxy(values: torch.Tensor) -> torch.Tensor:
                return ordinary(values)
        else:
            def proxy(values: torch.Tensor) -> torch.Tensor:
                # Hessian = K_ordinary + 2*c*X_M(current)*D2X_M.
                return ordinary(values) + alpha * x_value * parts(values)["X_M"]

    local, audit = assemble_coloured_local_hessian(
        proxy, zero, ny, nx, radius=radius, colour_batch=colour_batch,
        assembly_strategy=assembly_strategy)
    audit.update({
        "X_M": float(x_value),
        "xm_energy_coefficient_c": coefficient,
        "xm_multiplier": float(xm_multiplier),
        "rank_one_coefficient_2c": (
            alpha if variant in ("CB-M", "CB-Full") and alpha != 0.0 else 0.0),
        "residual_seconds": residual_seconds,
        "rank_one_vector_seconds": rank_vector_seconds,
    })
    return (
        LowRankTangent(local, g, rank_coefficient),
        residual.detach().cpu().numpy(),
        audit,
    )
