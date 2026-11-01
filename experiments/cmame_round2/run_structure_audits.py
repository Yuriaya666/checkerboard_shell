"""Derivative-locality and boundary-storage audits for Round 2."""

from __future__ import annotations

import json
from pathlib import Path
import time

import numpy as np
import torch

from checkerboard_shell.geometry.boundary import fast_boundary_geometry, koiter_targets
from checkerboard_shell.geometry.intrinsic import (
    intrinsic_connection_field,
    intrinsic_m_coordinate,
)
from checkerboard_shell.mechanics.energy import paper_components
from checkerboard_shell.mechanics.low_rank import xm_energy_coefficient
from experiments.cmame_round2.recovered_drivers import m3_problem


ROOT = Path(__file__).resolve().parent


def derivative_audit() -> dict:
    problem = m3_problem(6)
    reference = problem["reference"]
    phase = torch.arange(reference.numel(), dtype=reference.dtype).reshape_as(reference)
    displacement = 2.5e-3*torch.sin(0.23*phase)
    flat = displacement.reshape(-1).detach().requires_grad_(True)

    def xm(values):
        current = reference + values.reshape_as(reference)
        return intrinsic_m_coordinate(
            current, reference, problem["hx"], problem["hy"],
            problem["targets"].checkerboard.boundary.area_b)

    x = xm(flat)
    g, = torch.autograd.grad(x, flat, create_graph=True)
    hx = torch.autograd.functional.hessian(xm, flat, vectorize=True).detach().numpy()
    outer = np.outer(g.detach().numpy(), g.detach().numpy())
    tolerance = max(np.max(np.abs(hx)), 1.0)*1.0e-12
    nz = np.argwhere(np.abs(hx) > tolerance)
    nx = reference.shape[1]-1
    offsets = []
    for row, col in nz:
        node_r, node_c = int(row)//3, int(col)//3
        jr, ir = divmod(node_r, nx+1); jc, ic = divmod(node_c, nx+1)
        offsets.append((abs(jr-jc), abs(ir-ic)))

    field = intrinsic_connection_field(
        reference + displacement, reference, problem["hx"], problem["hy"])["mixed"]
    area = problem["targets"].checkerboard.boundary.area_b
    weights = 0.25*(area[:-1, :-1]+area[:-1, 1:]+area[1:, :-1]+area[1:, 1:])
    jj, ii = torch.meshgrid(torch.arange(field.shape[0]), torch.arange(field.shape[1]),
                            indexing="ij")
    parity = torch.where((ii+jj).remainder(2)==0, 1.0, -1.0).to(field)
    local_sum = (weights*parity*field/(4.0*torch.sqrt(weights.sum()))).sum()

    coefficient = xm_energy_coefficient(problem["material"])
    eigen_outer = np.linalg.eigvalsh(outer)
    payload = {
        "case": "finite deformed curved M3 n=6",
        "source_energy": "E_XM=(t/4)*(8*mu)*X_M^2=c*X_M^2",
        "t": problem["material"].thickness,
        "mu": problem["material"].young/(2*(1+problem["material"].poisson)),
        "c": coefficient,
        "rank_one_coefficient_2c": 2*coefficient,
        "X_M": float(x.detach()),
        "local_plaquette_sum_difference": abs(float(x.detach()-local_sum.detach())),
        "g_size": int(g.numel()),
        "g_nonzero_count": int(np.count_nonzero(np.abs(g.detach().numpy()) > 1e-14)),
        "outer_rank": int(np.linalg.matrix_rank(outer, tol=max(abs(eigen_outer[-1]),1.0)*1e-12)),
        "outer_nnz": int(np.count_nonzero(np.abs(outer) > 1e-24)),
        "outer_density": float(np.count_nonzero(np.abs(outer) > 1e-24)/outer.size),
        "H_X_nnz": int(nz.shape[0]),
        "H_X_density": float(nz.shape[0]/hx.size),
        "H_X_threshold": tolerance,
        "H_X_max_node_offset_j": max(pair[0] for pair in offsets),
        "H_X_max_node_offset_i": max(pair[1] for pair in offsets),
        "consistent_tangent_identity": "2*c*(g*g^T + X_M*H_X)",
        "local_term_retained": "2*c*X_M*H_X",
    }
    out = ROOT/"02_low_rank_xm"/"raw_data"
    out.mkdir(parents=True, exist_ok=True)
    (out/"xm_derivative_structure.json").write_text(
        json.dumps(payload, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return payload


def boundary_audit() -> dict:
    rows = []
    for n in (8, 16, 32, 64, 128):
        theta0 = np.deg2rad(18.0)
        theta = torch.linspace(theta0, np.pi/2, n+1, dtype=torch.float64)
        phi = torch.linspace(0.0, np.pi/2, n+1, dtype=torch.float64)
        tt, pp = torch.meshgrid(theta, phi, indexing="ij")
        reference = 10.0*torch.stack((
            torch.sin(tt)*torch.cos(pp), torch.sin(tt)*torch.sin(pp), torch.cos(tt)), -1)
        hx, hy = (np.pi/2)/n, (np.pi/2-theta0)/n
        targets = koiter_targets(reference, hx, hy)
        stencil = targets.reconstruction
        sites, faces = (n+1)**2, n*n
        dense_bytes = 4*sites*faces*8
        samples = []
        for _ in range(4):
            start = time.perf_counter()
            geometry = fast_boundary_geometry(
                reference, hx, hy, targets)
            float(geometry.second_form.sum())
            samples.append(time.perf_counter()-start)
        rows.append({
            "n": n, "sites": sites, "B_faces": faces,
            "old_dense_bytes_projected": dense_bytes,
            "old_dense_status": "projected_from_exact_shape_not_allocated",
            "new_stencil_bytes_measured": stencil.storage_bytes,
            "new_over_old_storage_ratio": stencil.storage_bytes/dense_bytes,
            "new_geometry_time_median_seconds": float(np.median(samples[1:])),
            "stencil_width": int(stencil.indices.shape[1]),
        })
    out = ROOT/"01_sparse_boundary"/"raw_data"
    out.mkdir(parents=True, exist_ok=True)
    (out/"boundary_storage.json").write_text(
        json.dumps(rows, indent=2)+"\n", encoding="utf-8")
    return {"rows": rows}


if __name__ == "__main__":
    print(json.dumps({"derivative": derivative_audit(),
                      "boundary": boundary_audit()}, indent=2))
