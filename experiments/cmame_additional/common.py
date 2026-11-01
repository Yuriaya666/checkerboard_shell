"""Shared, reproducible utilities for the CMAME additional experiments.

This module intentionally lives outside :mod:`checkerboard_shell`.  It does not
change the accepted paper formulation; formulation switches are assembled from
the already exposed energy components.
"""

from __future__ import annotations

from dataclasses import dataclass
import csv
import json
import math
import os
from pathlib import Path
import platform
import resource
import subprocess
import time
from typing import Callable, Iterable

import numpy as np
import torch

from checkerboard_shell.geometry.core import edge_midpoints
from checkerboard_shell.geometry.mesh import CheckerboardMesh
from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.solver.boundary_conditions import BoundaryConditions


ROOT = Path(__file__).resolve().parent
SEED = 20260921
DTYPE = torch.float64

BENCHMARKS = {
    "plate": {
        "E": 1.0, "nu": 0.3, "thickness": 0.01,
        "pressure_z": -1.0e-8, "reference": 4.4360891e-4,
        "observable": "downward center displacement",
    },
    "roof": {
        "R": 25.0, "L": 50.0, "half_angle_degrees": 40.0,
        "E": 4.32e8, "nu": 0.0, "thickness": 0.25,
        "pressure_z": -90.0, "reference": 0.3024,
        "modern_reference": 0.30059,
        "observable": "free-edge midpoint downward displacement",
    },
    "m3": {
        "R": 10.0, "cutout_degrees": 18.0,
        "E": 6.825e7, "nu": 0.3, "thickness": 0.04,
        "quarter_point_loads": [-1.0, 1.0],
        "reference": 0.09352155,
        "observable": "inward equator displacement at B",
    },
}


@dataclass
class BenchmarkCase:
    name: str
    mesh: CheckerboardMesh
    material: IsotropicMaterial
    boundary: BoundaryConditions
    load: torch.Tensor
    observable_weights: torch.Tensor
    reference_value: float | None
    metadata: dict


def ensure_tree() -> None:
    for name in (
        "01_mesh_distortion", "02_gauge_invariance", "03_scaling",
        "04_beta_sensitivity", "05_curved_spectrum", "06_boundary_modes",
    ):
        base = ROOT / name
        for child in ("config", "raw_data", "figures", "logs"):
            (base / child).mkdir(parents=True, exist_ok=True)


def write_json(path: Path, value) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(value, indent=2, sort_keys=True, allow_nan=True) + "\n")


def write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("")
        return
    keys: list[str] = []
    for row in rows:
        for key in row:
            if key not in keys:
                keys.append(key)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=keys, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def environment_record() -> dict:
    try:
        commit = subprocess.check_output(
            ["git", "rev-parse", "HEAD"], cwd=ROOT.parents[1], text=True).strip()
        worktree = subprocess.check_output(
            ["git", "status", "--short"], cwd=ROOT.parents[1], text=True).splitlines()
    except Exception:
        commit = "unavailable"
        worktree = ["unavailable"]
    return {
        "timestamp_local": time.strftime("%Y-%m-%dT%H:%M:%S%z"),
        "random_seed": SEED,
        "python": platform.python_version(),
        "platform": platform.platform(),
        "torch": torch.__version__,
        "numpy": np.__version__,
        "torch_threads": torch.get_num_threads(),
        "git_commit": commit,
        "git_worktree_status_at_recording": worktree,
        "dtype": str(DTYPE),
    }


def _topology(nx: int, ny: int, width: float, height: float) -> CheckerboardMesh:
    return CheckerboardMesh.rectangle(
        width=width, height=height, nx=nx, ny=ny, dtype=DTYPE)


def _distribution(count: int, ratio: float) -> np.ndarray:
    """Coordinates in [0, 1] with exactly requested max/min interval ratio."""
    if ratio <= 1.0 + 1.0e-14:
        return np.linspace(0.0, 1.0, count + 1)
    intervals = np.geomspace(1.0, float(ratio), count)
    values = np.concatenate(([0.0], np.cumsum(intervals)))
    return values / values[-1]


def parameter_grid(
    nx: int,
    ny: int,
    *,
    distortion: float = 0.0,
    distortion_kind: str = "none",
    nonuniform_ratio: float = 1.0,
) -> tuple[np.ndarray, np.ndarray]:
    x = _distribution(nx, nonuniform_ratio)
    y = _distribution(ny, nonuniform_ratio)
    yy, xx = np.meshgrid(y, x, indexing="ij")
    if distortion == 0.0 or distortion_kind == "none":
        return xx, yy
    interior = np.zeros_like(xx, dtype=bool)
    interior[1:-1, 1:-1] = True
    if distortion_kind == "smooth":
        dx = (distortion / nx) * np.sin(2.0*np.pi*xx) * np.sin(np.pi*yy)
        dy = (distortion / ny) * np.sin(np.pi*xx) * np.sin(2.0*np.pi*yy)
    elif distortion_kind == "random":
        rng = np.random.default_rng(SEED)
        dx = rng.normal(size=xx.shape)
        dy = rng.normal(size=yy.shape)
        dx /= math.sqrt(float(np.mean(dx[interior]**2)))
        dy /= math.sqrt(float(np.mean(dy[interior]**2)))
        # Match each component's *interior RMS* to the deterministic smooth
        # field at the same delta, instead of matching only its peak scale.
        smooth_x = np.sin(2.0*np.pi*xx)*np.sin(np.pi*yy)
        smooth_y = np.sin(np.pi*xx)*np.sin(2.0*np.pi*yy)
        dx *= (distortion / nx)*math.sqrt(float(np.mean(smooth_x[interior]**2)))
        dy *= (distortion / ny)*math.sqrt(float(np.mean(smooth_y[interior]**2)))
    else:
        raise ValueError(f"unknown distortion kind {distortion_kind!r}")
    xx = xx + np.where(interior, dx, 0.0)
    yy = yy + np.where(interior, dy, 0.0)
    return xx, yy


def _mesh_from_mapping(
    nx: int,
    ny: int,
    width: float,
    height: float,
    mapping: Callable[[torch.Tensor, torch.Tensor], torch.Tensor],
    **grid_options,
) -> CheckerboardMesh:
    base = _topology(nx, ny, width, height)
    xx, yy = parameter_grid(nx, ny, **grid_options)
    u = torch.as_tensor(width*xx, dtype=DTYPE)
    v = torch.as_tensor(height*yy, dtype=DTYPE)
    reference = mapping(u, v)
    return CheckerboardMesh(
        reference_vertices=reference, hx=width/nx, hy=height/ny,
        quads=base.quads, boundaries=base.boundaries)


def plate_mesh(nx: int, ny: int, **options) -> CheckerboardMesh:
    return _mesh_from_mapping(
        nx, ny, 1.0, 1.0,
        lambda u, v: torch.stack((u, v, torch.zeros_like(u)), -1),
        **options)


def roof_mesh(nx: int, ny: int, **options) -> CheckerboardMesh:
    data = BENCHMARKS["roof"]
    radius = data["R"]
    arc = radius*math.radians(2.0*data["half_angle_degrees"])
    half = math.radians(data["half_angle_degrees"])

    def mapping(u, v):
        theta = u/radius - half
        return torch.stack((radius*torch.sin(theta), v,
                            radius*torch.cos(theta)), -1)

    return _mesh_from_mapping(nx, ny, arc, data["L"], mapping, **options)


def m3_mesh(nx: int, ny: int, **options) -> CheckerboardMesh:
    data = BENCHMARKS["m3"]
    radius = data["R"]
    azimuth_arc = radius*math.pi/2.0
    meridian_arc = radius*math.radians(90.0-data["cutout_degrees"])
    cut = math.radians(data["cutout_degrees"])

    def mapping(u, v):
        phi = u/radius
        theta = cut + v/radius
        return torch.stack((radius*torch.sin(theta)*torch.cos(phi),
                            radius*torch.sin(theta)*torch.sin(phi),
                            radius*torch.cos(theta)), -1)

    return _mesh_from_mapping(
        nx, ny, azimuth_arc, meridian_arc, mapping, **options)


def curved_patch_mesh(kind: str, n: int) -> CheckerboardMesh:
    if kind == "sphere":
        radius = 3.0
        span = math.radians(50.0)

        def mapping(u, v):
            a = (u-0.5*radius*span)/radius
            b = (v-0.5*radius*span)/radius
            x = radius*torch.sin(a)*torch.cos(b)
            y = radius*torch.sin(b)
            z = radius*torch.cos(a)*torch.cos(b)
            return torch.stack((x, y, z), -1)
        return _mesh_from_mapping(n, n, radius*span, radius*span, mapping)
    if kind == "saddle":
        def mapping(u, v):
            x = u-0.5
            y = v-0.5
            return torch.stack((x, y, 0.6*(x*x-y*y)), -1)
        return _mesh_from_mapping(n, n, 1.0, 1.0, mapping)
    raise ValueError(kind)


def strip_mesh(kind: str, n: int) -> CheckerboardMesh:
    # ``n`` is elements per parametric direction; the physical 2:1 domain is a
    # strip without making the spectral state count four times larger.
    nx, ny = n, n
    if kind == "flat":
        return _mesh_from_mapping(
            nx, ny, 2.0, 1.0,
            lambda u, v: torch.stack((u, v, torch.zeros_like(u)), -1))
    if kind == "cylinder":
        radius = 2.0
        width = radius*math.radians(50.0)
        def mapping(u, v):
            theta = v/radius
            return torch.stack((u, radius*torch.sin(theta),
                                radius*(1.0-torch.cos(theta))), -1)
        return _mesh_from_mapping(nx, ny, 2.0, width, mapping)
    raise ValueError(kind)


def _empty_mask(mesh: CheckerboardMesh) -> torch.Tensor:
    return torch.zeros_like(mesh.reference_vertices, dtype=torch.bool)


def _pressure_load(mesh: CheckerboardMesh, pressure_z: float) -> torch.Tensor:
    """Lumped global-z dead load using physical reference quad areas."""
    q = mesh.reference_vertices
    values = torch.zeros_like(q)
    for j in range(mesh.ny):
        for i in range(mesh.nx):
            q00, q10 = q[j, i], q[j, i+1]
            q11, q01 = q[j+1, i+1], q[j+1, i]
            area = 0.5*torch.linalg.vector_norm(torch.linalg.cross(
                q10-q00, q11-q00, dim=-1))
            area += 0.5*torch.linalg.vector_norm(torch.linalg.cross(
                q11-q00, q01-q00, dim=-1))
            values[j:j+2, i:i+2, 2] += float(pressure_z)*area/4.0
    return values


def _center_weights(mesh: CheckerboardMesh) -> torch.Tensor:
    """Gauge-annihilating centered observable (partition of unity)."""
    weights = torch.zeros_like(mesh.reference_vertices)
    j, i = mesh.ny//2, mesh.nx//2
    weights[j, i, 2] = 0.5
    weights[j-1, i, 2] += 0.125
    weights[j+1, i, 2] += 0.125
    weights[j, i-1, 2] += 0.125
    weights[j, i+1, 2] += 0.125
    return weights


def _edge_midpoint_weights(mesh: CheckerboardMesh, component: int) -> torch.Tensor:
    weights = torch.zeros_like(mesh.reference_vertices)
    j = mesh.ny//2
    weights[j, 0, component] = 0.5
    weights[j-1, 0, component] += 0.25
    weights[j+1, 0, component] += 0.25
    return weights


def _corner_physical_weights(mesh: CheckerboardMesh, corner: str, component: int) -> torch.Tensor:
    """Average of the two incident physical edge midpoints at a corner."""
    weights = torch.zeros_like(mesh.reference_vertices)
    if corner == "lower_left":
        j, i, dj, di = mesh.ny, 0, -1, 1
    elif corner == "lower_right":
        j, i, dj, di = mesh.ny, mesh.nx, -1, -1
    else:
        raise ValueError(corner)
    weights[j, i, component] = 0.5
    weights[j+dj, i, component] += 0.25
    weights[j, i+di, component] += 0.25
    return weights


def benchmark_case(
    name: str,
    nx: int,
    ny: int | None = None,
    **mesh_options,
) -> BenchmarkCase:
    ny = nx if ny is None else ny
    if name == "plate":
        data = BENCHMARKS[name]
        mesh = plate_mesh(nx, ny, **mesh_options)
        mask = _empty_mask(mesh)
        mask[0, :, 2] = mask[-1, :, 2] = True
        mask[:, 0, 2] = mask[:, -1, 2] = True
        # Minimal physical in-plane rigid-motion removal.
        mask[0, 0, 0:2] = True
        mask[0, -1, 1] = True
        load = _pressure_load(mesh, data["pressure_z"])
        obs = -_center_weights(mesh)
    elif name == "roof":
        data = BENCHMARKS[name]
        mesh = roof_mesh(nx, ny, **mesh_options)
        mask = _empty_mask(mesh)
        # Rigid diaphragms on the two curved ends: global ux=uz=0.
        mask[0, :, 0] = mask[0, :, 2] = True
        mask[-1, :, 0] = mask[-1, :, 2] = True
        mask[0, 0, 1] = True  # remove axial rigid translation
        load = _pressure_load(mesh, data["pressure_z"])
        obs = -_edge_midpoint_weights(mesh, 2)
    elif name == "m3":
        data = BENCHMARKS[name]
        mesh = m3_mesh(nx, ny, **mesh_options)
        mask = _empty_mask(mesh)
        # Quarter-sector symmetry planes and one z anchor.
        mask[:, 0, 1] = True
        mask[:, -1, 0] = True
        mask[0, 0, 2] = True
        wa = _corner_physical_weights(mesh, "lower_left", 0)
        wb = _corner_physical_weights(mesh, "lower_right", 1)
        load = -wa + wb
        # B is the inward-loaded x-axis point in the conventional labeling.
        obs = -wa
    else:
        raise ValueError(name)
    material = IsotropicMaterial(
        young=data["E"], poisson=data["nu"], thickness=data["thickness"])
    return BenchmarkCase(
        name=name, mesh=mesh, material=material,
        boundary=BoundaryConditions(mask), load=load,
        observable_weights=obs, reference_value=data["reference"],
        metadata={"nx": nx, "ny": ny, **mesh_options})


def formulation_energy(
    positions: torch.Tensor,
    reference: torch.Tensor,
    mesh: CheckerboardMesh,
    targets,
    material: IsotropicMaterial,
    variant: str = "CB-Full",
    beta_multiplier: float = 1.0,
) -> tuple[torch.Tensor, dict]:
    parts = paper_components(
        positions, reference, mesh.hx, mesh.hy, targets, material)
    value = parts["membrane_b"] + parts["bending"]
    if variant in ("CB-M", "CB-Full"):
        value = value + float(beta_multiplier)*parts["membrane_xm"]
    if variant == "CB-Full":
        value = value + parts["bending_compatibility"]
    if variant not in ("CB-0", "CB-M", "CB-Full"):
        raise ValueError(variant)
    return value, parts


def internal_function(case: BenchmarkCase, variant="CB-Full", beta_multiplier=1.0):
    ref = case.mesh.reference_vertices
    targets = make_targets(ref, case.mesh.hx, case.mesh.hy)
    q0 = case.boundary.initial_coordinates(ref)

    def energy(q):
        value, _ = formulation_energy(
            case.boundary.positions(ref, q), ref, case.mesh, targets,
            case.material, variant, beta_multiplier)
        return value
    return energy, q0, targets


def dense_linear_response(
    case: BenchmarkCase,
    variant: str = "CB-Full",
    beta_multiplier: float = 1.0,
    eigen_tolerance: float = 1.0e-10,
) -> dict:
    """Reference-tangent solve with a spectral pseudoinverse.

    A pseudoinverse is deliberate: the raw checkerboard representation has an
    exact gauge.  Only load and observables that annihilate that gauge are used.
    """
    energy, q0, targets = internal_function(case, variant, beta_multiplier)
    started = time.perf_counter()
    tangent = torch.autograd.functional.hessian(energy, q0, vectorize=True)
    tangent = 0.5*(tangent+tangent.T)
    assembly = time.perf_counter()-started
    force = case.load[case.boundary.free_mask]
    started = time.perf_counter()
    eigenvalues, eigenvectors = torch.linalg.eigh(tangent)
    scale = max(float(eigenvalues.abs().max()), 1.0e-300)
    active = eigenvalues > eigen_tolerance*scale
    projected_force = eigenvectors[:, active].T @ force
    coordinates = eigenvectors[:, active] @ (projected_force/eigenvalues[active])
    solve_time = time.perf_counter()-started
    positions = case.boundary.positions(case.mesh.reference_vertices, coordinates)
    displacement = positions-case.mesh.reference_vertices
    observable = float((case.observable_weights*displacement).sum())
    _, parts = formulation_energy(
        positions, case.mesh.reference_vertices, case.mesh, targets,
        case.material, variant, beta_multiplier)
    residual = tangent@coordinates-force
    positive = eigenvalues[active]
    return {
        "coordinates": coordinates.detach(), "positions": positions.detach(),
        "tangent": tangent.detach(), "eigenvalues": eigenvalues.detach(),
        "observable": observable, "assembly_time_s": assembly,
        "solve_time_s": solve_time, "residual_l2": float(torch.linalg.norm(residual)),
        "condition_estimate": (float(positive[-1]/positive[0]) if len(positive) else math.inf),
        "nullity": int((~active).sum()), "parts": parts,
        "targets": targets,
    }


def mesh_quality(mesh: CheckerboardMesh) -> dict:
    q = mesh.reference_vertices.detach().cpu().numpy()
    scaled, aspects, angles = [], [], []
    inverted = 0
    for j in range(mesh.ny):
        for i in range(mesh.nx):
            points = [q[j, i], q[j, i+1], q[j+1, i+1], q[j+1, i]]
            normal = np.cross(points[1]-points[0], points[3]-points[0])
            normal_norm = np.linalg.norm(normal)
            if normal_norm == 0.0:
                inverted += 1
                continue
            normal /= normal_norm
            edges = [points[(k+1) % 4]-points[k] for k in range(4)]
            lengths = [np.linalg.norm(value) for value in edges]
            aspects.append(max(lengths)/min(lengths))
            for k in range(4):
                previous = points[(k-1) % 4]-points[k]
                following = points[(k+1) % 4]-points[k]
                cosine = np.dot(previous, following)/(np.linalg.norm(previous)*np.linalg.norm(following))
                angle = math.degrees(math.acos(float(np.clip(cosine, -1.0, 1.0))))
                angles.append(angle)
                # Polygon is counter-clockwise in parameter order. ``previous``
                # points from the corner back to its predecessor, so the
                # oriented incoming/outgoing pair is (-previous, following).
                jac = np.dot(np.cross(-previous, following), normal)
                scaled.append(jac/(np.linalg.norm(following)*np.linalg.norm(previous)))
                if jac <= 0.0:
                    inverted += 1
    from checkerboard_shell.geometry.core import torch_geometry
    geometry = torch_geometry(mesh.reference_vertices, mesh.hx, mesh.hy)
    determinants = torch.linalg.det(geometry.first_form)
    normals = geometry.normals_b
    adjacent = []
    if mesh.nx > 1:
        adjacent.append((normals[:, 1:]*normals[:, :-1]).sum(-1))
    if mesh.ny > 1:
        adjacent.append((normals[1:]*normals[:-1]).sum(-1))
    min_normal_dot = min(float(value.min()) for value in adjacent)
    return {
        "minimum_scaled_jacobian": float(min(scaled)) if scaled else math.nan,
        "maximum_aspect_ratio": float(max(aspects)) if aspects else math.nan,
        "minimum_interior_angle_deg": float(min(angles)) if angles else math.nan,
        "maximum_skew_angle_deg": float(max(abs(a-90.0) for a in angles)) if angles else math.nan,
        "inverted_corner_count": inverted,
        "minimum_b_metric_determinant": float(determinants.min()),
        "minimum_adjacent_normal_dot": min_normal_dot,
    }


def geometry_diagnostics(result: dict, case: BenchmarkCase) -> dict:
    geometry = result["parts"]["geometry"]
    det_b = torch.linalg.det(geometry.bulk.first_form)
    det_w = torch.linalg.det(geometry.first_form_w)
    ref_normals = result["targets"].checkerboard.boundary.geometry.bulk.normals_b
    normal_dot = (geometry.bulk.normals_b*ref_normals).sum(-1)
    total = float(result["parts"]["total"].detach())
    xm = float(result["parts"]["membrane_xm"].detach())
    xw = float(result["parts"]["bending_compatibility"].detach())
    return {
        "minimum_B_frame_metric_det": float(det_b.min()),
        "minimum_W_frame_metric_det": float(det_w.min()),
        "minimum_reference_normal_dot": float(normal_dot.min()),
        "normal_flip": bool(torch.any(normal_dot < 0.0)),
        "E_XM_over_E_int": xm/total if total else 0.0,
        "E_XW_over_E_int": xw/total if total else 0.0,
    }


def alternating_vector(shape: torch.Size) -> torch.Tensor:
    ny, nx, _ = shape
    jj, ii = torch.meshgrid(torch.arange(ny), torch.arange(nx), indexing="ij")
    parity = torch.where((ii+jj).remainder(2) == 0, 1.0, -1.0).to(DTYPE)
    return parity[..., None].expand(ny, nx, 3)


def invariant_metrics(
    reference: torch.Tensor,
    current: torch.Tensor,
    hx: float,
    hy: float,
    material: IsotropicMaterial,
    alphas: Iterable[float],
    direction: torch.Tensor,
) -> list[dict]:
    """Full-coordinate finite-configuration gauge check."""
    mesh = _topology(reference.shape[1]-1, reference.shape[0]-1,
                     hx*(reference.shape[1]-1), hy*(reference.shape[0]-1))
    mesh = CheckerboardMesh(reference, hx, hy, mesh.quads, mesh.boundaries)
    targets = make_targets(reference, hx, hy)
    parity = alternating_vector(reference.shape)
    gauge = parity*direction

    def quantities(vertices, need_tangent=True):
        variable = vertices.detach().clone().requires_grad_(True)
        parts = paper_components(variable, reference, hx, hy, targets, material)
        energy = parts["total"]
        residual, = torch.autograd.grad(energy, variable, create_graph=need_tangent)
        tangent = None
        if need_tangent:
            flat = variable.reshape(-1)
            tangent = torch.autograd.functional.hessian(
                lambda x: paper_components(
                    x.reshape_as(variable), reference, hx, hy, targets, material)["total"],
                flat, vectorize=True)
        mids = edge_midpoints(variable)
        return parts, residual.detach(), tangent.detach(), mids

    base, r0, k0, mids0 = quantities(current)
    characteristic_length = float(torch.linalg.vector_norm(
        reference.reshape(-1, 3).max(0).values
        - reference.reshape(-1, 3).min(0).values))
    norm_r = float(torch.linalg.norm(r0))
    norm_k = float(torch.linalg.norm(k0))
    base_energy = float(base["total"].detach())
    rows = []
    for alpha in alphas:
        shifted = current + float(alpha)*gauge
        parts, residual, tangent, mids = quantities(shifted)
        midpoint_change = math.sqrt(
            float(torch.linalg.norm(mids.horizontal-mids0.horizontal)**2
                  + torch.linalg.norm(mids.vertical-mids0.vertical)**2))
        a_change = torch.linalg.norm(parts["geometry"].bulk.first_form-base["geometry"].bulk.first_form)
        b_change = torch.linalg.norm(parts["geometry"].second_form-base["geometry"].second_form)
        energy_change = abs(float(parts["total"].detach())-base_energy)
        residual_change = float(torch.linalg.norm(residual-r0))
        rows.append({
            "alpha": float(alpha),
            "physical_midpoint_change": midpoint_change,
            "relative_a_B_change": float(a_change/max(float(torch.linalg.norm(base["geometry"].bulk.first_form)), 1.0)),
            "relative_b_W_change": float(b_change/max(float(torch.linalg.norm(base["geometry"].second_form)), 1.0)),
            "absolute_XM_change": abs(float(parts["X_M"].detach()-base["X_M"].detach())),
            "absolute_XWrel_change": float(torch.linalg.norm(parts["X_W_rel"]-base["X_W_rel"])),
            "absolute_internal_energy_change": energy_change,
            "relative_energy_change": (energy_change/abs(base_energy)
                                       if abs(base_energy)>1e-25 else math.nan),
            "absolute_internal_residual_change": residual_change,
            "relative_residual_change": (residual_change/norm_r
                                         if norm_r>1e-12 else math.nan),
            "residual_change_over_tangent_times_L": residual_change/
                max(norm_k*characteristic_length, 1e-300),
            "relative_tangent_change": float(torch.linalg.norm(tangent-k0))/max(norm_k, 1e-300),
        })
    return rows


def rss_megabytes() -> float:
    return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0


def timed_median(function: Callable, repeats: int = 3, warmups: int = 1) -> tuple[float, object]:
    result = None
    for _ in range(warmups):
        result = function()
    values = []
    for _ in range(repeats):
        started = time.perf_counter()
        result = function()
        values.append(time.perf_counter()-started)
    return float(np.median(values)), result


def save_figure(fig, base: Path) -> None:
    base.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(base.with_suffix(".png"), dpi=220, bbox_inches="tight")
    fig.savefig(base.with_suffix(".pdf"), bbox_inches="tight")
    fig.savefig(base.with_suffix(".svg"), bbox_inches="tight")


def plot_style() -> None:
    import matplotlib as mpl
    mpl.use("Agg")
    mpl.rcParams.update({
        "font.size": 8.0, "axes.labelsize": 8.0, "legend.fontsize": 8.0,
        "xtick.labelsize": 8.0, "ytick.labelsize": 8.0,
        "lines.linewidth": 1.25, "lines.markersize": 4.0,
        "figure.dpi": 120, "savefig.dpi": 220,
    })


def relative_error(value: float, reference: float | None) -> float:
    if reference is None or reference == 0.0:
        return math.nan
    return abs(value-reference)/abs(reference)


def state_counts(case: BenchmarkCase) -> tuple[int, int]:
    return case.mesh.num_vertices*3, int(case.boundary.free_mask.sum())


def method_unavailable_rows(method: str, contexts: Iterable[dict], reason: str) -> list[dict]:
    return [{**context, "method": method, "status": "UNAVAILABLE", "reason": reason}
            for context in contexts]
