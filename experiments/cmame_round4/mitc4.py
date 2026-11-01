"""Validated MITC4 baseline via OpenSeesPy ``ShellMITC4``.

The element is the standard 4-node MITC4 shell in OpenSees
(``ElasticMembranePlateSection``, 6 DOF per node).  It is used here as an
external, provenance-recorded standard-shell baseline, not as a replacement
for the Checkerboard formulation.

The three problems use the same geometry, material, load, and physical
boundary conditions as the recovered manuscript drivers.  The MITC4 mesh is
the same structured quad mesh; the element is a Mindlin-Reissner shell, so
its plate response contains a small transverse-shear contribution that a
Kirchhoff--Love shell does not have.  This is reported explicitly.
"""

from __future__ import annotations

import json
import math
import sys
import time
from pathlib import Path

import numpy as np

from experiments.cmame_round4.common import DATA, REFERENCES, write_csv, write_json


def _ops():
    try:
        import openseespy.opensees as ops
    except ImportError:
        sys.path.insert(0, "/tmp/openseespy_pkg")
        import openseespy.opensees as ops
    return ops


def _solve(ops, model_name: str) -> None:
    ops.system("SparseGeneral")
    ops.numberer("RCM")
    ops.constraints("Plain")
    ops.integrator("LoadControl", 1.0)
    ops.algorithm("Linear")
    ops.analysis("Static")
    status = ops.analyze(1)
    if status != 0:
        raise RuntimeError(f"{model_name}: OpenSees analyze returned {status}")


def _plate_load(n: int, q: float) -> np.ndarray:
    """Consistent nodal load for uniform pressure q on a unit square."""
    values = np.zeros((n + 1, n + 1))
    area = (1.0 / n) ** 2
    for j in range(n):
        for i in range(n):
            for ii, jj in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
                values[jj, ii] += q * area / 4.0
    return values


def run_plate(n: int, *, E=1.0, nu=0.3, thickness=0.01,
              pressure=1.0e-8) -> dict:
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)

    def tag(i, j):
        return j * (n + 1) + i + 1

    for j in range(n + 1):
        for i in range(n + 1):
            ops.node(tag(i, j), i / n, j / n, 0.0)
    ops.section("ElasticMembranePlateSection", 1, E, nu, thickness, 0.0)
    element = 1
    for j in range(n):
        for i in range(n):
            ops.element("ShellMITC4", element, tag(i, j), tag(i + 1, j),
                        tag(i + 1, j + 1), tag(i, j + 1), 1)
            element += 1

    fixed = 0
    for j in range(n + 1):
        for i in range(n + 1):
            if i in (0, n) or j in (0, n):
                ops.fix(tag(i, j), 0, 0, 1, 0, 0, 0)
                fixed += 1
    ops.fix(tag(0, 0), 1, 1, 0, 0, 0, 0); fixed += 2
    ops.fix(tag(n, 0), 0, 1, 0, 0, 0, 0); fixed += 1

    loads = _plate_load(n, pressure)
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(n + 1):
        for i in range(n + 1):
            if abs(loads[j, i]) > 0.0:
                ops.load(tag(i, j), 0.0, 0.0, -loads[j, i], 0.0, 0.0, 0.0)

    started = time.perf_counter()
    _solve(ops, f"plate-{n}")
    elapsed = time.perf_counter() - started
    response = float(-ops.nodeDisp(tag(n // 2, n // 2), 3))
    vertices = (n + 1) ** 2
    return {
        "benchmark": "plate", "n": int(n), "method": "MITC4",
        "response": response, "reference": REFERENCES["plate"],
        "relative_error": abs(response - REFERENCES["plate"])
        / REFERENCES["plate"],
        "raw_scalar_states": int(6 * vertices),
        "active_states": int(6 * vertices - fixed),
        "nodes": int(vertices),
        "elements": int(n * n),
        "runtime_seconds": elapsed,
        "status": "SUCCESS",
        "element": "OpenSees ShellMITC4",
        "section": "ElasticMembranePlateSection",
        "notes": "Mindlin-Reissner shell; transverse-shear floor at t/a=0.01",
    }


def _roof_geometry(n: int, radius=25.0, half_angle=math.radians(40.0),
                   half_length=25.0):
    nodes = {}
    for j in range(n + 1):
        theta = -half_angle + 2.0 * half_angle * j / n
        for i in range(n + 1):
            x = half_length * i / n
            nodes[(i, j)] = (x, radius * math.sin(theta),
                             radius * math.cos(theta))
    return nodes


def _quad_area(points) -> float:
    p = np.asarray(points)
    return (0.5 * np.linalg.norm(np.cross(p[1] - p[0], p[3] - p[0]))
            + 0.5 * np.linalg.norm(np.cross(p[1] - p[2], p[3] - p[2])))


def run_roof(n: int, *, E=4.32e8, nu=0.0, thickness=0.25,
             pressure=90.0) -> dict:
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nodes = _roof_geometry(n)

    def tag(i, j):
        return j * (n + 1) + i + 1

    for (i, j), (x, y, z) in nodes.items():
        ops.node(tag(i, j), x, y, z)
    ops.section("ElasticMembranePlateSection", 1, E, nu, thickness, 0.0)
    element = 1
    for j in range(n):
        for i in range(n):
            ops.element("ShellMITC4", element, tag(i, j), tag(i + 1, j),
                        tag(i + 1, j + 1), tag(i, j + 1), 1)
            element += 1

    fixed = 0
    for j in range(n + 1):
        # Symmetry plane x=0: ux=0 and the two bending rotations vanish.
        ops.fix(tag(0, j), 1, 0, 0, 0, 1, 1)
        fixed += 3
        # Rigid diaphragm x=L/2: uy=uz=0.
        ops.fix(tag(n, j), 0, 1, 1, 0, 0, 0)
        fixed += 2

    loads = np.zeros((n + 1, n + 1))
    for j in range(n):
        for i in range(n):
            area = _quad_area([
                nodes[(i, j)], nodes[(i + 1, j)],
                nodes[(i + 1, j + 1)], nodes[(i, j + 1)]])
            for ii, jj in ((i, j), (i + 1, j), (i + 1, j + 1), (i, j + 1)):
                loads[jj, ii] += pressure * area / 4.0
    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    for j in range(n + 1):
        for i in range(n + 1):
            if abs(loads[j, i]) > 0.0:
                ops.load(tag(i, j), 0.0, 0.0, -loads[j, i], 0.0, 0.0, 0.0)

    started = time.perf_counter()
    _solve(ops, f"roof-{n}")
    elapsed = time.perf_counter() - started
    response = float(-ops.nodeDisp(tag(0, n), 3))
    vertices = (n + 1) ** 2
    return {
        "benchmark": "roof", "n": int(n), "method": "MITC4",
        "response": response,
        "reference": REFERENCES["roof_historical"],
        "relative_error": abs(response - REFERENCES["roof_historical"])
        / REFERENCES["roof_historical"],
        "reference_modern": REFERENCES["roof_modern"],
        "relative_error_modern": abs(response - REFERENCES["roof_modern"])
        / REFERENCES["roof_modern"],
        "raw_scalar_states": int(6 * vertices),
        "active_states": int(6 * vertices - fixed),
        "nodes": int(vertices),
        "elements": int(n * n),
        "runtime_seconds": elapsed,
        "status": "SUCCESS",
        "element": "OpenSees ShellMITC4",
        "section": "ElasticMembranePlateSection",
        "notes": "symmetry x=0: ux=0, thy=thz=0; diaphragm x=L/2: uy=uz=0",
    }


def _m3_geometry(n: int, radius=10.0, theta0=math.radians(18.0)):
    nodes = {}
    for j in range(n + 1):
        phi = 0.5 * math.pi * j / n
        for i in range(n + 1):
            theta = theta0 + (0.5 * math.pi - theta0) * i / n
            nodes[(i, j)] = (
                radius * math.sin(theta) * math.cos(phi),
                radius * math.sin(theta) * math.sin(phi),
                radius * math.cos(theta))
    return nodes


_SMOOTHED_POINT_WEIGHTS = {0: 15.0/16.0, 1: 5.0/16.0,
                          2: -7.0/16.0, 3: 3.0/16.0}


def run_m3(n: int, *, E=6.825e7, nu=0.3, thickness=0.04,
           load_mode: str = "raw") -> dict:
    ops = _ops()
    ops.wipe()
    ops.model("basic", "-ndm", 3, "-ndf", 6)
    nodes = _m3_geometry(n)

    def tag(i, j):
        return j * (n + 1) + i + 1

    for (i, j), (x, y, z) in nodes.items():
        ops.node(tag(i, j), x, y, z)
    ops.section("ElasticMembranePlateSection", 1, E, nu, thickness, 0.0)
    element = 1
    for j in range(n):
        for i in range(n):
            ops.element("ShellMITC4", element, tag(i, j), tag(i + 1, j),
                        tag(i + 1, j + 1), tag(i, j + 1), 1)
            element += 1

    fixed = 0
    for i in range(n + 1):
        # phi=0 symmetry plane (y=0): uy=0, thx=thz=0.
        ops.fix(tag(i, 0), 0, 1, 0, 1, 0, 1)
        fixed += 3
        # phi=90 symmetry plane (x=0): ux=0, thy=thz=0.
        ops.fix(tag(i, n), 1, 0, 0, 0, 1, 1)
        fixed += 3
    # Remove the remaining rigid z translation.  The node already has
    # uy, thx, thz fixed, so only uz is added here.
    ops.fix(tag(n // 2, 0), 0, 0, 1, 0, 0, 0)
    fixed += 1

    ops.timeSeries("Linear", 1)
    ops.pattern("Plain", 1, 1)
    if load_mode == "raw":
        ops.load(tag(n, 0), 1.0, 0.0, 0.0, 0.0, 0.0, 0.0)
        ops.load(tag(n, n), 0.0, -1.0, 0.0, 0.0, 0.0, 0.0)
    elif load_mode == "smoothed":
        # Same one-sided physical-point distribution used by the Checkerboard
        # driver, mapped from vertical edge midpoints to the adjacent raw
        # nodes (each midpoint weight is split equally between its two nodes).
        for offset, weight in _SMOOTHED_POINT_WEIGHTS.items():
            ops.load(tag(n, offset), weight, 0.0, 0.0, 0.0, 0.0, 0.0)
            ops.load(tag(n, n - offset), 0.0, -weight, 0.0, 0.0, 0.0, 0.0)
    else:
        raise ValueError(load_mode)

    started = time.perf_counter()
    _solve(ops, f"m3-{n}")
    elapsed = time.perf_counter() - started
    response = float(-ops.nodeDisp(tag(n, n), 2))
    vertices = (n + 1) ** 2
    return {
        "benchmark": "M3", "n": int(n), "method": "MITC4",
        "response": response, "reference": REFERENCES["M3"],
        "relative_error": abs(response - REFERENCES["M3"])
        / REFERENCES["M3"],
        "raw_scalar_states": int(6 * vertices),
        "active_states": int(6 * vertices - fixed),
        "nodes": int(vertices),
        "elements": int(n * n),
        "runtime_seconds": elapsed,
        "status": "SUCCESS",
        "element": "OpenSees ShellMITC4",
        "section": "ElasticMembranePlateSection",
        "load_mode": load_mode,
        "notes": (
            "symmetry rotations constrained; "
            + ("raw equator corner point load (same discrete load as the "
               "LibShell MidedgeTan driver)" if load_mode == "raw"
               else "smoothed physical-point load (same one-sided "
                    "distribution as the Checkerboard driver)")),
    }


def run_plate_thin(n: int, thickness: float) -> dict:
    """Thin-limit plate check: w scales as 1/t^3, so rescale the reference."""
    row = run_plate(n, thickness=thickness)
    reference = REFERENCES["plate"] * (0.01 / thickness) ** 3
    row.update({
        "benchmark": "plate_thin",
        "response": row["response"],
        "reference": reference,
        "relative_error": abs(row["response"] - reference) / reference,
        "thickness": thickness,
        "notes": "Mindlin shear is negligible in this thin limit",
    })
    return row


def run() -> dict:
    rows = []
    for n in (8, 16, 32, 64):
        rows.append(run_plate(n))
    for thickness in (0.001, 0.0001):
        for n in (32, 64):
            rows.append(run_plate_thin(n, thickness))
    for n in (12, 16, 24, 32):
        rows.append(run_roof(n))
    for n in (6, 8, 12, 16, 24, 32, 48, 64):
        rows.append(run_m3(n, load_mode="raw"))
    for n in (6, 8, 12, 16, 24, 32):
        rows.append(run_m3(n, load_mode="smoothed"))
    write_csv(DATA/"mitc4.csv", rows)
    write_json(DATA/"mitc4.json", rows)
    lines = [
        "# MITC4 baseline (OpenSeesPy ShellMITC4)",
        "",
        "| benchmark | n | MITC4 response | error | raw states | active states |",
        "|---|---:|---:|---:|---:|---:|",
    ]
    for row in rows:
        if row["benchmark"] == "roof":
            error = f'{100.0*row["relative_error"]:.3f}% / {100.0*row["relative_error_modern"]:.3f}%'
        else:
            error = f'{100.0*row["relative_error"]:.3f}%'
        label = row["benchmark"]
        if row["benchmark"] == "plate_thin":
            label += f' (t={row["thickness"]:g})' 
        if row.get("load_mode"):
            label += f' ({row["load_mode"]})'
        lines.append(
            f'| {label} | {row["n"]} | {row["response"]:.8g} | '
            f'{error} | {row["raw_scalar_states"]} | {row["active_states"]} |')
    lines.extend([
        "",
        "Roof errors are relative to the historical `0.3024` and modern "
        "`0.3005924566` references, respectively.",
        "",
        "MITC4 is a Mindlin-Reissner shell with 6 DOF per node; its plate "
        "response contains a small transverse-shear contribution that the "
        "Kirchhoff--Love Checkerboard/MidedgeTan models do not have.",
    ])
    (DATA.parent/"MITC4_TABLE.md").write_text("\n".join(lines) + "\n",
                                             encoding="utf-8")
    return {"rows": len(rows), "failed": sum(row["status"] != "SUCCESS" for row in rows)}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
