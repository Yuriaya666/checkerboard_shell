"""Four-variant ablation and prescribed exact cylinder isometry."""
from __future__ import annotations

import csv
import json
import time
from pathlib import Path

import numpy as np
import scipy.linalg as la
import torch

from checkerboard_shell.material.isotropic import IsotropicMaterial
from checkerboard_shell.mechanics.energy import make_targets, paper_components
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.cmame_round2.recovered_drivers import roof_problem, m3_problem

from .run_scan import DATA, save_rows


def ablation() -> None:
    rows = []
    for benchmark, n in (("roof", 8), ("roof", 12), ("M3", 6), ("M3", 8)):
        problem = roof_problem(n) if benchmark == "roof" else m3_problem(n)
        if benchmark == "M3":
            force = np.zeros(problem["reference"].numel())
            force[3*n*(n+1)] = 1.0
            force[3*(n*(n+1)+n)+1] = -1.0
        else:
            force = problem["force"]
        for variant in ("CB-0", "CB-M", "CB-W", "CB-Full"):
            start = time.perf_counter()
            tangent, _, audit = assemble_paper_low_rank_tangent(
                torch.zeros_like(problem["reference"]), problem["reference"],
                problem["hx"], problem["hy"], problem["targets"],
                problem["material"], variant=variant)
            dense = tangent.local.toarray() + tangent.coefficient * np.outer(tangent.vector, tangent.vector)
            if benchmark == "roof":
                active = problem["active"]
                matrix = dense[np.ix_(active, active)]
                rhs = force
                solution = la.solve(matrix, rhs, assume_a="sym")
                full = np.zeros(dense.shape[0]); full[active] = solution
                response = -full[problem["output_raw"]]
            else:
                constraint = problem["constraints"].toarray()
                null = la.null_space(constraint)
                matrix = null.T @ dense @ null
                reduced_force = null.T @ force
                reduced = la.solve(matrix, reduced_force, assume_a="sym")
                full = null @ reduced
                response = -full[3*(n*(n+1)+n)+1]
            spectrum = la.eigvalsh(matrix)
            positive = spectrum[spectrum > max(abs(spectrum).max()*1e-10, 1e-12)]
            # Sample the quadratic energy of the reference tangent, rather
            # than interpreting a finite deformation of a linear solve.
            amplitude = 1.0e-4
            positions = problem["reference"] + amplitude * torch.tensor(full).reshape_as(problem["reference"])
            with torch.no_grad():
                parts = paper_components(positions, problem["reference"],
                    problem["hx"], problem["hy"], problem["targets"], problem["material"])
            row = {"benchmark": benchmark, "n": n, "variant": variant,
                   "response": float(response), "raw_states": dense.shape[0],
                   "reduced_states": matrix.shape[0],
                   "near_nullity": int(len(spectrum)-len(positive)),
                   "smallest_positive_eigenvalue": float(positive[0]) if len(positive) else None,
                   "second_positive_eigenvalue": float(positive[1]) if len(positive)>1 else None,
                   "membrane_base_quadratic": float(parts["membrane_b"])/amplitude**2,
                   "xm_quadratic": float(parts["membrane_xm"])/amplitude**2,
                   "bending_base_quadratic": float(parts["bending"])/amplitude**2,
                   "xw_quadratic": float(parts["bending_compatibility"])/amplitude**2,
                   "time_s": time.perf_counter()-start}
            rows.append(row)
            save_rows(DATA/"completion_ablation.csv", rows)
            print(benchmark, n, variant, response, row["near_nullity"], flush=True)


def cylinder() -> None:
    rows = []
    radius = 1.5
    for n in (8, 12, 16, 24, 32):
        s = torch.linspace(-0.5, 0.5, n+1, dtype=torch.float64)
        y = torch.linspace(-0.5, 0.5, n+1, dtype=torch.float64)
        ss, yy = torch.meshgrid(s, y, indexing="ij")
        zero = torch.zeros_like(ss)
        reference = torch.stack((ss, yy, zero), -1)
        positions = torch.stack((radius*torch.sin(ss/radius), yy,
                                 radius*(1-torch.cos(ss/radius))), -1)
        targets = make_targets(reference, 1/n, 1/n)
        for thickness in (0.01, 0.0025):
            parts = paper_components(positions, reference, 1/n, 1/n,
                                      targets, IsotropicMaterial(1.0, 0.3, thickness))
            first = parts["geometry"].bulk.first_form
            baseline = targets.checkerboard.boundary.first_form_b
            metric_error = float(torch.linalg.vector_norm(first-baseline) /
                                 torch.linalg.vector_norm(baseline))
            em = float(parts["membrane"])
            eb = float(parts["bending"]+parts["bending_compatibility"])
            total = em + eb
            rows.append({"n": n, "thickness": thickness, "radius": radius,
                         "metric_relative_frobenius": metric_error,
                         "membrane_energy": em, "bending_energy": eb,
                         "membrane_to_bending": em/eb,
                         "xm_fraction": float(parts["membrane_xm"])/total,
                         "xw_fraction": float(parts["bending_compatibility"])/total,
                         "interpretation": "prescribed exact continuum isometry; not equilibrium"})
    save_rows(DATA/"cylinder_isometry.csv", rows)


if __name__ == "__main__":
    ablation()
    cylinder()
