"""One isolated scaling measurement (invoked by :mod:`run_scaling`)."""

from __future__ import annotations

import argparse
import json
import os
import resource
import subprocess
import time

import numpy as np
import scipy.sparse.linalg as spla
import torch

from checkerboard_shell.mechanics.low_rank import LowRankTangent, assemble_paper_low_rank_tangent
from experiments.cmame_round2.recovered_drivers import (
    M3_REFERENCE,
    PLATE_REFERENCE,
    m3_problem,
    plate_stiffness,
)


THREAD_COUNT = int(os.environ.get("CB_SCALING_THREADS", "8"))
torch.set_num_threads(THREAD_COUNT)
try:
    torch.set_num_interop_threads(1)
except RuntimeError:
    pass


def _median(values):
    return float(np.median(np.asarray(values, dtype=float)))


def plate(n: int, variant: str) -> dict:
    include_xw = variant == "CB-Full"
    # Warm-up is deliberately excluded.
    plate_stiffness(n, include_xw=include_xw)
    assembly = []
    matrix = None
    for _ in range(3):
        start = time.perf_counter()
        matrix = plate_stiffness(n, include_xw=include_xw)
        assembly.append(time.perf_counter()-start)
    grid = np.arange((n+1)**2).reshape(n+1, n+1)
    active = grid[1:-1, 1:-1].reshape(-1)
    reduced = matrix[active][:, active].tocsc()
    force = np.full(active.size, 1.0e-8/(n*n))
    factor_times, solve_times = [], []
    solution = None
    for _ in range(3):
        start = time.perf_counter(); factor = spla.splu(reduced)
        factor_times.append(time.perf_counter()-start)
        start = time.perf_counter(); solution = factor.solve(force)
        solve_times.append(time.perf_counter()-start)
    center = float(solution[np.flatnonzero(active == grid[n//2, n//2])[0]])
    start = time.perf_counter(); residual = reduced@solution-force
    residual_time = time.perf_counter()-start
    return {
        "benchmark": "plate", "n": n, "method": variant, "status": "SUCCESS",
        "raw_scalar_states": 3*(n+1)**2, "active_states": int(active.size),
        "local_tangent_nnz": int(reduced.nnz), "rank_one_storage": 0,
        "augmented_system_nnz": int(reduced.nnz),
        "assembly_times_s": assembly,
        "assembly_time_median_s": _median(assembly),
        "residual_time_median_s": residual_time,
        "local_hessian_time_median_s": _median(assembly),
        "rank_vector_time_median_s": 0.0,
        "factorization_time_median_s": _median(factor_times),
        "factorization_times_s": factor_times,
        "linear_solve_time_median_s": _median(solve_times),
        "linear_solve_times_s": solve_times,
        "observable": center,
        "relative_error": abs(center-PLATE_REFERENCE)/PLATE_REFERENCE,
        "linear_residual_l2": float(np.linalg.norm(residual)),
        "repeat_count": 3, "warmup_count": 1,
        "torch_thread_count": THREAD_COUNT,
        "note": "X_M is identically inactive in the manuscript scalar transverse branch",
    }


def m3(n: int, variant: str) -> dict:
    problem = m3_problem(n)
    reference = problem["reference"]
    zero = torch.zeros_like(reference)
    # Warm-up is deliberately excluded.
    assemble_paper_low_rank_tangent(
        zero, reference, problem["hx"], problem["hy"], problem["targets"],
        problem["material"], variant=variant)
    audits, totals = [], []
    tangent = residual = None
    for _ in range(3):
        start = time.perf_counter()
        tangent, residual, audit = assemble_paper_low_rank_tangent(
            zero, reference, problem["hx"], problem["hy"], problem["targets"],
            problem["material"], variant=variant)
        totals.append(time.perf_counter()-start); audits.append(audit)
    factor_times, solve_times, augmented = [], [], []
    solution = None
    status = "SUCCESS"
    reason = ""
    try:
        for _ in range(3):
            solution, solve_audit = tangent.solve(
                problem["force"], constraints=problem["constraints"])
            factor_times.append(solve_audit["factorization_seconds"])
            solve_times.append(solve_audit["linear_solve_seconds"])
            augmented.append(solve_audit["augmented_nnz"])
    except Exception as error:
        status, reason = "FAILED", f"{type(error).__name__}: {error}"
        solution = np.full(reference.numel(), np.nan)
    observable = -float(problem["point_b"][1]@solution)
    return {
        "benchmark": "M3", "n": n, "method": variant, "status": status,
        "reason": reason,
        "raw_scalar_states": int(reference.numel()),
        "active_states": int(reference.numel()-problem["constraints"].shape[0]),
        "local_tangent_nnz": int(tangent.local.nnz),
        "rank_one_storage": int(tangent.rank_one_storage if tangent.coefficient else 0),
        "augmented_system_nnz": int(_median(augmented)) if augmented else None,
        "assembly_times_s": totals,
        "assembly_time_median_s": _median(totals),
        "residual_times_s": [a["residual_seconds"] for a in audits],
        "residual_time_median_s": _median([a["residual_seconds"] for a in audits]),
        "local_hessian_time_median_s": _median([
            a["graph_seconds"]+a["hvp_seconds"] for a in audits]),
        "rank_vector_time_median_s": _median([a["rank_one_vector_seconds"] for a in audits]),
        "local_hessian_times_s": [
            a["graph_seconds"]+a["hvp_seconds"] for a in audits],
        "rank_vector_times_s": [a["rank_one_vector_seconds"] for a in audits],
        "factorization_time_median_s": _median(factor_times) if factor_times else None,
        "factorization_times_s": factor_times,
        "linear_solve_time_median_s": _median(solve_times) if solve_times else None,
        "linear_solve_times_s": solve_times,
        "observable": observable,
        "relative_error": abs(observable-M3_REFERENCE)/M3_REFERENCE,
        "reference_residual_inf": float(np.max(np.abs(residual))),
        "repeat_count": 3, "warmup_count": 1,
        "torch_thread_count": THREAD_COUNT,
    }


def baseline(benchmark: str, n: int) -> dict:
    if benchmark == "plate":
        command = ["/tmp/platepinn_checkerboard_benchmark_build/plate_linear_driver",
                   "MidedgeTan", str(n)]
        reference, key = PLATE_REFERENCE, "center_deflection"
    else:
        command = ["/tmp/platepinn_checkerboard_thin_limit_build/libshell_linear_thickness_driver",
                   str(n), "0.04"]
        reference, key = M3_REFERENCE, "minus_u_B"
    _ = subprocess.run(command, check=True, text=True, capture_output=True)
    rows, wall_times = [], []
    for _ in range(3):
        started = time.perf_counter()
        completed = subprocess.run(command, check=True, text=True, capture_output=True)
        wall_times.append(time.perf_counter()-started)
        rows.append(json.loads(completed.stdout.strip().splitlines()[-1]))
    row = rows[-1]
    assembly_values = [r.get("assembly_seconds") for r in rows
                       if r.get("assembly_seconds") is not None]
    states = (row.get("total_scalar_dof") if benchmark == "plate"
              else 6*n*n+8*n+3)
    return {
        "benchmark": benchmark, "n": n, "method": "MidedgeTan",
        "status": "SUCCESS", "raw_scalar_states": int(states),
        "active_states": int(row.get("active_linear_dof", row.get("active_dof"))),
        "local_tangent_nnz": (int(row["hessian_nnz"])
                              if "hessian_nnz" in row else None),
        "rank_one_storage": 0,
        "augmented_system_nnz": (int(row["hessian_nnz"])
                                 if "hessian_nnz" in row else None),
        "assembly_time_median_s": (_median(assembly_values)
                                   if assembly_values else None),
        "assembly_times_s": [r.get("assembly_seconds") for r in rows],
        "residual_time_median_s": None, "local_hessian_time_median_s": None,
        "rank_vector_time_median_s": 0.0,
        "factorization_time_median_s": None,
        "linear_solve_time_median_s": (_median([r["solve_seconds"] for r in rows])
                                       if "solve_seconds" in row else None),
        "total_process_time_median_s": _median(wall_times),
        "total_process_times_s": wall_times,
        "observable": float(row[key]),
        "relative_error": abs(float(row[key])-reference)/reference,
        "repeat_count": 3, "warmup_count": 1,
        "torch_thread_count": THREAD_COUNT,
        "timing_comparability": "auxiliary only: C++ LibShell backend differs",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("benchmark", choices=("plate", "M3"))
    parser.add_argument("n", type=int)
    parser.add_argument("method", choices=("CB-0", "CB-M", "CB-Full", "MidedgeTan"))
    args = parser.parse_args()
    if args.method == "MidedgeTan":
        result = baseline(args.benchmark, args.n)
    elif args.benchmark == "plate":
        result = plate(args.n, args.method)
    else:
        result = m3(args.n, args.method)
    result["peak_rss_mb"] = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss/1024.0
    child = resource.getrusage(resource.RUSAGE_CHILDREN).ru_maxrss/1024.0
    result["child_peak_rss_mb"] = child if child else None
    print(json.dumps(result, allow_nan=True))


if __name__ == "__main__":
    main()
