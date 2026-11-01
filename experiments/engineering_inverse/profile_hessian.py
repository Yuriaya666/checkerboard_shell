"""Reproduce exact-tangent AD profile without changing the shell model."""
from __future__ import annotations

import cProfile
import io
import json
import pstats
import time
from pathlib import Path

import torch

from checkerboard_shell.mechanics.energy import make_targets
from checkerboard_shell.mechanics.low_rank import assemble_paper_low_rank_tangent
from experiments.engineering_inverse.aircraft_skin import build_case

OUT = Path(__file__).resolve().parent / "results"


def run(nx=20, ny=40, strategy="legacy"):
    case = build_case(nx, ny)
    mesh = case["mesh"]
    ref = mesh.reference_vertices
    target = make_targets(ref, mesh.hx, mesh.hy)
    displacement = torch.zeros_like(ref)
    profile = cProfile.Profile()
    start = time.perf_counter()
    profile.enable()
    tangent, force, audit = assemble_paper_low_rank_tangent(
        displacement, ref, mesh.hx, mesh.hy, target, case["material"],
        colour_batch=7 if strategy=="legacy" else 21,
        assembly_strategy=strategy)
    profile.disable()
    text = io.StringIO()
    pstats.Stats(profile, stream=text).sort_stats("cumtime").print_stats(25)
    path = OUT/f"hessian_profile_{nx}x{ny}_{strategy}.txt"
    path.write_text(text.getvalue())
    summary = dict(mesh=[nx,ny],strategy=strategy,
                   elapsed_seconds=time.perf_counter()-start,
                   audit=audit,local_nnz=int(tangent.local.nnz),
                   residual_inf=float(force.abs().max()) if hasattr(force,"abs")
                   else float(abs(force).max()),profile_file=str(path))
    (OUT/f"hessian_profile_{nx}x{ny}_{strategy}.json").write_text(
        json.dumps(summary,indent=2)+"\n")
    print(json.dumps(summary,indent=2))
    return summary


if __name__ == "__main__":
    import argparse
    parser=argparse.ArgumentParser()
    parser.add_argument("--nx",type=int,default=20)
    parser.add_argument("--ny",type=int,default=40)
    parser.add_argument("--strategy",choices=("legacy","cached"),default="legacy")
    args=parser.parse_args()
    run(args.nx,args.ny,args.strategy)
