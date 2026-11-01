"""Reproducible cProfile of old and shared-compliance aircraft operators."""
from __future__ import annotations

import cProfile
import itertools
import json
import pstats
import time

from experiments.engineering_inverse.aircraft_skin import (
    OUT,build_case,layout_operator,prepare_layout_compliance)
from experiments.engineering_inverse.common import checkerboard_displacement_operator


def _stats(profile):
    stats=pstats.Stats(profile).stats
    rows=[]
    for (filename,line,function),(primitive,calls,self_time,cumulative,callers) in stats.items():
        if any(word in function or word in filename for word in (
            "gstrf","splu","factorized","assemble_coloured_local_hessian",
            "low_rank_tangent","positional_beam_stiffness",
            "checkerboard_displacement_operator","prepare_layout_compliance",
            "layout_operator")):
            rows.append(dict(file=filename,line=line,function=function,
                             calls=calls,self_seconds=self_time,
                             cumulative_seconds=cumulative))
    return sorted(rows,key=lambda row:row["cumulative_seconds"],reverse=True)


def run(nx=20,ny=40):
    case=build_case(nx,ny)
    subsets=[subset for count in range(1,6)
             for subset in itertools.combinations(range(5),count)]
    result=dict(mesh=[nx,ny],layouts=len(subsets),methods={})
    for method in ("direct_per_layout", "shared_compliance"):
        profile=cProfile.Profile();start=time.perf_counter();profile.enable()
        base=checkerboard_displacement_operator(case["mesh"],case["material"],
            case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"])
        prepared=({"K":base["K"],"residual":base["residual"]}
                  if method=="direct_per_layout" else
                  prepare_layout_compliance(case,base))
        for subset in subsets:layout_operator(case,subset,prepared)
        profile.disable();elapsed=time.perf_counter()-start
        result["methods"][method]=dict(seconds=elapsed,profile=_stats(profile))
        print(method,elapsed,flush=True)
        (OUT/"profile_round2.json").write_text(json.dumps(result,indent=2)+"\n")
    return result


if __name__=="__main__":run()
