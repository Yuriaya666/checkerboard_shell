"""Same-control cold nonlinear forward costs, separate from optimizer paths."""
from __future__ import annotations

import json
import resource
import time

import numpy as np
import scipy.sparse as sp

from experiments.engineering_inverse.automotive_panel import build_case
from experiments.engineering_inverse.automotive_spring_design import OUT
from experiments.engineering_inverse.mount_springs import checkerboard_forward,mitc4_forward


def asd_shell_graph_nnz(case):
    n=6*case["mesh"].num_vertices
    rows=[];cols=[]
    for quad in case["mesh"].quads:
        ids=[6*int(node)+component for node in quad for component in range(6)]
        rows.extend(np.repeat(ids,len(ids)))
        cols.extend(np.tile(ids,len(ids)))
    graph=sp.coo_matrix((np.ones(len(rows),dtype=bool),(rows,cols)),
                        shape=(n,n)).tocsr()
    fixed=[6*int(node)+component for node in case["datum"] for component in range(3)]
    fixed.append(6*(case["gauge_x_dof"]//3))
    free=np.setdiff1d(np.arange(n),fixed)
    return int(graph[free][:,free].nnz)


def run(repeats=3,nx=12,ny=8):
    saved=json.loads((OUT/"automotive_efficiency_cb_consistent_12x8_zero.json").read_text())
    control=np.asarray(saved["adjustment_mm"])*.001
    case=build_case(nx,ny)
    result=dict(mesh=[nx,ny],input="final matched-locator Checkerboard consistent design",
        repeats=repeats,asd_structural_graph_nnz=asd_shell_graph_nnz(case),
        methods={})
    for method in ("cb","asd"):
        rows=[]
        for _ in range(repeats):
            t=time.perf_counter()
            if method=="cb":
                r=checkerboard_forward(case,control,initial=None)
                if r["status"]!="CONVERGED":raise RuntimeError(r["status"])
                row=dict(status=r["status"],newton_iterations=r["newton_iterations_total"],
                    tangent_assemblies=r["tangent_assemblies_total"],
                    tangent_assembly_seconds=r["tangent_assembly_seconds"],
                    numeric_factorization_seconds=r["numeric_factorization_seconds"],
                    local_nnz=r["final_shell_local_nnz"],
                    rank_one_storage=r["rank_one_storage"])
            else:
                r=mitc4_forward(case,control,corotational=True,
                                analysis_steps=2,norm_type=0,gauge_x_fix=True)
                if r["status"]!=0:raise RuntimeError(r["status"])
                row=dict(status=r["status"],newton_iterations=r["newton_iterations_total"],
                    model_construction_seconds=r["model_construction_seconds"],
                    nonlinear_analysis_seconds=r["analysis_seconds"],
                    system_size=r["system_size"])
            row.update(elapsed_seconds=time.perf_counter()-t,
                max_reaction_N=float(np.max(np.linalg.norm(r["mount_reaction"],axis=1))),
                peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss))
            rows.append(row)
        result["methods"][method]=rows
        (OUT/"automotive_efficiency_forward.json").write_text(
            json.dumps(result,indent=2)+"\n")
        print(method,[round(row["elapsed_seconds"],3) for row in rows],flush=True)
    return result


if __name__=="__main__":run()
