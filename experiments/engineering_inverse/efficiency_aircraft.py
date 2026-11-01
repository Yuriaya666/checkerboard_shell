"""Four-strategy, same-task aircraft linear layout efficiency experiment."""
from __future__ import annotations

import argparse
import itertools
import json
import os
import resource
import subprocess
import sys
import time
from pathlib import Path

import numpy as np
import scipy.sparse as sp

from experiments.engineering_inverse.aircraft_skin import (
    OUT, PARAMS, build_case, layout_operator as cb_layout,
    optimize_controls, prepare_layout_compliance)
from experiments.engineering_inverse.common import (
    checkerboard_displacement_operator,mitc4_prescribed_response)
from experiments.engineering_inverse.mitc4_shared import (
    build_shared,layout_operator as mt_layout)
from experiments.engineering_inverse.aircraft_sparse import sparse_base_operator


METHODS=("cb_direct","cb_shared","mitc4_direct","mitc4_shared",
         "cb_shared_sparse","cb_shared_sparse_legacy")


def _subsets():
    return [tuple(s) for count in range(1,6)
            for s in itertools.combinations(range(5),count)]


def _mitc4_direct_layout(case,selected):
    groups=[[3*node+2 for node in case["patches"][station]] for station in selected]
    fixed=sorted(set(case["spar"]+[d for group in groups for d in group]))
    n=len(groups)
    base,_,r0=mitc4_prescribed_response(case["mesh"],case["material"],
        fixed,groups,np.zeros(n),case["O"],case["load"],beam_lines=case["beam_lines"])
    responses=[];reactions=[]
    for k in range(n):
        d=np.zeros(n);d[k]=.001
        response,_,reaction=mitc4_prescribed_response(
            case["mesh"],case["material"],fixed,groups,d,
            case["O"],case["load"],beam_lines=case["beam_lines"])
        responses.append((response-base)/.001)
        reactions.append((reaction-r0)/.001)
    RM=np.column_stack(reactions)
    lookup={int(d):k for k,d in enumerate(fixed)}
    net0=np.array([r0[[lookup[d] for d in group]].sum() for group in groups])
    netmap=np.array([RM[[lookup[d] for d in group]].sum(axis=0)
                     for group in groups])
    return dict(response0=base,response_map=np.column_stack(responses),
                actuator_reaction0=net0,actuator_reaction_map=netmap)


def _mitc4_graph_nnz(case):
    """Potential nonzero pattern from element connectivity, not numeric nnz."""
    n=6*case["mesh"].num_vertices
    rows=[];cols=[]
    for quad in case["mesh"].quads:
        ids=[6*int(node)+c for node in quad for c in range(6)]
        rows.extend(np.repeat(ids,len(ids)));cols.extend(np.tile(ids,len(ids)))
    for line,*_ in case["beam_lines"]:
        for a,b in zip(line[:-1],line[1:]):
            ids=[6*int(node)+c for node in (a,b) for c in range(6)]
            rows.extend(np.repeat(ids,len(ids)));cols.extend(np.tile(ids,len(ids)))
    matrix=sp.coo_matrix((np.ones(len(rows),dtype=bool),(rows,cols)),
                         shape=(n,n)).tocsr()
    spar=[6*(d//3)+d%3 for d in case["spar"]]
    keep=np.setdiff1d(np.arange(n),spar)
    return int(matrix[keep][:,keep].nnz)


def worker(nx,ny,method):
    started=time.perf_counter()
    case=build_case(nx,ny)
    case_construction_seconds=time.perf_counter()-started
    record=dict(mesh=[nx,ny],method=method,counts=dict(layouts=31,
        optimized_cases=62,forward_analyses=0,numeric_factorizations=None))
    if method.startswith("cb"):
        base=(sparse_base_operator(case,
                hessian_assembly_strategy=("legacy" if method.endswith("_legacy")
                                           else "cached"))
              if method.startswith("cb_shared_sparse") else
            checkerboard_displacement_operator(case["mesh"],case["material"],
                case["spar"],[],case["O"],case["load"],beam_lines=case["beam_lines"]))
        if method in ("cb_shared","cb_shared_sparse","cb_shared_sparse_legacy"):
            t=time.perf_counter();prepared=prepare_layout_compliance(case,base)
            shared_seconds=time.perf_counter()-t
            record["counts"]["numeric_factorizations"]=1
        else:
            prepared={"K":base["K"],"residual":base["residual"]}
            shared_seconds=0.0
            record["counts"]["numeric_factorizations"]=32
        record.update(total_state_dofs=base["total_dofs"],
            free_state_dofs=base["free_dofs"],
            local_shell_nnz=base["shell_local_nnz"],
            full_local_nnz=base.get("full_local_nnz"),
            rank_one_storage=base["rank_one_storage"],
            matrix_stored_entries=(base.get("augmented_nnz") if method=="cb_shared_sparse"
                                   else base["total_dofs"]**2),
            base_seconds=base["assembly_seconds"],
            shell_tangent_seconds=base["shell_tangent_seconds"],
            shell_tangent_audit=base.get("shell_tangent_audit"),
            base_system_build_seconds=base.get("system_build_seconds",
                base.get("beam_build_seconds")),
            base_factor_seconds=base["factor_seconds"],
            augmented_nnz=base.get("augmented_nnz"),
            candidate_response_seconds=shared_seconds)
        make=lambda subset:cb_layout(case,subset,prepared)
    elif method=="mitc4_shared":
        shared=build_shared(case)
        record.update(total_state_dofs=shared["total_dofs"],
            free_state_dofs=shared["free_dofs"],
            base_seconds=shared["construction_seconds"]+
                         shared["first_analysis_seconds"],
            model_construction_seconds=shared["construction_seconds"],
            base_factor_and_solve_seconds=shared["first_analysis_seconds"],
            candidate_response_seconds=shared["candidate_response_seconds"],
            candidate_count=shared["candidate_count"],
            numFact_reported=shared["factorization_count_reported"])
        record["counts"]["forward_analyses"]=shared["candidate_count"]+1
        record["counts"]["numeric_factorizations"]="factorOnce; counter unavailable"
        make=lambda subset:mt_layout(case,shared,subset)
    else:
        record.update(total_state_dofs=6*case["mesh"].num_vertices,
            free_state_dofs=6*case["mesh"].num_vertices-len(case["spar"]),
            base_seconds=0.0,candidate_response_seconds=0.0)
        make=lambda subset:_mitc4_direct_layout(case,subset)
    rows=[];operator_seconds=0.;optimization_seconds=0.
    for subset in _subsets():
        t=time.perf_counter();operator=make(subset)
        operator_seconds+=time.perf_counter()-t
        if method=="mitc4_direct":
            record["counts"]["forward_analyses"]+=len(subset)+1
        t=time.perf_counter()
        cases=[optimize_controls(operator,target) for target in case["targets"]]
        optimization_seconds+=time.perf_counter()-t
        rows.append(dict(stations=list(subset),cases=cases))
    record.update(case_construction_seconds=case_construction_seconds,
                  operator_seconds=operator_seconds,
                  optimization_seconds=optimization_seconds,
                  total_seconds=time.perf_counter()-started,
                  peak_rss_kib=int(resource.getrusage(resource.RUSAGE_SELF).ru_maxrss),
                  layouts=rows)
    if method.startswith("mitc4"):
        record["structural_graph_nnz"]=_mitc4_graph_nnz(case)
    return record


def measure(repeats=3,meshes=((12,24),(16,32),(20,40)),methods=METHODS):
    summary=OUT/"aircraft_efficiency_four_methods.json"
    report=(json.loads(summary.read_text()) if summary.exists() else
        dict(scope="complete 31-layout, two-target linear screen",
             cpu="Intel Xeon Gold 6148; one BLAS/OpenMP thread",meshes={}))
    report["repeats"]=repeats
    report["methods"]=list(dict.fromkeys(list(report.get("methods",[]))+list(methods)))
    for nx,ny in meshes:
        key=f"{nx}x{ny}";report["meshes"].setdefault(key,{})
        for method in methods:
            records=[]
            for repeat in range(repeats):
                path=OUT/f"aircraft_efficiency_{key}_{method}_{repeat}.json"
                command=[sys.executable,"-m",
                         "experiments.engineering_inverse.efficiency_aircraft",
                         "--worker",str(nx),str(ny),method,str(path)]
                subprocess.run(command,check=True,env=os.environ.copy(),
                               stdout=subprocess.DEVNULL)
                row=json.loads(path.read_text());records.append(row)
                report["meshes"][key][method]=dict(
                    median_seconds=float(np.median([r["total_seconds"] for r in records])),
                    repeats=len(records),rows=[dict(total_seconds=r["total_seconds"],
                        peak_rss_kib=r["peak_rss_kib"],
                        base_seconds=r["base_seconds"],
                        candidate_response_seconds=r["candidate_response_seconds"],
                        operator_seconds=r["operator_seconds"],
                        optimization_seconds=r["optimization_seconds"])
                        for r in records])
                summary.write_text(
                    json.dumps(report,indent=2)+"\n")
                print(key,method,repeat,row["total_seconds"],flush=True)
    return report


if __name__=="__main__":
    parser=argparse.ArgumentParser()
    parser.add_argument("--worker",nargs=4)
    parser.add_argument("--repeats",type=int,default=3)
    parser.add_argument("--methods",nargs="+",choices=METHODS,default=list(METHODS))
    parser.add_argument("--meshes",nargs="+",default=["12x24","16x32","20x40"])
    args=parser.parse_args()
    if args.worker:
        nx,ny,method,path=args.worker
        Path(path).write_text(json.dumps(worker(int(nx),int(ny),method),indent=2)+"\n")
    else:measure(args.repeats,
                 [tuple(map(int,item.split("x"))) for item in args.meshes],
                 args.methods)
