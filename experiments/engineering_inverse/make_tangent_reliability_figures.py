"""Vector figures and flat tables for the final tangent/reliability audit."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
OUT=HERE/"results"
FIG=HERE/"figures"
FIG.mkdir(exist_ok=True)


def save(name):
    plt.tight_layout()
    for suffix in ("pdf","svg"):
        plt.savefig(FIG/f"{name}.{suffix}",bbox_inches="tight")
    plt.close()


def run():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,
                         "pdf.fonttype":42,"svg.fonttype":"none",
                         "axes.spines.top":False,"axes.spines.right":False})
    aircraft=json.loads((OUT/"aircraft_efficiency_four_methods.json").read_text())
    meshes=("12x24","16x32","20x40")
    methods=("cb_shared_sparse_legacy","cb_shared_sparse","mitc4_shared")
    labels=("CB original","CB cached","MITC4 shared")
    colors=("#8ca6b8","#0b4969","#ae543f")
    rows=[]
    for mesh in meshes:
        for method in methods:
            series=aircraft["meshes"][mesh][method]
            records=[json.loads((OUT/f"aircraft_efficiency_{mesh}_{method}_{i}.json").read_text())
                     for i in range(3)]
            med=lambda key:float(np.median([r[key] for r in records]))
            rows.append(dict(mesh=mesh,method=method,total_s=series["median_seconds"],
                             tangent_assembly_s=med("shell_tangent_seconds") if method.startswith("cb") else None,
                             factorization_s=med("base_factor_seconds") if method.startswith("cb") else None,
                             candidate_response_s=med("candidate_response_seconds"),
                             peak_rss_MiB=med("peak_rss_kib")/1024,
                             state_dof=records[0]["total_state_dofs"]))
    with (OUT/"tangent_aircraft_final.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=rows[0].keys())
        writer.writeheader();writer.writerows(rows)
    x=np.arange(3)
    fig,ax=plt.subplots(figsize=(6.7,3.2))
    for j,method in enumerate(methods[:2]):
        vals=[next(r["tangent_assembly_s"] for r in rows
                   if r["mesh"]==mesh and r["method"]==method) for mesh in meshes]
        ax.bar(x+(j-.5)*.3,vals,.3,label=labels[j],color=colors[j])
    ax.set_xticks(x,meshes)
    ax.set_ylabel("Exact shell tangent assembly (s)")
    ax.set_xlabel("Aircraft mesh")
    ax.legend(frameon=False)
    save("final_aircraft_tangent_assembly")
    fig,ax=plt.subplots(figsize=(6.7,3.2))
    for j,method in enumerate(methods):
        vals=[next(r["total_s"] for r in rows
                   if r["mesh"]==mesh and r["method"]==method) for mesh in meshes]
        ax.bar(x+(j-1)*.23,vals,.23,label=labels[j],color=colors[j])
    ax.set_xticks(x,meshes)
    ax.set_ylabel("Full 31-layout linear screen (s)")
    ax.set_xlabel("Aircraft mesh")
    ax.legend(frameon=False)
    save("final_aircraft_shared_screen")
    fig,ax=plt.subplots(figsize=(6.7,3.2))
    pair=[r for r in rows if r["mesh"]=="20x40"]
    contributions=np.array([
        [r["tangent_assembly_s"] or 0.,
         r["factorization_s"] or 0.,
         r["candidate_response_s"]] for r in pair])
    other=np.array([max(0.,r["total_s"]-sum(parts))
                    for r,parts in zip(pair,contributions)])
    contributions=np.column_stack((contributions,other))
    bottom=np.zeros(3)
    for j,(label,color) in enumerate(zip(("CB tangent","CB factorization",
                                          "Candidate responses","Other"),
                                         ("#0b4969","#8ca6b8","#ae543f","#cccccc"))):
        ax.bar(np.arange(3),contributions[:,j],bottom=bottom,
               color=color,label=label)
        bottom+=contributions[:,j]
    ax.set_xticks(np.arange(3),labels)
    ax.set_ylabel("20×40 full screen (s)")
    ax.legend(frameon=False,ncol=2)
    save("final_aircraft_cost_attribution")
    forward=json.loads((OUT/"automotive_tangent_before_after.json").read_text())
    frows=[]
    for strategy in ("legacy","cached"):
        for repeat,row in enumerate(forward[strategy]):
            frows.append(dict(strategy=strategy,repeat=repeat,
                              total_s=row["total_seconds"],
                              tangent_assembly_s=row["assembly_seconds"],
                              factorization_s=row["factor_seconds"],
                              newton_iterations=row["iterations"]))
    with (OUT/"tangent_automotive_forward_final.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=frows[0].keys())
        writer.writeheader();writer.writerows(frows)
    fig,ax=plt.subplots(figsize=(5.4,3.0))
    for j,strategy in enumerate(("legacy","cached")):
        r=[v for v in frows if v["strategy"]==strategy]
        assembly=np.median([v["tangent_assembly_s"] for v in r])
        total=np.median([v["total_s"] for v in r])
        ax.bar(j,assembly,color=colors[j],label="Tangent assembly" if j==0 else None)
        ax.bar(j,total-assembly,bottom=assembly,color="#c6ccd0",
               label="Other forward cost" if j==0 else None)
    ax.set_xticks((0,1),("Original","Cached"))
    ax.set_ylabel("Cold nonlinear forward (s)")
    ax.legend(frameon=False)
    save("final_automotive_forward_cost")
    matched=json.loads((OUT/"automotive_efficiency_forward.json").read_text())
    fig,ax=plt.subplots(figsize=(5.4,3.0))
    forward_times=[float(np.median([row["elapsed_seconds"]
                   for row in matched["methods"][method]]))
                   for method in ("cb","asd")]
    ax.bar(("CB","ASDShellQ4"),forward_times,color=("#0b4969","#ae543f"))
    ax.set_ylabel("Same final controls: cold forward (s)")
    save("final_automotive_same_input_forward")
    methods_auto=("cb_consistent","cb_fd","asd_fd")
    auto=[json.loads((OUT/f"automotive_efficiency_{m}_12x8_zero.json").read_text())
          for m in methods_auto]
    with (OUT/"tangent_automotive_optimization_final.csv").open("w",newline="") as file:
        writer=csv.writer(file)
        writer.writerow(("method","total_s","iterations","equilibria","objective",
                         "max_error_mm","max_force_N","reference_objective",
                         "reference_max_error_mm","reference_max_force_N",
                         "projected_kkt_gradient_per_mm"))
        references=json.loads((OUT/"automotive_efficiency_common_reference_zero.json").read_text())
        for a in auto:
            ref=references["methods"][a["method"]]
            stationarity=(json.loads((OUT/f"automotive_reliability_{a['method']}_zero.json").read_text())
                          ["kkt_projected_gradient_infinity_norm"]
                          if a["method"].startswith("cb") else
                          json.loads((OUT/"automotive_asd_stationarity.json").read_text())
                          ["gradients"]["0.01"]["projected_infinity_norm"])
            writer.writerow((a["method"],a["total_seconds"],a["iterations"],
                             a["counts"]["equilibrium_solves"],a["final"]["objective"],
                             a["final"]["max_abs_error_mm"],
                             a["final"]["max_mount_reaction_N"],
                             ref["reference_objective"],
                             ref["reference_quality"]["max_abs_error_mm"],
                             ref["reference_max_reaction_N"],stationarity))
    fig,ax=plt.subplots(figsize=(6.0,3.0))
    ax.bar(("CB consistent","CB FD","ASD FD"),[a["total_seconds"] for a in auto],
           color=("#0b4969","#8ca6b8","#ae543f"))
    ax.set_ylabel("Complete nonlinear optimization (s)")
    save("final_automotive_optimization_cost")
    fig,ax=plt.subplots(figsize=(6.0,3.0))
    labels_auto=("CB consistent","CB FD","ASD FD")
    for a,label,color in zip(auto,labels_auto,("#0b4969","#8ca6b8","#ae543f")):
        history=a["history"]
        ax.plot([row["elapsed_seconds"] for row in history],
                [row["objective"] for row in history],
                label=label,color=color)
    ax.set_xlabel("Elapsed time (s)")
    ax.set_ylabel("Design-grid objective (mm²)")
    ax.legend(frameon=False)
    save("final_automotive_convergence")
    before=json.loads((OUT/"automotive_efficiency_branch_probe_pre_gauge.json").read_text())
    after=json.loads((OUT/"automotive_efficiency_branch_probe.json").read_text())
    fig,ax=plt.subplots(figsize=(5.8,3.0))
    names=("Cold","ASD initial")
    for j,record in enumerate((before,after)):
        values=[v["maximum_edge_error_mm"] for v in record["states"].values()]
        ax.bar(np.arange(2)+(j-.5)*.3,values,.3,
               color=colors[j],label=("Ungauged" if j==0 else "Gauge fixed"))
    ax.axhline(.6,color="black",ls="--",lw=.8)
    ax.set_xticks(np.arange(2),names)
    ax.set_ylabel("Same controls: max edge error (mm)")
    ax.legend(frameon=False)
    save("final_automotive_gauge_diagnostic")


if __name__=="__main__":
    run()
