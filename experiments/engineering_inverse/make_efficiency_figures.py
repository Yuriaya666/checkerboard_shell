"""Rebuild paper-ready efficiency figures and flat CSV tables from raw JSON."""
from __future__ import annotations

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

HERE=Path(__file__).resolve().parent
RESULTS=HERE/"results"
FIGURES=HERE/"figures"
FIGURES.mkdir(exist_ok=True)


def save(name):
    plt.tight_layout()
    for suffix in ("svg","pdf"):
        plt.savefig(FIGURES/f"{name}.{suffix}",bbox_inches="tight")
    plt.close()


def run():
    plt.rcParams.update({"font.family":"DejaVu Sans","font.size":9,
        "axes.spines.top":False,"axes.spines.right":False,
        "pdf.fonttype":42,"svg.fonttype":"none"})
    aircraft=json.loads((RESULTS/"aircraft_efficiency_four_methods.json").read_text())
    meshes=list(aircraft["meshes"])
    methods=[m for m in ("cb_direct","cb_shared","cb_shared_sparse",
                         "mitc4_direct","mitc4_shared")
             if all(m in aircraft["meshes"][grid] for grid in meshes)]
    labels={"cb_direct":"CB direct","cb_shared":"CB shared, dense",
            "cb_shared_sparse":"CB shared, sparse+rank 1",
            "mitc4_direct":"MITC4 direct","mitc4_shared":"MITC4 shared"}
    colors={"cb_direct":"#8ca6b8","cb_shared":"#2c6c8c",
            "cb_shared_sparse":"#0b4969","mitc4_direct":"#d29877",
            "mitc4_shared":"#ae543f"}
    fig,ax=plt.subplots(figsize=(7.2,3.5))
    x=np.arange(len(meshes));width=.16
    for j,method in enumerate(methods):
        values=[aircraft["meshes"][mesh][method]["median_seconds"] for mesh in meshes]
        ax.bar(x+(j-(len(methods)-1)/2)*width,values,width,label=labels[method],
               color=colors[method])
    ax.set_xticks(x,meshes);ax.set_xlabel("Aircraft mesh (chord × span cells)")
    ax.set_ylabel("Full 31-layout linear screen (s)")
    ax.legend(ncol=2,frameon=False,fontsize=8)
    save("efficiency_aircraft_four_strategies")

    fig,ax=plt.subplots(figsize=(6.5,3.1))
    cb_reuse=[];mt_reuse=[];shared_ratio=[]
    for mesh in meshes:
        d=aircraft["meshes"][mesh]
        cb_reuse.append(d["cb_direct"]["median_seconds"]/
                        d["cb_shared_sparse"]["median_seconds"])
        mt_reuse.append(d["mitc4_direct"]["median_seconds"]/
                        d["mitc4_shared"]["median_seconds"])
        shared_ratio.append(d["cb_shared_sparse"]["median_seconds"]/
                            d["mitc4_shared"]["median_seconds"])
    ax.plot(meshes,cb_reuse,"o-",color="#0b4969",label="CB reuse gain")
    ax.plot(meshes,mt_reuse,"s-",color="#ae543f",label="MITC4 reuse gain")
    ax.plot(meshes,shared_ratio,"^-",color="#627469",
            label="MITC4 speedup over CB (shared)")
    ax.set_ylabel("Measured time ratio (×)")
    ax.set_xlabel("Aircraft mesh (chord × span cells)")
    ax.legend(frameon=False,fontsize=8)
    save("efficiency_aircraft_reuse_ratios")

    fig,ax=plt.subplots(figsize=(7.2,3.5))
    pair=("cb_shared_sparse","mitc4_shared")
    names=[];data=[]
    for mesh in meshes:
        for method in pair:
            series=aircraft["meshes"][mesh][method]
            record=json.loads((RESULTS/f"aircraft_efficiency_{mesh}_{method}_0.json").read_text())
            rows=series["rows"]
            med=lambda field:float(np.median([row[field] for row in rows]))
            if method.startswith("cb"):
                parts=[record["shell_tangent_seconds"],
                       record["base_system_build_seconds"],
                       record["base_factor_seconds"],
                       med("candidate_response_seconds"),
                       med("operator_seconds")+med("optimization_seconds")]
            else:
                parts=[record["model_construction_seconds"],
                       record["base_factor_and_solve_seconds"],0.,
                       med("candidate_response_seconds"),
                       med("operator_seconds")+med("optimization_seconds")]
            parts.append(max(0.,series["median_seconds"]-sum(parts)))
            data.append(parts)
            names.append(f"{mesh}\n{'CB' if method.startswith('cb') else 'MITC4'}")
    bottom=np.zeros(len(data));data=np.asarray(data)
    for j,(label,color) in enumerate(zip(("CB tangent / MITC4 model build",
            "CB beam build / MITC4 first solve","CB factorization",
            "Candidate responses","Layout + optimizer","Other"),
            ("#336b89","#7ea4b9","#244b68","#c17b61","#d9ab8d","#b9bdc0"))):
        ax.bar(np.arange(len(data)),data[:,j],bottom=bottom,color=color,label=label)
        bottom+=data[:,j]
    ax.set_xticks(np.arange(len(data)),names)
    ax.set_ylabel("Full 31-layout screen (s)")
    ax.legend(ncol=2,frameon=False,fontsize=7)
    save("efficiency_aircraft_cost_breakdown")

    with (RESULTS/"efficiency_aircraft_summary.csv").open("w",newline="") as handle:
        writer=csv.writer(handle)
        writer.writerow(["mesh","method","median_total_s","median_peak_rss_MiB",
            "total_state_dof","free_state_dof","numeric_local_nnz_or_none",
            "structural_graph_nnz_or_none","factorization_s_or_none",
            "candidate_response_s","layout_operator_s","optimization_s"])
        for mesh in meshes:
            for method in methods:
                series=aircraft["meshes"][mesh][method]
                record=json.loads((RESULTS/f"aircraft_efficiency_{mesh}_{method}_0.json").read_text())
                rows=series["rows"]
                med=lambda field:float(np.median([row[field] for row in rows]))
                writer.writerow([mesh,method,series["median_seconds"],
                    med("peak_rss_kib")/1024,record["total_state_dofs"],
                    record["free_state_dofs"],record.get("full_local_nnz"),
                    record.get("structural_graph_nnz"),record.get("base_factor_seconds"),
                    med("candidate_response_seconds"),med("operator_seconds"),
                    med("optimization_seconds")])

    auto=[]
    for method in ("cb_consistent","cb_fd","asd_fd"):
        path=RESULTS/f"automotive_efficiency_{method}_12x8_zero.json"
        if path.exists():auto.append(json.loads(path.read_text()))
    if len(auto)==3:
        fig,(ax1,ax2)=plt.subplots(1,2,figsize=(7.2,3.1))
        names=[{"cb_consistent":"CB consistent","cb_fd":"CB central FD",
                "asd_fd":"ASD central FD"}[a["method"]] for a in auto]
        bars=ax1.bar(names,[a["total_seconds"] for a in auto],
                     color=["#0b4969","#7199ad","#ae543f"])
        ax1.set_ylabel("Complete nonlinear optimization (s)")
        ax1.tick_params(axis="x",rotation=20)
        ax2.bar(names,[a["counts"]["equilibrium_solves"] for a in auto],
                color=["#0b4969","#7199ad","#ae543f"])
        ax2.set_ylabel("Equilibrium analyses, including FD probes")
        ax2.tick_params(axis="x",rotation=20)
        save("efficiency_automotive_complete")
        fig,ax=plt.subplots(figsize=(6.7,3.2))
        direct=[];fd=[];implicit=[];other=[]
        for a in auto:
            t=a["times"]
            fd.append(t["finite_difference_seconds"])
            implicit.append(t["implicit_gradient_seconds"])
            direct.append(max(0.,t["forward_seconds"]-t["finite_difference_seconds"]))
            other.append(max(0.,a["total_seconds"]-direct[-1]-fd[-1]-implicit[-1]))
        bottom=np.zeros(3)
        for values,label,color in ((direct,"Direct forward analyses","#2c6c8c"),
                                    (fd,"Centered FD probes","#a9bac4"),
                                    (implicit,"Consistent derivative","#c17b61"),
                                    (other,"Optimizer + overhead","#d9ab8d")):
            ax.bar(names,values,bottom=bottom,label=label,color=color)
            bottom+=np.asarray(values)
        ax.set_ylabel("Complete nonlinear optimization (s)")
        ax.tick_params(axis="x",rotation=20)
        ax.legend(frameon=False,fontsize=7)
        save("efficiency_automotive_cost_breakdown")
        fig,ax=plt.subplots(figsize=(6.7,3.2))
        for a in auto:
            h=a["history"]
            ax.plot([row["elapsed_seconds"] for row in h],
                    [row["objective"] for row in h],label=a["method"])
        ax.set_xlabel("Elapsed time (s)");ax.set_ylabel("Weighted objective (mm²)")
        ax.legend(frameon=False)
        save("efficiency_automotive_convergence")
        with (RESULTS/"efficiency_automotive_summary.csv").open("w",newline="") as handle:
            writer=csv.writer(handle)
            writer.writerow(["start","method","mesh","optimizer_success","engineering_feasible",
                "iterations","equilibria","gradients","total_s","forward_s",
                "gradient_s","state_dof","objective","max_error_mm","max_reaction_N"])
            rows=auto[:]
            for method in ("cb_consistent","cb_fd","asd_fd"):
                path=RESULTS/f"automotive_efficiency_{method}_12x8_cb_feasible.json"
                if path.exists():rows.append(json.loads(path.read_text()))
            for a in rows:
                writer.writerow([a["start"],a["method"],"12x8",a["optimizer_success"],
                    a["final"]["engineering_feasible"],a["iterations"],
                    a["counts"]["equilibrium_solves"],
                    a["counts"]["linearizations"]+a["counts"]["finite_difference_gradients"],
                    a["total_seconds"],a["times"]["forward_seconds"],
                    a["times"]["implicit_gradient_seconds"]+
                    a["times"]["finite_difference_seconds"],a["state_dofs"],
                    a["final"]["objective"],a["final"]["max_abs_error_mm"],
                    a["final"]["max_mount_reaction_N"]])

        common=json.loads((RESULTS/"automotive_efficiency_common_reference_zero.json").read_text())
        fig,(a1,a2,a3)=plt.subplots(1,3,figsize=(8.5,2.9))
        methods=("cb_consistent","cb_fd","asd_fd")
        color=["#0b4969","#7199ad","#ae543f"]
        refs=[common["methods"][m] for m in methods]
        for ax,values,limit,ylabel in (
            (a1,[r["reference_quality"]["max_abs_error_mm"] for r in refs],.6,
             "Reference max edge error (mm)"),
            (a2,[r["reference_max_reaction_N"] for r in refs],300,
             "Reference max bracket force (N)"),
            (a3,[r["reference_objective"] for r in refs],None,
             "Reference weighted objective (mm²)")):
            ax.bar(names,values,color=color)
            if limit is not None:
                ax.axhline(limit,color="black",lw=.9,ls="--")
                ax.set_ylim(0,1.1*limit)
            ax.set_ylabel(ylabel);ax.tick_params(axis="x",rotation=30)
        save("efficiency_automotive_common_quality")

        with (RESULTS/"efficiency_integrated_summary.csv").open("w",newline="") as handle:
            writer=csv.writer(handle)
            writer.writerow(["task","method","quality_scope","reference_accepted",
                "state_dof","total_seconds","common_acceptance_seconds",
                "design_plus_acceptance_seconds","speed_ratio_to_reference_method",
                "speed_ratio_reference","reference_max_error_mm",
                "reference_max_reaction_N","reference_objective"])
            for mesh in meshes:
                d=aircraft["meshes"][mesh]
                for method in ("cb_direct","cb_shared","cb_shared_sparse",
                               "mitc4_direct","mitc4_shared"):
                    record=json.loads((RESULTS/f"aircraft_efficiency_{mesh}_{method}_0.json").read_text())
                    writer.writerow([f"aircraft_31_layout_linear_{mesh}",method,
                        "screening only; saved [0,2,4] separately nonlinear-checked",
                        "not applicable to every screened layout",record["total_state_dofs"],
                        d[method]["median_seconds"],"","",
                        d["mitc4_shared"]["median_seconds"]/
                        d[method]["median_seconds"],"mitc4_shared",
                        "","",""])
            asd=next(a for a in auto if a["method"]=="asd_fd")
            for a in auto:
                ref=common["methods"][a["method"]]
                writer.writerow(["automotive_complete_nonlinear",a["method"],
                    "common ASDShellQ4 24x20 engineering limits",
                    ref["reference_accepted"],a["state_dofs"],a["total_seconds"],
                    ref["reference_seconds"],a["total_seconds"]+ref["reference_seconds"],
                    asd["total_seconds"]/a["total_seconds"],"asd_fd",
                    ref["reference_quality"]["max_abs_error_mm"],
                    ref["reference_max_reaction_N"],ref["reference_objective"]])


if __name__=="__main__":run()
