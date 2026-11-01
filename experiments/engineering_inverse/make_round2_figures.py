"""Vector figures for the mechanics-consistency follow-up (saved data only)."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

HERE=Path(__file__).resolve().parent
OUT=HERE/"results"
FIG=HERE/"figures"
FIG.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.spines.top":False,
                     "axes.spines.right":False,"savefig.bbox":"tight"})


def _save(fig,name):
    for ext in ("svg","pdf"):
        fig.savefig(FIG/f"{name}.{ext}")
    plt.close(fig)


def _read(name):return json.loads((OUT/name).read_text())


def automotive():
    reference=_read("automotive_spring_reference.json")
    hard=_read("automotive_fixed_input_diagnostics.json")
    mesh_keys=list(reference["meshes"])
    size=np.array([int(k.split("x")[0])*int(k.split("x")[1]) for k in mesh_keys])
    fig,axes=plt.subplots(1,2,figsize=(8.5,3.2))
    labels=[("checkerboard","Checkerboard","#1f77b4"),
            ("mitc4_linear","ShellMITC4 linear","#ff7f0e"),
            ("asd_corotational","ASDShellQ4 corotational","#2ca02c")]
    rows=[]
    for key,label,color in labels:
        results=[reference["meshes"][mesh][key] for mesh in mesh_keys]
        axes[0].plot(size,[r["max_abs_error_mm"] for r in results],"o-",label=label,color=color)
        axes[1].plot(size,[r["max_mount_reaction_N"] for r in results],"o-",label=label,color=color)
        for mesh,r in zip(mesh_keys,results):
            rows.append(dict(mesh=mesh,method=key,
                             max_edge_error_mm=r["max_abs_error_mm"],
                             gap_rms_mm=r["gap_error_rms_mm"],
                             flush_rms_mm=r["flush_rms_mm"],
                             max_mount_reaction_N=r["max_mount_reaction_N"],
                             status=r["status"]))
    axes[0].axhline(.6,color=".4",ls="--",lw=.8)
    axes[1].axhline(300,color=".4",ls="--",lw=.8)
    axes[0].set(xlabel="Shell cells",ylabel="Maximum edge error (mm)",
                title="Same optimized bracket adjustments")
    axes[1].set(xlabel="Shell cells",ylabel="Largest bracket resultant (N)",
                title="Physical-area bracket connection")
    axes[0].legend(frameon=False,fontsize=7)
    _save(fig,"round2_automotive_common_reference")
    with (OUT/"round2_automotive_common_reference.csv").open("w",newline="") as stream:
        writer=csv.DictWriter(stream,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)

    fig,axes=plt.subplots(1,2,figsize=(8.5,3.2))
    hard_keys=list(hard["meshes"])
    x=np.array([int(k.split("x")[0])*int(k.split("x")[1]) for k in hard_keys])
    for method,label,color in (("checkerboard","Checkerboard","#1f77b4"),
                               ("mitc4","ShellMITC4","#ff7f0e")):
        val=[]
        for mesh in hard_keys:
            r=hard["meshes"][mesh]["cases"]["one_z"][method]["mount_reactions_N"]
            val.append(np.max(np.linalg.norm(np.asarray(r),axis=1)))
        axes[0].plot(x,val,"o-",color=color,label=label)
    for method,label,color in (("checkerboard","Checkerboard","#1f77b4"),
                               ("mitc4","ShellMITC4","#ff7f0e")):
        val=[]
        for mesh in hard_keys:
            r=_read("automotive_spring_audit.json")["meshes"][mesh]["cases"]["one_z"][method]["reaction_N"]
            val.append(np.max(np.linalg.norm(np.asarray(r),axis=1)))
        axes[1].plot(x,val,"o-",color=color,label=label)
    axes[0].set(xlabel="Shell cells",ylabel="Largest mount reaction (N)",
                title="Hard multi-node patch; 0.1 mm input")
    axes[1].set(xlabel="Shell cells",ylabel="Largest mount reaction (N)",
                title="Area-integrated elastic bracket; 0.1 mm input")
    axes[0].legend(frameon=False)
    _save(fig,"round2_automotive_connection_diagnosis")

    design=_read("automotive_spring_design_12x8.json")
    controls=np.asarray(design["adjustment_mm"]).reshape(4,2)
    reactions=np.linalg.norm(np.asarray(design["final"]["mount_reaction_N"]),axis=1)
    fig,axes=plt.subplots(1,2,figsize=(7.8,3.1))
    x=np.arange(4)
    axes[0].bar(x-.16,controls[:,0],width=.32,label="y",color="#1f77b4")
    axes[0].bar(x+.16,controls[:,1],width=.32,label="z",color="#ff7f0e")
    axes[0].axhline(1.5,color=".4",ls="--",lw=.8)
    axes[0].axhline(-1.5,color=".4",ls="--",lw=.8)
    axes[0].set(xticks=x,xticklabels=[f"M{i+1}" for i in x],
                ylabel="Adjustment (mm)",title="Optimized bracket positions")
    axes[0].legend(frameon=False)
    axes[1].bar(x,reactions,color="#1f77b4")
    axes[1].axhline(300,color=".4",ls="--",lw=.8)
    axes[1].set(xticks=x,xticklabels=[f"M{i+1}" for i in x],
                ylabel="Resultant reaction (N)",title="Checkerboard nonlinear balance")
    _save(fig,"round2_automotive_controls_reactions")

    fig,ax=plt.subplots(figsize=(6.0,2.8))
    ax.add_patch(plt.Rectangle((.05,.30),.82,.15,facecolor=".82",edgecolor=".35"))
    ax.add_patch(plt.Rectangle((.40,.37),.13,.33,facecolor="#a6c8e2",edgecolor="#1f77b4"))
    ax.add_patch(plt.Rectangle((.39,.74),.15,.10,facecolor=".55",edgecolor=".3"))
    ax.plot([.455,.455],[.70,.74],color=".2",lw=2)
    ax.annotate("y adjustment",xy=(.455,.81),xytext=(.63,.88),
                arrowprops=dict(arrowstyle="->",color="#1f77b4"),color="#1f77b4")
    ax.annotate("z adjustment (out of plane)",xy=(.45,.53),xytext=(.59,.60),
                arrowprops=dict(arrowstyle="->",color="#ff7f0e"),color="#ff7f0e")
    ax.text(.05,.24,"sheet flange, 0.10 × 0.12 m integrated footprint",fontsize=8)
    ax.text(.05,.93,"schematic effective bracket; two orthogonal flexure leaves",fontsize=9)
    ax.set(xlim=(0,1),ylim=(.15,1),xticks=[],yticks=[])
    _save(fig,"round2_automotive_bracket_schematic")

    from experiments.engineering_inverse.automotive_panel import (
        build_case,membrane_principal_strain)
    case=build_case()
    q=np.asarray(design["final"]["displacement_m"])
    strain=np.asarray(membrane_principal_strain(case,q))
    fig,ax=plt.subplots(figsize=(6,3.1))
    image=ax.pcolormesh(np.linspace(0,.8,13),np.linspace(0,.45,9),
                        strain*1e6,shading="flat",cmap="viridis")
    fig.colorbar(image,ax=ax,label="Largest absolute membrane principal strain (µε)")
    ax.set(xlabel="Material x (m)",ylabel="Material y (m)",
           title="Optimized panel; bending stress excluded")
    _save(fig,"round2_automotive_membrane_strain")

    finest=reference["meshes"][mesh_keys[-1]]
    coordinate=np.linspace(0,.8,len(finest["checkerboard"]["gap_mm"]))
    fig,axes=plt.subplots(2,1,figsize=(7.3,4.7),sharex=True)
    for key,label,color in labels:
        row=finest[key]
        axes[0].plot(coordinate,np.asarray(row["gap_mm"])-3,"o-",ms=2,
                     color=color,label=label)
        axes[1].plot(coordinate,row["flush_mm"],"o-",ms=2,color=color)
    axes[0].axhline(0,color=".4",lw=.7);axes[1].axhline(0,color=".4",lw=.7)
    axes[0].set(ylabel="Gap error (mm)",title=f"Inspection edge, {mesh_keys[-1]} grid")
    axes[1].set(xlabel="Edge x (m)",ylabel="Flush (mm)")
    axes[0].legend(frameon=False,fontsize=7)
    _save(fig,"round2_automotive_edge_distribution")
    cross=OUT/"automotive_design_crosscheck.json"
    if cross.exists():
        data=json.loads(cross.read_text())
        if "24x20" in data["meshes"]:
            mesh=data["meshes"]["24x20"]
            designs=list(mesh)
            methods=[("checkerboard","Checkerboard","#1f77b4"),
                     ("asd_corotational","ASDShellQ4 corotational","#2ca02c")]
            fig,axes=plt.subplots(1,2,figsize=(8.3,3.2))
            x=np.arange(len(designs))
            for k,(method,label,color) in enumerate(methods):
                dx=(k-.5)*.3
                axes[0].bar(x+dx,[mesh[d][method]["max_abs_error_mm"] for d in designs],
                            width=.3,color=color,label=label)
                axes[1].bar(x+dx,[mesh[d][method]["max_mount_reaction_N"] for d in designs],
                            width=.3,color=color)
            axes[0].axhline(.6,color=".4",ls="--",lw=.8)
            axes[1].axhline(300,color=".4",ls="--",lw=.8)
            labels=["CB nonlinear\ndesign","ASD seeded\nrefinement"]
            for ax in axes:
                ax.set(xticks=x,xticklabels=labels)
            axes[0].set(ylabel="Maximum edge error (mm)",
                        title="Both designs at 24×20")
            axes[1].set(ylabel="Largest bracket reaction (N)",
                        title="Same physical connector")
            axes[0].legend(frameon=False,fontsize=7)
            _save(fig,"round2_automotive_two_designs")


def aircraft():
    reference=OUT/"aircraft_corotational_reference.json"
    if not reference.exists():return
    report=json.loads(reference.read_text())
    meshes=list(report["meshes"])
    if not meshes:return
    fig,axes=plt.subplots(1,2,figsize=(8,3.2))
    for target_id in (0,1):
        rows=[report["meshes"][mesh][target_id] for mesh in meshes]
        x=np.arange(len(meshes))
        axes[0].plot(x,[r["rms_error_mm"] for r in rows],"o-",label=f"Case {target_id+1}")
        axes[1].plot(x,[r["max_error_mm"] for r in rows],"o-",label=f"Case {target_id+1}")
    for ax,limit,label in ((axes[0],.25,"RMS error (mm)"),
                           (axes[1],.8,"Maximum error (mm)")):
        ax.axhline(limit,color=".4",ls="--",lw=.8)
        ax.set(xticks=np.arange(len(meshes)),xticklabels=meshes,
               xlabel="Shell mesh",ylabel=label)
    axes[0].legend(frameon=False)
    _save(fig,"round2_aircraft_corotational_reference")
    perf=OUT/"performance_round2.json"
    if perf.exists():
        p=json.loads(perf.read_text())
        fig,ax=plt.subplots(figsize=(5.4,3.2))
        x=np.arange(len(p["meshes"]))
        ax.plot(x,[p["meshes"][m]["checkerboard_median_seconds"] for m in p["meshes"]],
                "o-",label="Checkerboard shared factor")
        ax.plot(x,[p["meshes"][m]["mitc4_median_seconds"] for m in p["meshes"]],
                "s-",label="ShellMITC4 rebuild workflow")
        ax.set(xticks=x,xticklabels=list(p["meshes"]),xlabel="Shell mesh",
               ylabel="Complete 31-layout linear screen (s)")
        ax.legend(frameon=False)
        _save(fig,"round2_aircraft_linear_screen_cost")
    audit=OUT/"aircraft_nonlinear_audit_12x24.json"
    if audit.exists():
        entries=json.loads(audit.read_text())["statuses"]
        table=[]
        fig,ax=plt.subplots(figsize=(6.5,3.5))
        legend_done=set()
        for entry in entries:
            cases=[c for c in entry["cases"] if c is not None
                   and c.get("rms_error_mm") is not None]
            if not cases:continue
            complete=len(cases)==2
            table.append(dict(stations="-".join(map(str,entry["stations"])),
                actuator_count=entry["count"],status=entry["status"],
                both_cases_evaluated=complete,
                worst_observed_rms_mm=max(c["rms_error_mm"] for c in cases),
                worst_observed_max_mm=max(c["max_error_mm"] for c in cases)))
            label=("Both cases evaluated" if complete else
                   "One case evaluated; shown value is a lower bound")
            ax.scatter(entry["count"],max(c["rms_error_mm"] for c in cases),
                       marker="o" if complete else "^",s=45,
                       facecolor="#1f77b4" if complete else "white",
                       edgecolor="#1f77b4",
                       label=label if label not in legend_done else None)
            legend_done.add(label)
        final=json.loads((OUT/"aircraft_nonlinear_final_12x24.json").read_text())
        ax.scatter(3,max(c["rms_error_mm"] for c in final["cases"]),
                   marker="*",s=105,color="#2ca02c",label="[0,2,4] verified")
        ax.axhline(.25,color=".4",ls="--",lw=.8,label="RMS study tolerance")
        ax.set(xlim=(.5,3.5),xticks=[1,2,3],xlabel="Actuator count",
               ylabel="Worst observed case RMS error (mm)",
               title="Nonlinear lower-cardinality audit; local optimizations")
        ax.legend(frameon=False,fontsize=7)
        _save(fig,"round2_aircraft_nonlinear_layout_audit")
        if table:
            with (OUT/"round2_aircraft_nonlinear_layout_audit.csv").open("w",newline="") as stream:
                writer=csv.DictWriter(stream,fieldnames=list(table[0]))
                writer.writeheader();writer.writerows(table)


if __name__=="__main__":automotive();aircraft()
