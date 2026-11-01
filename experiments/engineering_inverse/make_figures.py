"""Rebuild compact vector figures and CSV summaries from saved results."""

import csv
import json
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from experiments.engineering_inverse.aircraft_skin import (
    build_case as aircraft_case, PARAMS as AIRCRAFT_PARAMS)
from experiments.engineering_inverse.automotive_panel import (
    build_case as automotive_case, membrane_principal_strain)

HERE=Path(__file__).resolve().parent
RESULTS=HERE/"results"
FIGURES=HERE/"figures"
FIGURES.mkdir(exist_ok=True)
plt.rcParams.update({"font.size":9,"axes.spines.top":False,
                     "axes.spines.right":False,"savefig.bbox":"tight"})


def save(fig,name):
    for ext in ("pdf","svg"):
        fig.savefig(FIGURES/f"{name}.{ext}")
    plt.close(fig)


def aircraft():
    data=json.loads((RESULTS/"aircraft_12x24.json").read_text())
    final=json.loads((RESULTS/"aircraft_nonlinear_final_12x24.json").read_text())
    case=aircraft_case()
    xyz=case["mesh"].reference_vertices.numpy()
    fig,ax=plt.subplots(figsize=(3.6,6.0))
    for j in range(0,25,4): ax.plot(xyz[j,:,0]+.42,xyz[j,:,1],color=".75",lw=.6)
    for i in range(0,13,3): ax.plot(xyz[:,i,0]+.42,xyz[:,i,1],color=".75",lw=.6)
    for line in case["beam_lines"]:
        ids=np.asarray(line[0]);pts=xyz.reshape(-1,3)[ids]
        ax.plot(pts[:,0]+.42,pts[:,1],color=".25",lw=1.3)
    for k,patch in enumerate(case["patches"]):
        p=xyz.reshape(-1,3)[patch].mean(axis=0)
        ax.scatter(p[0]+.42,p[1],color="tab:blue" if k in final["stations"] else "white",
                   edgecolor="tab:blue",s=60,zorder=5)
        ax.text(p[0]+.425,p[1],str(k+1))
    ax.set(xlabel="Chord x (m)",ylabel="Span y (m)",
           title="Upper skin, bonded beam lines and candidate stroke attachments")
    ax.set_aspect("equal");save(fig,"aircraft_structure")
    fig,ax=plt.subplots(figsize=(6,3.3))
    for method,color in (("checkerboard","tab:blue"),("mitc4","tab:orange")):
        rows=data[method]["search"]
        best=[min(r["worst_rms_mm"] for r in rows if r["count"]==n)
              for n in range(1,6)]
        ax.plot(range(1,6),best,"o-",label=method,color=color)
    ax.axhline(.25,color=".3",ls="--",label="linear RMS tolerance")
    ax.set(xlabel="Actuator count",ylabel="Best worst-case RMS error (mm)",
           title="Complete 31-layout linear screening")
    ax.legend(frameon=False);save(fig,"aircraft_layout_curve")
    fig,axes=plt.subplots(1,2,figsize=(9,3.3),sharey=True)
    for k,ax in enumerate(axes):
        q=np.asarray(final["cases"][k]["displacement_m"]).reshape(xyz.shape)
        j=12
        x=xyz[j,:,0]+.42
        ax.plot(x,xyz[j,:,2]*1e3,label="reference",color=".5")
        xi=np.linspace(0,1,13)
        target=(xyz[j,:,2]+AIRCRAFT_PARAMS["uniform_trailing_edge_increment_m"]
                *(1.5*xi**2-.5*xi**3))
        ax.plot(x,target*1e3,label="target",color="black",ls="--")
        ax.plot(x,(xyz[j,:,2]+q[j,:,2])*1e3,label="optimized",color="tab:blue")
        ax.set(xlabel="Chord x (m)",title=f"Case {k+1}, midspan; true deformation")
    axes[0].set_ylabel("Upper surface z (mm)")
    axes[0].legend(frameon=False);save(fig,"aircraft_profiles")
    fig,ax=plt.subplots(figsize=(6,3.1))
    q=np.asarray(final["cases"][1]["displacement_m"]).reshape(xyz.shape)
    y=xyz[:,-1,1]
    target=xyz[:,-1,2]+case["targets"][1].reshape(25,11)[:,-1]
    ax.plot(y,xyz[:,-1,2]*1e3,label="reference",color=".5")
    ax.plot(y,target*1e3,label="target",color="black",ls="--")
    ax.plot(y,(xyz[:,-1,2]+q[:,-1,2])*1e3,
            label="optimized",color="tab:blue")
    ax.set(xlabel="Span y (m)",ylabel="Trailing-edge z (mm)",
           title="Case 2 spanwise camber gradient; true deformation")
    ax.legend(frameon=False);save(fig,"aircraft_trailing_edge_span")
    fig=plt.figure(figsize=(6,4));ax=fig.add_subplot(111,projection="3d")
    q=np.asarray(final["cases"][1]["displacement_m"]).reshape(xyz.shape)
    ax.plot_surface(xyz[...,0]+.42,xyz[...,1],(xyz[...,2]+q[...,2])*1e3,
                    color="steelblue",alpha=.85,edgecolor=".5",lw=.2)
    ax.set(xlabel="Chord x (m)",ylabel="Span y (m)",zlabel="z (mm)",
           title="Case 2: actual scale deformation")
    save(fig,"aircraft_deformed_3d")
    fig,ax=plt.subplots(figsize=(6,3.2))
    strain=np.asarray(membrane_principal_strain(
        {"manufactured":xyz},np.asarray(final["cases"][1]["displacement_m"])))
    image=ax.pcolormesh(np.linspace(.42,.60,13),np.linspace(0,1,25),
                        strain*1e6,shading="flat",cmap="viridis")
    fig.colorbar(image,ax=ax,label="Max absolute membrane principal strain (µε)")
    ax.set(xlabel="Chord x (m)",ylabel="Span y (m)",
           title="Case 2 membrane strain; excludes surface bending strain")
    save(fig,"aircraft_membrane_strain")
    timing=json.loads((RESULTS/"timing_repeats.json").read_text())
    fig,ax=plt.subplots(figsize=(5,3))
    sizes=["12×24","16×32","20×40"]
    ac=[timing["median"]["aircraft"]["checkerboard_seconds"]]
    mt=[timing["median"]["aircraft"]["mitc4_seconds"]]
    for size in ("16x32","20x40"):
        result=json.loads((RESULTS/f"aircraft_{size}.json").read_text())
        ac.append(result["checkerboard"]["total_seconds"])
        mt.append(result["mitc4"]["search_seconds"])
    x=np.arange(3)
    ax.plot(x,ac,"o-",color="tab:blue",label="Checkerboard")
    ax.plot(x,mt,"s-",color="tab:orange",label="MITC4")
    ax.set(xticks=x,xticklabels=sizes,ylabel="Full 31-layout wall time (s)",
           title="Linear screening cost; excludes nonlinear audit")
    ax.legend(frameon=False)
    save(fig,"aircraft_screen_cost")
    workflows=json.loads((RESULTS/"three_workflows.json").read_text())["methods"]
    fig,axes=plt.subplots(1,2,figsize=(8,3.1))
    names=["CB exact","CB FD","MITC4 FD"]
    colors=["tab:blue","steelblue","tab:orange"]
    axes[0].bar(names,[r["design_seconds"] for r in workflows],color=colors)
    axes[0].set(ylabel="Inner design time (s)",title="Same three-station twist case")
    axes[1].bar(names,[r["forward_solves"] for r in workflows],color=colors)
    axes[1].set(ylabel="Forward solves or analyses",title="Linear workflow work")
    save(fig,"aircraft_three_workflows")
    rows=[]
    for method in ("checkerboard","mitc4"):
        for r in data[method]["search"]:
            rows.append(dict(method=method,stations="-".join(map(str,r["stations"])),
                             count=r["count"],worst_rms_mm=r["worst_rms_mm"],
                             worst_max_mm=r["worst_max_mm"],linear_feasible=r["feasible"]))
    write_csv(RESULTS/"aircraft_layout_summary.csv",rows)


def automotive():
    data=json.loads((RESULTS/"automotive_12x8.json").read_text())
    reference=json.loads((RESULTS/"automotive_reference.json").read_text())
    case=automotive_case();nominal=case["nominal"];manufactured=case["manufactured"]
    fig=plt.figure(figsize=(6,4));ax=fig.add_subplot(111,projection="3d")
    ax.plot_surface(nominal[...,0],nominal[...,1],nominal[...,2]*1e3,
                    alpha=.65,color="silver",edgecolor=".55",lw=.25)
    m=manufactured.reshape(-1,3)[case["mounts"]]
    for patch in case["mount_patches"]:
        p=manufactured.reshape(-1,3)[patch]
        ax.scatter(p[:,0],p[:,1],p[:,2]*1e3,color="tab:blue",
                   alpha=.35,s=14)
    ax.scatter(m[:,0],m[:,1],m[:,2]*1e3,color="tab:blue",s=35,label="adjustable mounts")
    ax.plot(nominal[-1,:,0],nominal[-1,:,1],nominal[-1,:,2]*1e3,
            color="tab:orange",lw=2,label="nominal inspection edge")
    ax.plot(manufactured[-1,:,0],manufactured[-1,:,1],
            manufactured[-1,:,2]*1e3,color="tab:red",lw=1.4,
            label="manufactured inspection edge")
    ax.set(xlabel="x (m)",ylabel="y (m)",zlabel="z (mm)",
           title="Representative fender patch and mount locations")
    ax.legend(frameon=False);save(fig,"automotive_structure")
    x=np.linspace(0,case["width_m"],13)
    fig,axes=plt.subplots(2,1,figsize=(7,5),sharex=True)
    before=data["checkerboard"]["initial"]
    after=data["checkerboard_nonlinear_feasible"]["final"]
    candidate=json.loads((RESULTS/"automotive_crosscheck.json").read_text())["checkerboard_nonlinear"]
    for ax,key,label in zip(axes,("gap_mm","flush_mm"),("Gap (mm)","Flush (mm)")):
        ax.plot(x,before[key],"o--",label="manufactured",color=".5")
        ax.plot(x,after[key],"s-",label="CB own optimum",color="tab:blue")
        ax.plot(x,candidate[key],"^-",label="cross-feasible candidate",color="tab:orange")
        if key=="gap_mm":ax.axhline(3,color="black",lw=.7)
        else:ax.axhline(0,color="black",lw=.7)
        ax.set_ylabel(label)
    axes[-1].set_xlabel("Inspection-edge x (m)")
    axes[0].legend(frameon=False,ncol=2)
    save(fig,"automotive_gap_flush")
    fig,ax=plt.subplots(figsize=(6,3.2))
    strain=np.asarray(membrane_principal_strain(case,
        np.asarray(data["checkerboard_nonlinear_feasible"]["displacement_m"])))
    image=ax.pcolormesh(np.linspace(0,case["width_m"],13),
                        np.linspace(0,case["height_m"],9),
                        strain*1e6,shading="flat",cmap="viridis")
    fig.colorbar(image,ax=ax,label="Max absolute membrane principal strain (µε)")
    ax.set(xlabel="Material x (m)",ylabel="Material y (m)",
           title="Nonlinear adjusted panel; excludes surface bending strain")
    save(fig,"automotive_membrane_strain")
    fig,axes=plt.subplots(1,2,figsize=(8,3.2))
    controls=np.asarray(data["checkerboard_nonlinear_feasible"]["adjustment_mm"]).reshape(-1,2)
    labels=[f"M{i+1}" for i in range(len(controls))]
    axes[0].bar(np.arange(len(labels))-.16,controls[:,0],width=.32,label="y")
    axes[0].bar(np.arange(len(labels))+.16,controls[:,1],width=.32,label="z")
    axes[0].set(xticks=np.arange(len(labels)),xticklabels=labels,
                ylabel="Mount adjustment (mm)")
    axes[0].legend(frameon=False)
    reactions=np.asarray(data["checkerboard_nonlinear_feasible"]["mount_reaction_N"]).reshape(-1,3)
    axes[1].bar(labels,np.linalg.norm(reactions,axis=1))
    axes[1].set(ylabel="Nonlinear bracket resultant (N)")
    save(fig,"automotive_adjustments_reactions")
    fig,axes=plt.subplots(1,2,figsize=(9,3.3))
    names=["CB own\nnonlinear","CB → MITC4\nfine mesh",
           "MITC4 → MITC4\nfine mesh","Rigid datum"]
    values=[data["checkerboard_nonlinear_feasible"]["final"]["max_abs_error_mm"],
        reference["meshes"]["20x16"]["checkerboard"]["max_abs_error_mm"],
        reference["meshes"]["20x16"]["mitc4"]["max_abs_error_mm"],
        data["rigid_datum_baseline"]["final"]["max_abs_error_mm"]]
    axes[0].bar(names,values,color=["tab:blue","tab:blue","tab:orange",".6"])
    axes[0].axhline(.6,color="black",ls="--",label="research tolerance")
    axes[0].set(ylabel="Maximum edge error (mm)",title="Shape criterion")
    reactions=[data["checkerboard_nonlinear_feasible"]["maximum_mount_reaction_N"],
        reference["meshes"]["20x16"]["checkerboard"]["maximum_mount_reaction_N"],
        reference["meshes"]["20x16"]["mitc4"]["maximum_mount_reaction_N"]]
    axes[1].bar(names[:3],reactions,color=["tab:blue","tab:blue","tab:orange"])
    axes[1].axhline(300,color="black",ls="--",label="research limit")
    axes[1].set(ylabel="Largest bracket resultant (N)",title="Reaction criterion")
    axes[0].legend(frameon=False);save(fig,"automotive_cross_model")
    timing=json.loads((RESULTS/"timing_repeats.json").read_text())
    fig,ax=plt.subplots(figsize=(5,3))
    ac=[timing["median"]["automotive"]["checkerboard_seconds"]]
    mt=[timing["median"]["automotive"]["mitc4_seconds"]]
    for size in ("16x12","20x16"):
        result=json.loads((RESULTS/f"automotive_{size}.json").read_text())
        ac.append(result["checkerboard"]["total_seconds"])
        mt.append(result["mitc4"]["total_seconds"])
    x=np.arange(3)
    ax.plot(x,ac,"o-",color="tab:blue",label="Checkerboard")
    ax.plot(x,mt,"s-",color="tab:orange",label="MITC4")
    ax.set(xticks=x,xticklabels=["12×8","16×12","20×16"],
           ylabel="Linear design wall time (s)",
           title="Different implementations; no acceleration claim")
    ax.legend(frameon=False)
    save(fig,"automotive_cost")
    rows=[]
    for name,entry in reference["meshes"].items():
        for method,values in entry.items():
            rows.append(dict(mesh=name,design=method,
                max_error_mm=values["max_abs_error_mm"],
                gap_rms_mm=values["gap_error_rms_mm"],
                flush_rms_mm=values["flush_rms_mm"],
                max_reaction_N=values["maximum_mount_reaction_N"]))
    write_csv(RESULTS/"automotive_reference_summary.csv",rows)


def write_csv(path,rows):
    with path.open("w",newline="") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]))
        writer.writeheader();writer.writerows(rows)


if __name__=="__main__":
    aircraft();automotive()
