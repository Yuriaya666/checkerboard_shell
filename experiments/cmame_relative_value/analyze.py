"""Create derived CSV/JSON tables and vector figures from measured results."""
from __future__ import annotations

import csv
import json
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from .run_scan import HERE, DATA, save_rows, THICKNESS_RATIOS


FIGURES = HERE / "figures"
STANDARD = {"plate": 4.4360891e-4, "roof": 0.3005924566, "M3": 0.09352155}
METHODS = ("CB-Full", "MidedgeTan", "MITC4")
COLORS = {"CB-Full": "#087e8b", "MidedgeTan": "#bc5f04", "MITC4": "#6748a5"}


def comparable_response(benchmark: str, method: str, n: int, raw: float) -> float:
    if benchmark == "roof" and method == "MITC4":
        # Existing MITC4 driver integrates chordal quad area; CB and MT use
        # the same parametric cylinder area. For this uniform cylinder each
        # chordal area differs by exactly sinc(40 degrees / n). The analysis
        # is linear, so this exact uniform load correction needs no new solve.
        angle = math.radians(40.0) / n
        return raw * angle / math.sin(angle)
    return raw


def read(path: Path) -> list[dict]:
    return list(csv.DictReader(path.open()))


def measured_quality() -> tuple[list[dict], list[dict]]:
    scan = read(DATA / "thickness_scan.csv")
    references = read(DATA / "reference_convergence.csv")
    quality, errors = [], []
    for benchmark in ("roof", "M3"):
        for ratio in THICKNESS_RATIOS:
            family = [r for r in references if r["benchmark"] == benchmark
                      and float(r["t_over_t0"]) == ratio and r["status"] == "SUCCESS"]
            get = lambda m, n: next(comparable_response(benchmark, m, n, float(r["response"])) for r in family
                                    if r["method"] == m and int(r["n"]) == n)
            m64, m96 = get("MITC4", 64), get("MITC4", 96)
            t64, t96 = get("MidedgeTan", 64), get("MidedgeTan", 96)
            self_change = abs(m96-m64)/abs(m96)
            mt_change = abs(t96-t64)/abs(t96)
            cross = abs(t96-m96)/abs(m96)
            extended = benchmark == "M3" and ratio == .0625 and any(
                r["method"] == "MITC4" and int(r["n"]) == 128 for r in family)
            mref = get("MITC4", 128) if extended else m96
            tref = get("MidedgeTan", 128) if extended else t96
            final_step = abs(mref-m96)/abs(mref) if extended else self_change
            final_cross = abs(tref-mref)/abs(mref)
            # This operational check is a study reporting rule, not a proof
            # of convergence or a claim of exact continuum reference.
            quality.append({"benchmark": benchmark, "t_over_t0": ratio,
                "MITC4_n64": m64, "MITC4_n96_surrogate": m96,
                "MidedgeTan_n64": t64, "MidedgeTan_n96": t96,
                "MITC4_64_to_96_rel": self_change,
                "MidedgeTan_64_to_96_rel": mt_change,
                "cross_method_n96_rel": cross,
                "MITC4_reference_n": 128 if extended else 96,
                "MITC4_reference_surrogate": mref,
                "MidedgeTan_reference_n": 128 if extended else 96,
                "MidedgeTan_final": tref,
                "last_refinement_rel": final_step,
                "cross_method_final_rel": final_cross,
                "reference_flag": "screened" if final_step < .01 and final_cross < .02 else "provisional"})
            for r in scan:
                if r["benchmark"] != benchmark or float(r["t_over_t0"]) != ratio or r["status"] != "SUCCESS":
                    continue
                response = comparable_response(benchmark, r["method"], int(r["n"]), float(r["response"]))
                errors.append({"benchmark": benchmark, "t_over_t0": ratio,
                    "n": int(r["n"]), "method": r["method"],
                    "response": response, "surrogate_reference": mref,
                    "relative_discrepancy": abs(response-mref)/abs(mref),
                    "response_ratio": response/mref,
                    "reference_flag": quality[-1]["reference_flag"],
                    "raw_states": int(r["raw_states"]),
                    "active_states": int(r["active_states"]),
                    "time_s": float(r.get("total_wall_s") or r["time_s"])})
    save_rows(DATA/"reference_quality.csv", quality)
    save_rows(DATA/"thickness_discrepancy.csv", errors)
    return quality, errors


def accuracy_states() -> tuple[list[dict], list[dict]]:
    original = read(HERE.parent / "cmame_round4/raw_data/combined_three_methods.csv")
    scan = read(DATA/"thickness_scan.csv")
    rows = []
    for old in original:
        b, m, n = old["benchmark"], old["method"], int(old["n"])
        if b == "M3" and n <= 16:
            # The original CB M3 used a reconstructed point load. The new
            # series uses raw point loads in all three methods.
            candidate = next((x for x in scan if x["benchmark"] == "M3"
                              and float(x["t_over_t0"]) == 1 and int(x["n"]) == n
                              and x["method"] == ("CB-Full" if m == "Checkerboard" else m)), None)
            if candidate is not None:
                response = float(candidate["response"])
            else:
                response = float(old["response"])
        else:
            response = float(old["response"])
        response = comparable_response(b, "CB-Full" if m.startswith("Checkerboard") else m, n, response)
        rows.append({"benchmark": b, "n": n,
            "method": "CB-Full" if m.startswith("Checkerboard") else m,
            "response": response, "reference": STANDARD[b],
            "relative_error": abs(response-STANDARD[b])/abs(STANDARD[b]),
            "raw_states": int(old["raw_scalar_states"]),
            "active_states": int(old["active_states"]),
            "M3_load": "common raw vertex" if b == "M3" and n <= 16 else "original validated"})
    save_rows(DATA/"accuracy_states.csv", rows)
    thresholds = (.05, .02, .01)
    targets = []
    for benchmark in ("plate", "roof", "M3"):
        for method in METHODS:
            family = sorted([r for r in rows if r["benchmark"] == benchmark and r["method"] == method],
                            key=lambda r: r["raw_states"])
            for target in thresholds:
                crossings = [i for i, r in enumerate(family) if r["relative_error"] <= target]
                if not crossings:
                    targets.append({"benchmark": benchmark, "method": method, "target_error": target,
                                    "first_measured_raw_states": "", "persistent": False,
                                    "status": "not attained in measured range"})
                    continue
                first = crossings[0]
                persistent = all(r["relative_error"] <= target for r in family[first:])
                targets.append({"benchmark": benchmark, "method": method, "target_error": target,
                    "first_measured_raw_states": family[first]["raw_states"],
                    "first_measured_active_states": family[first]["active_states"],
                    "first_measured_n": family[first]["n"], "persistent": persistent,
                    "status": "persistent measured attainment" if persistent else "nonmonotone crossing"})
    save_rows(DATA/"target_states.csv", targets)
    return rows, targets


def draw(errors: list[dict], accuracy: list[dict]) -> None:
    FIGURES.mkdir(exist_ok=True)
    plt.rcParams.update({"font.size": 9, "font.family": "DejaVu Sans",
                         "axes.spines.top": False, "axes.spines.right": False})
    fig, axes = plt.subplots(2, 3, figsize=(12, 6), constrained_layout=True)
    for i, benchmark in enumerate(("roof", "M3")):
        for j, method in enumerate(METHODS):
            family = [r for r in errors if r["benchmark"] == benchmark and r["method"] == method]
            ns = sorted(set(r["n"] for r in family))
            grid = np.array([[100*next(r["relative_discrepancy"] for r in family
                              if r["n"] == n and r["t_over_t0"] == ratio)
                              for n in ns] for ratio in THICKNESS_RATIOS])
            ax = axes[i,j]
            im = ax.imshow(np.log10(np.maximum(grid, .01)), cmap="viridis", aspect="auto", vmin=-2, vmax=2)
            ax.set_xticks(range(len(ns)), ns)
            ax.set_yticks(range(len(THICKNESS_RATIOS)), [f"{r:g}" for r in THICKNESS_RATIOS])
            ax.set_title(f"{benchmark}: {method}")
            ax.set_xlabel("cells per direction")
            if j == 0: ax.set_ylabel("thickness / original")
            for a in range(grid.shape[0]):
                for b in range(grid.shape[1]):
                    ax.text(b,a,f"{grid[a,b]:.1f}",ha="center",va="center",fontsize=6.7,
                            color="white" if grid[a,b]>5 else "black")
    fig.suptitle("Scalar response discrepancy (%) to finest MITC4 surrogate; see reference quality table")
    for ext in ("svg", "pdf"): fig.savefig(FIGURES/f"thickness_mesh.{ext}")
    plt.close(fig)

    fig, axes = plt.subplots(1,3,figsize=(12,3.5),constrained_layout=True)
    for ax, b in zip(axes,("plate","roof","M3")):
        for m in METHODS:
            rows = sorted([r for r in accuracy if r["benchmark"] == b and r["method"] == m],
                          key=lambda r:r["raw_states"])
            ax.loglog([r["raw_states"] for r in rows],
                      [100*r["relative_error"] for r in rows], "o-", color=COLORS[m], label=m)
        ax.axhline(1,color="0.6",ls=":",lw=.8)
        ax.axhline(5,color="0.7",ls=":",lw=.8)
        ax.set_title(b); ax.set_xlabel("raw scalar states")
        ax.set_ylabel("standard-thickness scalar error (%)")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    for ext in ("svg","pdf"): fig.savefig(FIGURES/f"accuracy_states.{ext}")
    plt.close(fig)

    ablation = read(DATA/"completion_ablation.csv")
    fig, axes = plt.subplots(1,2,figsize=(8,3.3),constrained_layout=True)
    for ax, (b,n) in zip(axes,(("roof",12),("M3",8))):
        rows = [r for r in ablation if r["benchmark"] == b and int(r["n"]) == n]
        full = float(next(r["response"] for r in rows if r["variant"] == "CB-Full"))
        ax.bar([r["variant"] for r in rows],
               [100*(float(r["response"])/full-1) for r in rows], color="0.4")
        ax.set_title(f"{b}, n={n}"); ax.set_ylabel("response shift from CB-Full (%)")
    for ext in ("svg","pdf"): fig.savefig(FIGURES/f"completion_ablation.{ext}")
    plt.close(fig)

    cylinder = read(DATA/"cylinder_isometry.csv")
    fig, axes = plt.subplots(1,2,figsize=(8,3.3),constrained_layout=True)
    for t in (0.01,0.0025):
        rows = [r for r in cylinder if float(r["thickness"]) == t]
        h = [1/int(r["n"]) for r in rows]
        axes[0].loglog(h,[float(r["metric_relative_frobenius"]) for r in rows],"o-",label=f"t={t:g}")
        axes[1].loglog(h,[float(r["membrane_to_bending"]) for r in rows],"o-",label=f"t={t:g}")
    axes[0].set_ylabel("relative first-form error")
    axes[1].set_ylabel("membrane / bending energy")
    for ax in axes: ax.set_xlabel("h / width"); ax.grid(alpha=.2); ax.legend(fontsize=7)
    for ext in ("svg","pdf"): fig.savefig(FIGURES/f"prescribed_isometry.{ext}")
    plt.close(fig)

    prior = read(HERE.parent/"cmame_additional/05_curved_spectrum/raw_data/spectra.csv")
    fig, axes = plt.subplots(1,2,figsize=(8.4,3.3),constrained_layout=True)
    for ax, geometry in zip(axes,("sphere","saddle")):
        for operator, color, marker in (("J_bW","0.4","o"),
                                        ("J_bW+DXWrel",COLORS["CB-Full"],"s")):
            rows = sorted([r for r in prior if r["curvature"] == geometry
                           and int(r["n"]) == 32 and r["operator"] == operator],
                          key=lambda r:int(r["mode_index"]))
            ax.semilogy([int(r["mode_index"])+1 for r in rows],
                        [float(r["relative_eigenvalue"]) for r in rows],
                        marker+"-", ms=3, color=color, label=operator)
        ax.set_title(f"{geometry}, n=32")
        ax.set_xlabel("ordered sampled mode")
        ax.set_ylabel("eigenvalue / largest")
        ax.grid(alpha=.2)
    axes[0].legend(fontsize=7)
    fig.suptitle("Curvature observation: first 20 modes only, not global nullity")
    for ext in ("svg","pdf"): fig.savefig(FIGURES/f"completion_low_spectrum.{ext}")
    plt.close(fig)


def main() -> None:
    quality, errors = measured_quality()
    accuracy, targets = accuracy_states()
    raw_scan = read(DATA/"thickness_scan.csv")
    runtime = []
    for benchmark in ("roof", "M3"):
        for method in METHODS:
            for n in sorted({int(r["n"]) for r in raw_scan if r["benchmark"] == benchmark}):
                values = [float(r.get("total_wall_s") or r["time_s"]) for r in raw_scan
                          if r["benchmark"] == benchmark and r["method"] == method
                          and int(r["n"]) == n and r["status"] == "SUCCESS"]
                if values:
                    runtime.append({"benchmark": benchmark, "method": method, "n": n,
                                    "thickness_cases": len(values),
                                    "median_setup_to_output_wall_s": float(np.median(values)),
                                    "scope": "driver entry to scalar output; distinct implementations"})
    save_rows(DATA/"runtime_summary.csv", runtime)
    prior = []
    flat_path = HERE.parent/"cmame_round3/01_beta_revalidated/raw_data/summary.json"
    flat = json.loads(flat_path.read_text())["flat_membrane_quotient_diagnostic"]
    prior.extend([
        {"geometry": "flat_plate", "n": 16, "operator": "J_B", "observed_nullity": flat["nullity_c_beta_0"],
         "scope": "complete free in-plane membrane operator", "source": str(flat_path.relative_to(HERE.parents[1]))},
        {"geometry": "flat_plate", "n": 16, "operator": "J_B+DXM", "observed_nullity": flat["nullity_c_beta_1"],
         "scope": "complete free in-plane membrane operator", "source": str(flat_path.relative_to(HERE.parents[1]))}])
    spectrum_path = HERE.parent/"cmame_additional/05_curved_spectrum/raw_data/Table_K1.csv"
    for row in read(spectrum_path):
        if int(row["n"]) == 32 and row["operator"] in ("J_bW", "J_bW+DXWrel"):
            prior.append({"geometry": row["curvature"], "n": 32, "operator": row["operator"],
                          "observed_nullity": row["observed_near_nullity_first_20"],
                          "smallest_observed_positive": row["smallest_observed_positive"],
                          "scope": "near-nullity among first 20, not global nullity",
                          "source": str(spectrum_path.relative_to(HERE.parents[1]))})
    save_rows(DATA/"completion_prior_spectrum_summary.csv", prior)
    draw(errors, accuracy)
    (DATA/"derived_summary.json").write_text(json.dumps({
        "reference_quality": quality, "target_states": targets,
        "notes": ["Only scalar observables are compared, not full-field errors.",
                  "M3 CB original-thickness values are recomputed with raw point loads.",
                  "The 96-grid MITC4 response is a surrogate, qualified by 64-to-96 and cross-method differences."]}, indent=2)+"\n")


if __name__ == "__main__": main()
