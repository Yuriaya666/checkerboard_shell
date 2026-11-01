"""Collect the validated manuscript figure set without regenerating old studies."""

from __future__ import annotations

import json
from pathlib import Path
import shutil


ROOT = Path(__file__).resolve().parent
REPOSITORY = ROOT.parents[1]
OUT = ROOT/"03_final_figures"

SELECTIONS = {
    "Figure_A1_gauge_physical_quotient": (
        REPOSITORY/"experiments/cmame_additional/02_gauge_invariance/figures/Figure_G1",
        "Round 1; benchmark-driver-independent and rerun by Round-2 regression"),
    "Figure_A2_gauge_invariance_errors": (
        REPOSITORY/"experiments/cmame_additional/02_gauge_invariance/figures/Figure_G2",
        "Round 1; benchmark-driver-independent and rerun by Round-2 regression"),
    "Figure_B_accuracy_per_state": (
        ROOT/"02_final_comparison/figures/Figure_Comparison_accuracy_per_state",
        "Round 3; recovered manuscript drivers and official LibShell baseline"),
    "Figure_C_sparse_rank_one_scaling": (
        REPOSITORY/"experiments/cmame_round2/03_scaling/figures/Figure_C3_nnz_scaling",
        "Round 2; exact local sparse plus rank-one production path"),
    "Figure_D_surface_preserving_distortion": (
        REPOSITORY/"experiments/cmame_round2/04_distortion_revalidated/figures/Figure_D1_distortion_comparison",
        "Round 2; recovered drivers and identical surface-preserving grids"),
    "Figure_E_beta_sensitivity": (
        ROOT/"01_beta_revalidated/figures/Figure_Beta_revalidated",
        "Round 3; recovered manuscript drivers"),
}


def run() -> dict:
    OUT.mkdir(parents=True, exist_ok=True)
    manifest = []
    for destination, (source, provenance) in SELECTIONS.items():
        for suffix in ("pdf", "svg", "png"):
            origin = source.with_suffix(f".{suffix}")
            target = OUT/f"{destination}.{suffix}"
            if not origin.is_file():
                raise FileNotFoundError(origin)
            shutil.copy2(origin, target)
        manifest.append({
            "figure": destination,
            "source_stem": str(source.relative_to(REPOSITORY)),
            "provenance": provenance,
            "formats": ["pdf", "svg", "png"],
        })
    (OUT/"figure_manifest.json").write_text(
        json.dumps(manifest, indent=2)+"\n", encoding="utf-8")
    return {"figure_sets": len(manifest)}


if __name__ == "__main__":
    print(json.dumps(run(), indent=2))
