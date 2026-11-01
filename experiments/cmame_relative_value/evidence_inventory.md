# Existing evidence and missing discriminating comparisons

Source audit: the original manuscript [arXiv:2608.30195v1](https://arxiv.org/pdf/2608.30195), `docs/CMAME_PROJECT_SYNTHESIS.md`, Round 4 unified comparison, and the supplementary β/spectrum experiments. The original paper's thickness–mesh M3 figure establishes a *within-CB* trend; it does not provide a per-thickness, cross-method reference.

| Status before this study | Evidence | Limit |
|---|---|---|
| Established analytically | B-face/W geometry, quotient gauge, flat observation blind directions, scale-correct (X_M) and (X_W^{rel}), and (O(h^2)) metric approximation for aligned generalized cylinders | The theorem covers a restricted isometry family, not all curved shells. |
| Established analytically | Modal factor (t^2/(\lambda+t^2)) at fixed mesh | It predicts residual locking for load-bearing modes with \(\lambda>0\); no general locking-free result follows. |
| Verified numerically | Flat membrane operator nullity drops from 6 to 5 upon adding (X_M); β sweep shows target-direction stiffness changes linearly | Curved roof/M3 boundary-reduced *full* tangents have no extra zero eigenvalue when (X_M) is removed. |
| Verified numerically | Curved sphere/saddle observation spectra: (J_{bW}+DX_W^{rel}) stiffens missing soft channels; M3/roof β sweep has small ordinary-response changes | Only first 20 eigenpairs were scanned at some resolutions; no external modal baseline. |
| Verified numerically | Plate, modern-reference roof, and M3 standard-thickness CB/MidedgeTan/MITC4 scalar-observable sequences in Round 4 | M3 CB originally used reconstructed point loads, while MT/MITC4 used raw points; the methods also implement symmetry constraints differently. |
| Numerical indication | M3 CB raw-observable error at standard thickness is small at n=12–16 with fewer states; roof and plate show no uniform CB accuracy advantage | One favorable point is not a general accuracy-per-state or locking result. |
| Missing before this study | Two curved benchmarks over several thicknesses and meshes with all three implementations and per-thickness reference-quality check | Added in `run_scan.py`, with reference convergence flags. |
| Missing before this study | Four-way (X_M/X_W^{rel}) ablation under identical loads and boundaries | Added in `run_completion.py`; full-spectrum nullity must still be interpreted separately from observation-rank defects. |
| Still missing | Verified BSQ or other same-3V quadrilateral implementation under identical conditions | The published BSQ description does not provide a reproducible adapter here; no quantitative BSQ ranking is made. |
| Still missing | Shared full-field displacement errors and a certified continuum limit at every ultra-thin point | This study reports scalar response discrepancies and flags unresolved reference convergence. |

The new study retains the existing shell energy and all original benchmark definitions. The only production-code extension is an ablation switch that omits (X_M) while retaining (X_W^{rel}); it does not change CB-Full.
