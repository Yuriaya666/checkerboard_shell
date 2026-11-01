# Checkerboard thin-shell relative numerical value study

**Scope.** This study tests the accepted CB-Full model on the existing plate, Scordelis–Lo roof, and MacNeal–Harder M3 benchmarks. It adds two curved-shell thickness scans, common scalar-observable comparisons, a four-way completion ablation, and one prescribed exact-cylinder isometry check. It does not alter the CB-Full energy or claim an inverse-design CPU acceleration.

## 1. Prior evidence and experiment definition

The [original manuscript](https://arxiv.org/pdf/2608.30195) already establishes the flat observation-rank defects, quotient gauge, the (X_M/X_W^{rel}) definitions, aligned-cylinder (O(h^2)) near-isometry approximation, and fixed-grid residual-locking relation (t^2/(\lambda+t^2)). Its M3 thickness figure is a within-CB trend. `evidence_inventory.md` lists proven, indicative, and missing evidence before this comparison.

The present curved scans use four meshes each (roof (n=8,12,16,24); M3 (n=6,8,12,16)) and (t/t_0=1,1/2,1/4,1/8,1/16), where (t_0=0.25) for roof ((R=25)) and (t_0=0.04) for M3 ((R=10)). Dead loads scale as ((t/t_0)^3); this makes a pure-bending limit thickness independent without assuming the roof or M3 is pure bending. M3 uses identical raw equator point loads and the raw equator observable for all three methods. Roof pressure uses the same parametric cylinder area after an exact chord-to-arc load correction for the linear MITC4 output. No material, completion coefficient, or structural boundary was tuned to improve CB.

Reference surrogates were computed separately at 48, 64, and 96 cells per direction with official LibShell MidedgeTan and OpenSees ShellMITC4. The thinnest M3 point was extended to 128. Our *operational reference screen* requires the final MITC4 refinement change below 1% and the two external methods' finest responses to differ below 2%. All ten thickness points pass, but this is not a certified continuum limit. The standard-thickness target-state table instead uses the published plate/M3 reference values and the modern roof value 0.3005924566. Only scalar displacement errors are reported; a common projected full-field norm was not available. At original thickness, the standard M3 reference 0.09352155 and 96-grid surrogate 0.09354998 differ by 0.0304%; the modern roof reference 0.30059246 and 96-grid MITC4 surrogate 0.30153166 differ by 0.312%, so the standard roof accuracy table and thickness-scan surrogate answer different, explicitly labeled questions.

## 2. Thickness and mesh: the relative result

The following percentages are scalar-response discrepancies to the finest MITC4 surrogate at **the same thickness**, rather than errors to a reference carried over from (t_0). A 5% discrepancy can be interpreted more strongly than a 1% discrepancy because the external finest-grid spread is up to 1.7% at (n=96) and 0.97% after the M3 extension.

| Benchmark and fixed mesh | (t/t_0) | CB-Full | MidedgeTan | MITC4 |
|---|---:|---:|---:|---:|
| Roof, (n=24) | 1 | 3.96% | 1.92% | 0.21% |
| Roof, (n=24) | 1/4 | 6.17% | 6.55% | 0.42% |
| Roof, (n=24) | 1/16 | 11.96% | 11.55% | 0.51% |
| M3, (n=16) | 1 | 0.87% | 4.98% | 3.64% |
| M3, (n=16) | 1/2 | 1.66% | 3.51% | 11.62% |
| M3, (n=16) | 1/4 | 4.10% | 0.32% | 33.89% |
| M3, (n=16) | 1/8 | 9.34% | 6.21% | 67.00% |
| M3, (n=16) | 1/16 | 19.73% | 17.28% | 88.92% |

The complete 2 × 5 × 4 × 3 scan is in `results/thickness_scan.csv` and the per-thickness comparisons in `results/thickness_discrepancy.csv`. The M3 (n=16) normalized responses at the thinnest point are 0.07315 (CB), 0.07538 (MidedgeTan), and 0.01010 (MITC4); the corrected 128-grid MITC4 surrogate is 0.09113. The result establishes substantially weaker **observed fixed-grid response suppression than this MITC4 implementation** on M3. It does not show a CB advantage over MidedgeTan in the ultra-thin M3 regime: MidedgeTan is slightly closer at (t/t_0\leq1/4). CB also retains substantial residual suppression at fixed (n), so “locking-free” is false.

Roof does not support a CB locking advantage. Its initial curvature and diaphragm induce membrane–bending coupling; in the reference-tangent quadratic energy at (n=24), CB's (E_m/E_b) is 0.878 at (t_0) and 1.419 at (t_0/16). The roof's (t^3)-normalized response changes substantially even in the refined external solutions. A shrinking normalized roof response alone is therefore not a valid locking diagnosis.

## 3. Accuracy per independent state

The table gives the **first measured** raw-state count that attains each standard-thickness scalar target. It does not extrapolate beyond the tested meshes. All-DOF active counts and every original data point are in `results/accuracy_states.csv` and `results/target_states.csv`.

| Benchmark | Target | CB-Full | MidedgeTan | MITC4 |
|---|---:|---:|---:|---:|
| Plate | 5% | 867 | 451 | 486 |
| Plate | 1% | 3267 | 6403 | 486 |
| Roof | 5% | 1875 | 1667 | 1014 |
| Roof | 1% | not attained through (n=32) | 6403 | 1014 |
| M3 | 5% | 507 | 963, nonmonotone crossing | 1734 |
| M3 | 1% | 867 | not attained through (n=16) | 6534 |

At the M3 1% target CB uses 867 raw scalar states and MITC4 first attains the target with 6534, a **7.54-fold state-count difference** in this measured range. That is a problem-specific accuracy–state advantage, not a CPU acceleration and not proof that CB outperforms a same-3V shell. On the roof MITC4 is clearly better; on the plate it also reaches 1% on the coarsest measured grid. MidedgeTan's M3 5% crossing at (n=12) is not persistent because its (n=16) error grows again; the table does not treat it as a converged threshold. The M3 CB values here are recomputed with raw loads, correcting the load mismatch in the older Round 4 comparison.

## 4. Completion selectivity and near-isometry

The four variants CB-0, CB-M, CB-W, and CB-Full share the same geometry, material, loads, and boundaries. On the roof (n=12), removing both terms shifts the scalar response by −0.15% relative to Full; on M3 (n=8), it shifts it by +0.94%. At these meshes the lowest constrained physical eigenvalue changes by only +0.26% (roof) and +1.03% (M3) from CB-0 to Full. For the Full response, the quadratic (X_M) energy fraction is about (8.1\times10^{-7}) on the roof and (8.9\times10^{-14}) on M3; (X_W^{rel}) contributes about 0.18% and 0.93%, respectively. The existing flat membrane-operator test shows nullity 6→5 when (X_M) is added. The existing spherical/saddle observation spectra show strong lifting of sampled soft curvature channels by (X_W^{rel}), but first-20 near-null counts are **not** a global nullity theorem. In the constrained curved roof/M3 full tangent, all four variants have zero additional numerical nullity; one must not claim that (X_M) removes a curved-shell full-tangent zero mode.

For the prescribed flat-to-cylinder exact continuum isometry ((R=1.5)), the CB first-form relative error falls from (4.091\times10^{-4}) at (n=8) to (2.557\times10^{-5}) at (n=32), consistent with (O(h^2)). At (t=0.0025), the membrane/bending energy ratio falls from 0.362 to 0.00141 across the same refinement. Both completion energies are negligible for this aligned family. At fixed (n=8), reducing thickness from 0.01 to 0.0025 raises (E_m/E_b) from 0.0226 to 0.362, consistent with residual (t^{-2}) amplification. This is **prescribed-deformation** evidence for geometric and energy consistency, not a loaded equilibrium comparison or a proof of uniform no-locking behavior. Stored mode shapes and rendered mode views from the earlier curved-spectrum experiment remain in `experiments/cmame_additional/05_curved_spectrum/raw_data/mode_shapes/` and `experiments/cmame_additional/05_curved_spectrum/figures/Figure_K3.svg`; this study's `figures/completion_low_spectrum.svg` plots their associated first-20 spectrum.

## 5. Cost and model comparability

The latest setup-to-output median wall times over the five thicknesses are 1.50 s (CB), 0.137 s (MidedgeTan), and 0.230 s (MITC4) for roof (n=24); for M3 (n=16) they are 0.815, 0.065, and 0.080 s. CB is slower in these implementations despite fewer states. Its larger PyTorch/automatic-differentiation assembly cost offsets the sparse-factor size reduction. Roof (n=24) has 1875/1800 raw/active CB states versus 3651/3552 MidedgeTan and 3750/3625 MITC4. The exposed sparse tangent has 136031 CB and 175784 MidedgeTan structural entries; OpenSees does not expose an equivalent assembled matrix count in this adapter. Timings include each driver invocation and necessary setup, but use different language stacks and solver libraries; they do not isolate intrinsic element cost. `results/runtime_summary.csv` records every mesh.

Three comparability limits remain material. (1) MITC4 is a Mindlin–Reissner shell while CB and MidedgeTan use Kirchhoff–Love/Koiter kinematics; the plate MITC4 has a known shear contribution at finite thickness. (2) Physical symmetry is imposed through midpoint-position constraints in CB and through rotational or edge-variable constraints externally; this is especially consequential on M3. Raw loading and output are now common, but the discrete boundary kinematics are not identical. (3) A reliable BSQ or other same-3V quadrilateral implementation was not available with matching geometry, load and boundary definitions, so this study **cannot** rank CB against the broader position-only class. The published [BSQ element paper](https://www.sciencedirect.com/science/article/pii/S0045-7825%2807%2900020-5) describes a translational-only element, but no external BSQ data are mixed into the quantitative plots.

## 6. Direct answers and manuscript positioning

- **Q1, near-isometric bending:** CB has a reproducible M3 advantage over MITC4 on the tested coarse and moderately thin grids, and comparable order to MidedgeTan. The roof provides a counterexample to broad superiority.
- **Q2, residual membrane locking:** M3 CB suppresses the thin response much less than tested MITC4 at equal (n), but at (t_0/16,n=16) CB still underresponds by 19.7%, slightly more than MidedgeTan's 17.3%. No general or even M3-specific “locking-free” statement is justified.
- **Q3, accuracy/state:** At M3 standard thickness and 1% scalar error, CB uses 7.54× fewer raw states than the first measured MITC4 crossing. The plate and roof reverse or erase that advantage.
- **Q4, completion:** The flat observation-rank defect is repaired, sampled curvature soft channels are lifted, and tested ordinary responses/eigenvalues shift modestly. This is the strongest distinctive, theoretically interpretable numerical evidence; it is not a universal guarantee against stiffening.
- **Q5, Checkerboard geometry:** The aligned cylinder test connects (O(h^2)) near-isometry to small membrane energy, while also showing fixed-grid thin-limit amplification. No direct same-3V competitor exists here to prove exclusivity.
- **Q6, other 3V shells:** Quantitative superiority remains unresolved without a verified same-condition BSQ or comparable baseline.
- **Q7, CMAME claim:** Position the contribution as a compact 3V B/W shell with quotient-aware, scale-correct blind-mode completion, verified M3-specific accuracy/state efficiency and limited thin-regime improvement over the tested MITC4. State openly that roof/plate accuracy and measured CPU time do not show a general CB advantage.

**Integrated judgment.** Extending Checkerboard geometry to shell mechanics has produced repeatable *selective* value: an aligned near-isometry/energy mechanism, targeted completion with small normal-response disturbance, and unusually favorable M3 accuracy per state at the original thickness. The same experiments rule out a universal thin-shell accuracy, membrane-locking, or runtime advantage. Full-field error, a certified ultra-thin continuum reference, and a verified same-3V comparison remain substantive limits on stronger claims; new aircraft or automotive studies are not needed to state this conclusion.

## Reproduction and validation

`README.md` gives environment-independent commands and source provenance. The modified MidedgeTan adapter reproduced its original plate/roof/M3 outputs exactly (`results/midedgetan_gate.json`). The added CB-W switch preserves CB-Full and passes a deformed-state exact additivity check for residuals and tangent actions; original roof and M3 benchmark responses also pass fixed-value regression gates. The complete current test suite passed: **35 tests** (`full_regression.log`). All plots in `figures/` have SVG and PDF versions and are generated from the result CSV by `analyze.py`. Original Round 4 raw data and benchmarks remain in their original locations.
