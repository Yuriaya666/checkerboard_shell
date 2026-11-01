# Engineering inverse design, mechanics and performance follow-up

Date: 2026-10-29. All dimensions are SI unless explicitly labelled. This report supersedes the **automotive acceptance and performance conclusions** in `RESEARCH_REPORT.md`; that file remains as the earlier reproducible record. The two structures are representative research specimens, not vehicle/aircraft production designs.

## Executive findings

1. The original automotive hard-patch model has a mesh-dependent connection: each refined mesh constrains more panel nodes to the same prescribed translation, while also fixing their x displacement. For a 0.1 mm single-patch z adjustment, the largest Checkerboard mount reaction rises from 245 N (12×8) to 675 N (20×16); ShellMITC4 rises from 502 to 1105 N. The identical 0.1 mm adjustment of *all* patches produces only about 15 N on the fine mesh. This isolates relative patch deformation and connection overconstraint as the dominant source of the reaction blow-up; it does not prove the two shell formulations have identical local bending stiffness.
2. An explicitly hypothetical, area-integrated finite-stiffness bracket model has a fixed physical footprint and total stiffness. A nonlinear Checkerboard optimization directly enforces adjustment, reaction and shape constraints. Its controls pass the 24×20 Checkerboard and independent geometrically nonlinear ASDShellQ4 checks: maximum edge error 0.256/0.256 mm and largest mount reaction 69.1/71.7 N, against 0.6 mm/300 N research limits. The bracket model needs physical calibration and full local stress verification before use for a real car.
3. Aircraft layout `[0,2,4]` remains feasible for both targets in a 20×40 corotational shell/frame check (RMS 0.181/0.206 mm, maxima 0.412/0.644 mm). Isolated position-line versus rotational-frame cantilevers differ by 12.2% in tip compliance at 12 elements; the difference decays with refinement. This prevents a strict element-only comparison of the coupled aircraft models. A three-actuator **minimum** has not been proved.
4. The prior 20×40 Checkerboard 31-layout linear screen spent approximately 29 of 39 seconds in 32 repeated sparse factorizations. Reusing one spar-clamped factorization and a candidate-point compliance matrix reduces the measured 20×40 median screen to 4.67 s. The ShellMITC4 implementation's 37.26 s median rebuilds OpenSees per analysis. This is a workflow result, **not** an intrinsic element speed ratio or a matched-quality nonlinear design speedup.

## Physical and numerical model

The original aircraft and automotive geometries, manufacturing-error field, loads, material inputs, targets and study tolerances are defined in `config.json` and the two case builders. For their provenance and limitations see `RESEARCH_REPORT.md`. The new automotive connector is a **research assumption**: two orthogonal effective flexure leaves per adjustable bracket, each 30 mm long, 30 mm wide, 1 mm thick, with study modulus 210 GPa. The resulting Euler–Bernoulli transverse stiffness per adjustable direction is `3EI/L³ = 58,333 N/m`. It is distributed by bilinear shape-function integration over a fixed 0.10×0.12 m flange footprint. The resulting nodal weights sum to one and preserve the prescribed centroid at every mesh resolution. This is a homogenized connection to a single adjustable lug; it does not resolve actual bolts, adhesive, leaf geometry or contact. The panel's two datum translations remain fixed. The bracket x direction is free in this model, unlike the original patch. Those are **different physical connection hypotheses**, so timing and reaction comparisons across them are diagnostic, not same-model element comparisons.

The flexure is intended to represent a heat-treated high-strength spring-steel part, separate from the 0.7 mm sheet panel. The [Hamilton Precision Metals 17-7 PH strip datasheet](https://www.hpmetals.com/-/media/ametekhpmetals/files/technical-data/stainless-steel/ss-17-7-ph.pdf?la=en&revision=00b0417f-b3e8-4744-bde0-ab1753882605) lists 0.2% yield strengths of roughly 200–260 ksi for cold-rolled/heat-treated conditions and a tensile modulus of about 200 GPa. These are **typical material values, not certified properties of this hypothetical bracket**; the model uses 210 GPa as an engineering approximation. At the Checkerboard optimum the largest effective leaf-root stress estimate `6FL/(bt²)` is about 413 MPa, below those typical high-strength values. The estimate assumes a single resultant per leaf and ignores attachment stress concentrations, fatigue and buckling. A generic mild-steel 1 mm cantilever would not be justified by this estimate.

ShellMITC4 uses small-displacement kinematics. The independent nonlinear check uses OpenSees [ASDShellQ4 with `-corotational`](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/elements/ASDShellQ4.html), which has MITC4 plate bending but different membrane/drilling interpolation. It is **not** a geometrically nonlinear version of the exact ShellMITC4 discretization. All three solvers share panel geometry, thickness, effective bracket energies, design inputs and edge measurements. The aircraft nonlinear reference also uses OpenSees [corotational frame transformations](https://opensees.github.io/OpenSeesDocumentation/user/manual/model/geomTransf/Corotational.html); its frame rotations are coupled to shell rotations, whereas the Checkerboard reinforcement derives curvature from positions. That remaining connector difference is measured below.

## Automotive diagnosis and verified design

The original fixed-input experiment applies zero, one y, one z, opposite z, all-patch z and saved old-optimum inputs on 12×8, 16×12 and 20×16 meshes. The region nodal counts are 3, 6 and 15 per patch. Reactions are **component sums across the patch followed by a vector norm**, excluding the two datum reactions; summing nodal norms would be a different and usually larger quantity. OpenSees' rotational nodal reactions are not included in the reported translational mount resultant. The force balances, spatial displacements and energy/work diagnostics are saved in `results/automotive_fixed_input_diagnostics.json`.

| Fixed 0.1 mm input | 12×8 CB / MITC4 (N) | 16×12 (N) | 20×16 (N) |
|---|---:|---:|---:|
| One patch, z | 245 / 502 | 367 / 690 | 675 / 1105 |
| Opposite patch z | 449 / 1001 | 683 / 1323 | 1271 / 2202 |
| All patches, z | 14.7 / 26.1 | 14.2 / 15.4 | 15.2 / 14.6 |

With the finite bracket, the same one-patch z input produces 0.521/1.009 N at 12×8, 0.486/0.809 N at 16×12, and 0.462/0.703 N at 20×16. Both methods trend toward mesh-stable reaction, although their finite-mesh values remain different. The old hard patch has no such trend at these resolutions. The 0.5×/1×/2× bracket-stiffness sensitivity at 12×8 gives Checkerboard maximum edge errors 0.217/0.150/0.358 mm and maximum reactions 35.3/69.6/136.9 N under the *same* optimum controls. This demonstrates that the result depends on bracket calibration; it was not tuned by silently relaxing the 300 N limit. Raw results are in `results/automotive_spring_audit.json` and `results/automotive_spring_reference.json`.

The direct nonlinear optimization uses eight y/z bracket displacements in millimetres, SLSQP bounds ±1.5 mm, four explicit resultant constraints ≤300 N, and all inspected gap/flush components within ±0.6 mm. The objective is the original weighted squared edge error plus a small stroke penalty. The consistent shell tangent, bracket stiffness, implicit state derivative and bracket-reaction derivative are used; a bordered sparse solve avoids cancellation in the local-plus-rank-one tangent. Failed direct Newton trials are retried with adaptive proportional input continuation. The final result required 42 optimizer iterations, 72 equilibria/linearizations and 328 s in the recorded run. All final constraints are verified from a converged nonlinear state, with free residual `5.95e-7 N`; there is **no post-optimization stroke scaling**. The force limit is inactive at this optimum. The first recorded run predates the added per-evaluation history field; rerunning the script writes a new history alongside the summary metrics. Repeated nonlinear optimization may select another local minimum, so the saved design is not assumed bit-for-bit reproducible across solver changes.

| Bracket | y adjustment (mm) | z adjustment (mm) | Checkerboard resultant (N) |
|---|---:|---:|---:|
| 1 | +0.780 | −0.405 | 48.5 |
| 2 | −1.252 | −0.581 | 69.6 |
| 3 | −1.007 | +0.093 | 53.5 |
| 4 | +0.948 | +0.132 | 59.6 |

Checkerboard 12×8 gap RMS is 0.080 mm, flush RMS 0.095 mm, maximum edge error 0.150 mm and maximum displacement 1.221 mm. The minimum gap is 2.850 mm, so the specified inspection edge does not cross its nominal adjacent reference edge. The largest sampled panel membrane principal strain is `7.76e-6`, but no complete bending-stress or flange-joint fatigue assessment has been made. The spring energy is 0.118 J and reported shell energy 0.00174 J; the former dominates this particular connector model.
At 24×20, the same panel with zero bracket adjustment has maximum edge error 0.633 mm, gap RMS 0.125 mm and flush RMS 0.463 mm. The optimized controls therefore improve the common-grid maximum and RMS contour errors, rather than merely passing a tolerance on the coarse design mesh. The previous fixed-datum rigid adjustment baseline in `automotive_12x8.json` has maximum error 0.679 mm; that baseline assumes only a rotation about the line through two datums and is deliberately limited, so it is evidence of deformation benefit for this fixture assumption rather than a universal rigid-registration lower bound.

**Same-design independent check:**

| Mesh | Checkerboard error / force | ShellMITC4 linear error / force | ASDShellQ4 corotational error / force |
|---|---:|---:|---:|
| 12×8 | 0.150 mm / 69.6 N | 0.376 mm / 74.9 N | 0.152 mm / 72.7 N |
| 16×12 | 0.215 / 69.2 | 0.321 / 74.3 | 0.193 / 72.0 |
| 20×16 | 0.247 / 68.9 | 0.278 / 74.0 | 0.239 / 71.7 |
| 24×20 | 0.256 / 69.1 | 0.258 / 73.7 | 0.256 / 71.7 |

Numbers are maximum absolute gap/flush error and largest area-resultant mount reaction. All shown solves report convergence, and force balances are near numerical zero. The common-input 24×20 error difference between Checkerboard and corotational ASDShellQ4 is below 0.001 mm, while resultant reaction differs by about 3.8%. This does **not** establish convergence to the exact real vehicle response. A moment audit shows Checkerboard's current-configuration force moment residual below `4e-9 N·m` for the saved designs, but the OpenSees summed *translational* reactions leave 0.065–0.126 N·m at 24×20. The OpenSees shell also has internal rotational moments; the two datum rotational reactions were near zero in a single-patch probe. A complete shell moment ledger has not been assembled, so **global angular-momentum equivalence remains unresolved** and the reference is accepted here for common-input shape/resultant-force checks only.

To separate element formulation from kinematics, the same 12×8/16×12/20×16/24×20 controls were also evaluated with **ASDShellQ4 without corotation**. Its maximum edge errors are 0.185/0.239/0.274/0.284 mm. The corresponding corotational values are 0.152/0.193/0.239/0.256 mm; ShellMITC4 linear values are 0.376/0.321/0.278/0.258 mm. Thus the difference is jointly affected by geometry updates and by the two OpenSees shell formulations, especially on coarse meshes. The added table is in `results/automotive_spring_reference.json`.

The separate 12×8 ShellMITC4 *linear* optimum reaches a 0.242 mm edge error and 86.9 N on its own mesh, but its same controls yield 1.403 mm maximum edge error in the 20×16 corotational reference. It therefore fails the common nonlinear shape criterion despite satisfying the linear model. Its 0.315 s linear design time cannot be compared as an equal-quality competitor with the 328 s nonlinear Checkerboard run. `results/automotive_mitc4_spring_design.json` records this counterexample. A separate **corotational** local refinement starts from the Checkerboard-feasible controls. At 12×8 its independently verified maximum error is 0.143 mm and resultant reaction 82.5 N. At 24×20 the same controls give Checkerboard/corotational maximum edge errors 0.294/0.321 mm and largest reactions 74.2/80.8 N; both pass the research constraints. This establishes that the independent nonlinear shell can improve a feasible design, but it uses a **different starting point** and its runtime is not an equal-start acceleration experiment. The full crosscheck is in `results/automotive_design_crosscheck.json`.

## Aircraft layout and connector checks

The existing full 31-layout **linear** screen remains reproducible. The three-position `[0,2,4]` layout at span fractions 0.05/0.50/0.95 is a demonstrated feasible design, not a proven cardinality minimum. Its two study targets, 4 mm stroke bounds and 80 N net actuator-reaction limits are unchanged. On the 12×24 Checkerboard nonlinear model, uniform and span-gradient target RMS/maxima are 0.146/0.412 mm and 0.168/0.599 mm. Their strokes are `(1.591,1.630,1.588)` mm and `(0.692,1.616,2.501)` mm. All net reactions are below 3.5 N. The 20×40 independent corotational shell/frame results for the **same** strokes are 0.181/0.412 mm and 0.206/0.644 mm, with largest net reactions 1.26 and 1.66 N. Both targets pass 0.25 mm RMS/0.8 mm maximum study tolerances.

The single-point and two-point nonlinear audit now retries a failed direct equilibrium with proportional pressure/stroke continuation and adaptive step size. All **5 single-point and 10 two-point** layouts obtained at least one converged nonlinear design; none of the tested local designs passed both target conditions. The previous `[1,2]` twist-case solver failure was resolved; the obtained equilibrium design remains far outside tolerance (RMS 0.539 mm, maximum 2.391 mm). `[2,4]` uniform case now resolves at RMS 0.283 mm and maximum 1.186 mm, also outside tolerance. The best observed one-point twist RMS is 0.308 mm (`[4]`), and best observed two-point twist RMS is 0.164 mm (`[2,4]`); the latter fails the uniform case. Some other uniform cases were not run after their twist case failed. The audit records status `not_accepted_local` rather than asserting global infeasibility. These are *local optimization outcomes*; no certified nonlinear global lower bound was computed. The full audit JSON and CSV record each layout status and must be interpreted accordingly.

For the same `[0,2,4]` strokes, Checkerboard nonlinear mesh checks at 12×24, 16×32 and 20×40 give uniform RMS 0.146/0.138/0.148 mm and twist RMS 0.168/0.159/0.169 mm. The maximum twist error is 0.599/0.624/0.639 mm; all remain within the study tolerance. Net reactions and solve histories are in `results/aircraft_checkerboard_meshcheck.json`. This is stability of the *chosen design*, not convergence of the discrete optimal layout.

For an isolated 0.18 m chordwise cantilever with matched E, A and I and a 1 N tip load, the position-line tip compliance is 0.878 of the rotational-frame value with 12 segments; the ratio is 0.908/0.926/0.938/0.963 for 16/20/24/40 segments. The two position nodes used to approximate a clamped slope differ from one frame node with fixed rotations. This confirms a measurable finite-mesh connector difference that diminishes with refinement. The full coupled aircraft comparison cannot attribute its remaining differences solely to the shell element.

## Performance and attribution

The original 20×40 profile for **base assembly plus 31 layout operators** took 39.17 s: sparse factorization 29.33 s, dense/sparse conversions several seconds, and the one shell tangent AD assembly 2.44 s. The original `common.py` built a dense global rank-one outer product and re-factorized each layout. The revised linear layout operator solves the *unchanged* spar-clamped system once for all candidate DOFs and uses the selected compliance submatrix for exact Lagrange multipliers/reactions. At 12×24, fast and direct elimination differ by at most `7.4e-10` in the response map and `4.8e-11` relative in the actuator-reaction map. No shell energy, Blind-Mode Completion or constitutive definition was changed.

| Aircraft mesh | Prior CB (s) | Revised CB median (s) | ShellMITC4 median (s) | Scope |
|---|---:|---:|---:|---|
| 12×24 | 3.62 | 1.24 | 9.46 | all 31 linear layouts |
| 16×32 | 11.87 | 2.15 | 23.15 | all 31 linear layouts |
| 20×40 | 39.92 | 4.67 | 37.26 | all 31 linear layouts |

The revised Checkerboard totals include tangent assembly, one spar-clamped factorization/candidate compliance solve and all layout optimizations. The MITC4 path performs 111 separate OpenSees analyses and rebuilds the model each time. Both paths use one CPU-thread setting in `performance_round2.py`; raw repeat times and environment are in `results/performance_round2.json`. The measured host has an Intel Xeon Gold 6148 CPU at 2.40 GHz, Python 3.10.22, OpenBLAS one thread and PyTorch one thread. The prior 12×24 number was a three-run median; prior finer numbers were single runs. The revised comparison identifies a successful **Checkerboard implementation optimization**, while the remaining OpenSees rebuild overhead prevents an element-level speed conclusion. The automotive nonlinear design was substantially more expensive than its invalid ShellMITC4 linear surrogate, and no matched-quality full nonlinear runtime ratio is reported.

The separately rerun 20×40 `cProfile` audit measured 44.94 s for the old base-plus-31-layout operator and 3.91 s for the shared-compliance operator. The old path spent 33.75 s in 32 SuperLU numeric factorizations; the new path spent 1.01 s in a single factorization, 1.80 s in local Hessian AD assembly, 0.178 s preparing the compliance basis and 0.102 s across the 31 layout operators. `results/profile_round2.json` contains the raw function table. These are operator-only profile timings, not the optimization-inclusive table above.

The original nine-force flat-plate Round 5 adjoint check was rerun: five sampled components have maximum relative gradient error `6.07e-12` against central differences. The same-model nine-force design loop was also rerun: Adjoint required 68 linear solves/1.076 s total, central FD 391 solves/1.908 s total, with nearly equal final target errors (`6.67e-5`/`6.73e-5` relative). This supports same-model gradient correctness and efficiency; it does not explain the whole cross-method runtime gap.

For the new automotive model at the final design, one nonlinear baseline equilibrium took 4.40 s, the exact implicit response/reaction linearization 0.675 s, and a centered eight-variable finite-difference objective gradient required 16 additional equilibria and 34.59 s. Their maximum absolute objective-gradient difference was `4.31e-5`, about 0.60% of the largest gradient scale at the chosen 0.01 mm difference step. This measures a **within-Checkerboard gradient-method benefit**, separate from any shell-element comparison (`results/gradient_cost_round2.json`).

## Verification and reproducibility

The full current test suite reports **29 passed**. New tests check area-integrated bracket weight sums/centroids at two meshes, exact response/reaction equivalence of shared-compliance and direct aircraft elimination, and automotive shape/reaction derivatives against central differences. The bracket reaction derivative differs from centered differences by less than 0.1 N/m absolute tolerance under the tested perturbation; response sensitivities agree within `3e-4` dimensionless absolute tolerance. The production tangent itself was checked against an independent automatic-differentiation Hessian-vector product at the optimized automotive state, with relative norm differences around `4e-16`.

Run from the repository root with an environment containing PyTorch, NumPy, SciPy, Matplotlib, pytest and OpenSeesPy:

```bash
export PYTHONPATH=src:.
export OMP_NUM_THREADS=1 MKL_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1
python -m pytest -q
python -m experiments.engineering_inverse.diagnose_reactions
python -m experiments.engineering_inverse.mount_springs
python -m experiments.engineering_inverse.automotive_spring_design
python -m experiments.engineering_inverse.automotive_spring_reference
python -m experiments.engineering_inverse.automotive_mitc4_spring_design
python -m experiments.engineering_inverse.automotive_corotational_design
python -m experiments.engineering_inverse.automotive_design_crosscheck
python -m experiments.engineering_inverse.gradient_cost_round2
python -m experiments.engineering_inverse.aircraft_nonlinear_reference
python -m experiments.engineering_inverse.aircraft_meshcheck
python -m experiments.engineering_inverse.nonlinear_audit
python -m experiments.engineering_inverse.beam_kinematics_audit
python -m experiments.engineering_inverse.performance_round2
python -m experiments.engineering_inverse.profile_round2
python -m experiments.engineering_inverse.make_round2_figures
```

The modules write JSON/CSV under `results/`, PDF/SVG under `figures/`, and can be rerun independently. `python -m experiments.engineering_inverse.run_round2` is the one-command full entry; `--quick` rebuilds fixed-input diagnostics, saved-design references and figures, provided the two saved designs exist. **The quick entry was executed end-to-end**; its stage times are in `results/round2_execution_log.json`. The full entry's long stages were executed separately and their outputs saved, but the full entry has **not** been run as one continuous job. `performance_round2.py` overwrites the saved aircraft linear-screen JSON with the same configuration and records all repeats. `nonlinear_audit.py` and `aircraft_meshcheck.py` resume from their partial JSON output. The older `run_all.py` continues to reproduce the **prior** hard-patch experiment; the commands above reproduce the present connector study. The bracket dimensions and objective limits are in `config.json`.

## What the evidence supports

- A representative automotive elastic-assembly inverse design **passes the stated research tolerances in two independently implemented geometrically nonlinear shell models** under a clearly specified effective bracket hypothesis. It is not a calibrated front-fender prediction, and the OpenSees moment ledger remains unresolved.
- A three-position aircraft layout meets both stated targets under the current Checkerboard model and a corotational shell/frame reference. A cardinality minimum and a fully identical reinforcement kinematics comparison remain unproved.
- The Checkerboard 31-layout **linear workflow** can be made substantially faster by exact factorization reuse. The data do not establish an intrinsic position-only element speed advantage, nor an equal-quality complete nonlinear task speedup over an optimized MITC4 implementation.
- The old hard automotive patch and the old ShellMITC4 linear automotive optimum are counterexamples to unqualified cross-method acceptance claims. They remain in the raw data and figures.

Outstanding work for a paper claim of generality is to calibrate or test the bracket stiffness and local stress, obtain a strictly matched nonlinear MITC4/Checkerboard reinforcement connection, resolve the aircraft minimum-layout question with a certified or much stronger global nonlinear feasibility search, and perform equal-start nonlinear design timing with comparably optimized implementations.
