# Experiment 3: dense and low-rank costs

`run_primary.py --only scaling` benchmarks plate and M3, three formulation
ablations, three post-warmup repeats at n=16/32 and explicit reference
tangents at n=16.  Requested n=64/128 are present as `RESOURCE_GUARD` rows,
with exact byte-count projections of the four dense boundary reconstruction
matrices and one dense tangent; these are **not measured timings**.  Timing
data and sparsity are in `scaling.csv`, rank-one Hessian checks in
`rank_one_audit.json`, exact augmented sparse versus naive dense reference
linear-solve checks in `low_rank_equivalence.csv`, and derived Table C1 in
`Table_C1.csv`.

The sparse/augmented variant changes neither the energy nor its reference
linear tangent: it stores a sparse local matrix and one global gradient vector.
Its current prototype *still obtains the local matrix via dense AD*, so
prototype end-to-end assembly and memory are not sparse-solver production
performance.  CPU timings cannot be compared to absent external baselines.
Figures C1/C2 distinguish measured tangent actions from predicted storage;
C3 uses a declared relative numerical nnz threshold.

