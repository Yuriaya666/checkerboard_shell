# X_M derivative structure audit

The production source evaluates

`E_XM = (t/4) (8 mu) X_M^2 = c X_M^2`, with `c = 2 t mu`.

Consequently its exact finite-configuration residual and tangent are

`R_XM = 2 c X_M g`,

`K_XM = 2 c [g g^T + X_M H_X]`,

where `g = D X_M` and `H_X = D2 X_M`.  The implemented rank-one coefficient
is therefore `2c = 4 t mu`; no coefficient was inferred from the prompt.

On a finitely deformed curved M3 n=6 state, the direct dependency audit found:

- `g` has 147 entries and is globally supported;
- `g g^T` has rank exactly one and all 21,609 entries are nonzero;
- `H_X` has 6,669 numerical nonzeros out of 21,609 (density 0.3086 on this
  very small patch);
- every `H_X` coupling lies within two raw-node steps in each logical
  direction;
- summing the individually weighted plaquette observations reproduces `X_M`
  to `1.36e-20`.

Thus `X_M = sum_p x_p(q_local)` and `D2X_M = sum_p D2x_p` are local.  Only
the outer product is globally dense.  The production tangent stores
`2c X_M H_X` in the coloured sparse local Hessian and stores `g` plus `2c`
for the global correction.  It never forms `g g^T`.

The full machine-readable audit is
`02_low_rank_xm/raw_data/xm_derivative_structure.json`.

