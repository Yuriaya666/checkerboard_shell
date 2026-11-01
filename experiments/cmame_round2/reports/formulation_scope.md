# Structured-grid and parameter-spacing scope

The current formulation requires a structured logical `(i,j)` topology and
constant parameter steps `hx,hy` within one discretization.  It supports
non-rectangular physical quads and curved embeddings obtained by sampling the
same smooth surface at redistributed logical nodes, but the present formulas
do **not** constitute a variable-step finite-difference derivation.

Formula-by-formula:

1. The B diagonals and their normalization use the two fixed logical
   directions `d1=(hx,hy)`, `d2=(-hx,hy)` and one fixed
   `ell=sqrt(hx^2+hy^2)`.  A variable grid would require face- or edge-local
   direction maps, not a single substituted scalar pair.
2. W normal differences use centered differences divided by the same fixed
   `ell`.  On unequal adjacent intervals, the centered coefficients and
   evaluation location change; replacing `h` locally without rederiving the
   weights is not second-order equivalent.
3. The boundary closure uses moment-exact samples at fixed cell-center
   locations 1/2, 3/2, 5/2, producing `(15/8,-5/4,3/8)` and
   `(-2,3,-1)/h`.  Nonuniform centers require new Vandermonde weights and a
   consistent corner tensor product.
4. `gamma_W` contains the fixed aspect-ratio factors `hy/hx` and `hx/hy`
   obtained by exact auxiliary-Q1 integration on one reference patch.
   Variable steps require patch-local integration using its four actual
   interval lengths.
5. `X_M` is a global alternating character on the structured lattice and its
   area weights may vary geometrically, but the connection differences are
   still the fixed logical mixed-difference operator.  Area weighting alone
   does not turn it into a variable-step construction.
6. `X_W^rel` uses the same fixed mixed reconstruction and is accumulated only
   on complete interior W patches.  That patch definition remains meaningful
   on a structured topology, but its variable-step coefficient has not been
   derived.

Therefore the Round-1 geometric-spacing-ratio study is marked
`OUTSIDE CURRENT FORMULATION SCOPE`.  In particular, the ratio-8 plate error
must not be presented as a failure of the claimed constant-spacing method.
A non-empirical extension is possible in principle by deriving each B/W
operator from actual local parameter coordinates, recomputing moment weights,
and re-integrating the local Q1 coefficient, but that is a formulation change
and is only a future-work proposal here.

