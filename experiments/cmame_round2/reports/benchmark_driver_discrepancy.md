# Benchmark-driver discrepancy audit

## Recovery result

The manuscript drivers were found in the clean historical repository
`/home/song/projects/CB-shell` at tag `v1.0-arxiv` (HEAD `093ad89`).  The
benchmark modules and their outputs are tracked Git files, not reconstructed
from target values.  They were ported to the current package and rerun through
the new sparse/local plus rank-one path.  The complete rerun is in
`00_driver_recovery/raw_data/driver_recovery.json`.

The rerun reproduces all requested manuscript Checkerboard values:

| benchmark | n | recovered response | recovered error |
|---|---:|---:|---:|
| plate | 8,16,32,64 | 4.68433868e-4, 4.50210236e-4, 4.45377092e-4, 4.44069708e-4 | 5.5961%, 1.4881%, 0.3986%, 0.1039% |
| roof | 12,16,24,32 | 0.33703868, 0.32322295, 0.31347934, 0.31008252 | 11.4546%, 6.8859%, 3.6638%, 2.5405% |
| M3 | 6,8,12,16 | 0.054027889, 0.078700913, 0.091462717, 0.093418247 | 42.2295%, 15.8473%, 2.2015%, 0.11046% |

## Plate

| item | manuscript driver | Round-1 reconstructed driver | consequence |
|---|---|---|---|
| geometry | unit flat square, n by n logical cells | same continuum square | no material difference at zero distortion |
| unknown branch | scalar linear transverse branch | full 3V reference tangent with in-plane anchors | transverse branch is decoupled, so this is not the main discrepancy |
| load | `pressure/n^2` on every active interior transverse state | physical-area lumping to raw nodes | identical on the active uniform flat grid |
| BC | `w=0` on all four edges | same transverse BC plus in-plane rigid anchors | no transverse difference |
| observable | raw transverse center vertex | five-node gauge-annihilating smoothed center | **main discrepancy** |
| mesh convention | n cells, n+1 raw vertices | same | none |

For n=8, the manuscript raw center value is `4.684338681969e-4`, giving
5.5961%.  Round 1 reported the smoothed center value `4.505145594119e-4`,
giving 1.5567%.  No material, load, beta, or reference value changed.  The
apparently better Round-1 number was caused by measuring a different discrete
observable.

## Scordelis--Lo roof

| item | manuscript driver | Round-1 reconstructed driver |
|---|---|---|
| geometry | half-length 25 and full 80-degree arc on the same R=25 cylinder | full length 50, with two physical ends |
| parameter order | columns = longitudinal x; rows = cylinder arc | columns = arc; rows = full longitudinal coordinate |
| longitudinal BC | symmetry `u_x=0` at midspan x=0 | no half-model symmetry plane |
| diaphragm | at x=25, `u_y=u_z=0` | `u_x=u_z=0` at both full-length ends; not the manuscript diaphragm components |
| load | tensor-product trapezoidal integration of global `-z q=90` in the isometric arc/length parameters | polygon-area lumping on the reconstructed full model |
| observable | raw z displacement at the positive-40-degree free-edge node on the longitudinal symmetry line | smoothed three-node free-edge midpoint observable |
| reference | 0.3024 primary; 0.3005924566 retained separately | 0.3024 |

These are different boundary-value problems.  The Round-1 n=12 value
0.32663423 cannot be compared with the manuscript n=12 value 0.33703868.
The discrepancy is driven by the full/half model and diaphragm/symmetry
conditions, with an additional observable difference.

## MacNeal--Harder M3

| item | manuscript driver | Round-1 reconstructed driver |
|---|---|---|
| geometry | quarter R=10 sphere, polar angle 18 to 90 degrees and azimuth 0 to 90 degrees | same continuum sphere sector |
| spacing variables | constant angular `hx=dphi`, `hy=dtheta` | arc-length-scaled parameters; continuum-equivalent on the uniform sphere but not byte-identical |
| symmetry | physical theta-edge midpoint reflection constraints plus second-order even-component derivative constraints | raw-node component pins only |
| gauge | three explicit alternating gauge-slice equations, separate from mechanical BC | gauge was implicitly/partially selected by raw pins |
| rigid mode | physical mean-z equation | one raw z anchor |
| point A/B | third-order equator extrapolation `(15/8,-5/4,3/8)` of physical theta-edge midpoints | incident-edge corner average |
| load | `+1 e_x` at physical A and `-1 e_y` at physical B (the quarter loads for canonical P=2) | signs were reversed (`-A_x,+B_y`) |
| observable | `-Delta y_B` at extrapolated physical B | the Round-1 code measured `-A_x`, despite labeling it B |

The last three rows are decisive.  Round 1 n=6 gave 0.08707038, while the
recovered driver gives 0.054027889.  This was not a formulation change: the
Round-1 problem used different symmetry constraints, different point
operators, reversed load signs, and the wrong labeled output point.

## Indexing, midpoint interpolation, and quadrature conclusion

- All manuscript mesh labels are numbers of logical cells per direction;
  there is no hidden n-versus-n+1 convention.
- Plate mechanics is evaluated in the raw scalar transverse branch, while the
  manuscript figure may reconstruct the gauge-fixed physical surface only for
  visualization.
- Roof and M3 manuscript point operators are deliberately tied to their
  physical benchmark definitions.  Replacing them by generic midpoint
  averages changes the discrete observable and, for M3, the applied load.
- The recovered drivers use the accepted moment boundary closure and complete
  interior `X_W^rel` policy unchanged.

