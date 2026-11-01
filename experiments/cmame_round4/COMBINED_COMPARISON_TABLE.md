# Combined Checkerboard / MidedgeTan / MITC4 comparison

MITC4: OpenSeesPy `ShellMITC4` with `ElasticMembranePlateSection`, 6 DOF per node.  MidedgeTan: official LibShell `MidedgeAngleTanFormulation`.  Checkerboard: current production energy.  Roof errors are historical / modern.

| benchmark | n | Checkerboard | MidedgeTan | MITC4 | CB / MT / MITC4 raw states |
|---|---:|---:|---:|---:|---:|
| plate | 8 | 5.547% | 4.953% | 0.389% | 243 / 451 / 486 |
| plate | 16 | 1.448% | 1.248% | 0.073% | 867 / 1667 / 1734 |
| plate | 32 | 0.361% | 0.313% | 0.303% | 3267 / 6403 / 6534 |
| plate | 64 | 0.066% | 0.078% | 0.535% | 12675 / 25091 / 25350 |
| roof | 12 | 11.455% / 12.125% | 9.099% / 8.553% | 0.294% / 0.306% | 507 / 963 / 1014 |
| roof | 16 | 6.886% / 7.529% | 5.504% / 4.936% | 0.497% / 0.101% | 867 / 1667 / 1734 |
| roof | 24 | 3.664% / 4.287% | 2.200% / 1.612% | 0.516% / 0.083% | 1875 / 3651 / 3750 |
| roof | 32 | 2.541% / 3.157% | 0.856% / 0.260% | 0.481% / 0.118% | 3267 / 6403 / 6534 |
| M3 | 6 | 43.301% / 42.229% | 43.122% | 59.379% | 147 / 267 / 294 |
| M3 | 8 | 16.816% / 15.847% | 13.215% | 32.787% | 243 / 451 / 486 |
| M3 | 12 | 2.794% / 2.201% | 3.741% | 9.478% | 507 / 963 / 1014 |
| M3 | 16 | 0.479% / 0.110% | 5.009% | 3.615% | 867 / 1667 / 1734 |
| M3 | 24 | — | — | 1.038% | — / — / 3750 |
| M3 | 32 | — | — | 0.472% | — / — / 6534 |
| M3 | 48 | — | — | 0.155% | — / — / 14406 |
| M3 | 64 | — | — | 0.053% | — / — / 25350 |

M3 Checkerboard column: raw equator vertex / physical-midpoint reconstruction.  MidedgeTan and MITC4 use the raw equator vertex.  MITC4 M3 n=24,32,48,64 are convergence points beyond the manuscript range.

MITC4 is a Mindlin-Reissner shell; at t/a=0.01 its plate response has a ~0.5% transverse-shear floor relative to the thin-plate reference.  At t/a=0.001 the error drops to 0.03%/0.001%, confirming the source.
