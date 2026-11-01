# MITC4 baseline (OpenSeesPy ShellMITC4)

| benchmark | n | MITC4 response | error | raw states | active states |
|---|---:|---:|---:|---:|---:|
| plate | 8 | 0.00044188424 | 0.389% | 486 | 451 |
| plate | 16 | 0.00044393098 | 0.073% | 1734 | 1667 |
| plate | 32 | 0.00044495292 | 0.303% | 6534 | 6403 |
| plate | 64 | 0.00044598195 | 0.535% | 25350 | 25091 |
| plate_thin (t=0.001) | 32 | 0.44348481 | 0.028% | 6534 | 6403 |
| plate_thin (t=0.001) | 64 | 0.44360286 | 0.001% | 25350 | 25091 |
| plate_thin (t=0.0001) | 32 | 443.46935 | 0.031% | 6534 | 6403 |
| plate_thin (t=0.0001) | 64 | 443.57433 | 0.008% | 25350 | 25091 |
| roof | 12 | 0.30151202 | 0.294% / 0.306% | 1014 | 949 |
| roof | 16 | 0.30089602 | 0.497% / 0.101% | 1734 | 1649 |
| roof | 24 | 0.30084107 | 0.516% / 0.083% | 3750 | 3625 |
| roof | 32 | 0.30094654 | 0.481% / 0.118% | 6534 | 6369 |
| M3 (raw) | 6 | 0.03798935 | 59.379% | 294 | 251 |
| M3 (raw) | 8 | 0.062858423 | 32.787% | 486 | 431 |
| M3 (raw) | 12 | 0.084657467 | 9.478% | 1014 | 935 |
| M3 (raw) | 16 | 0.090141058 | 3.615% | 1734 | 1631 |
| M3 (raw) | 24 | 0.092551032 | 1.038% | 3750 | 3599 |
| M3 (raw) | 32 | 0.093080467 | 0.472% | 6534 | 6335 |
| M3 (raw) | 48 | 0.093376617 | 0.155% | 14406 | 14111 |
| M3 (raw) | 64 | 0.093472416 | 0.053% | 25350 | 24959 |
| M3 (smoothed) | 6 | 0.038292532 | 59.055% | 294 | 251 |
| M3 (smoothed) | 8 | 0.062905775 | 32.737% | 486 | 431 |
| M3 (smoothed) | 12 | 0.084517454 | 9.628% | 1014 | 935 |
| M3 (smoothed) | 16 | 0.089997072 | 3.769% | 1734 | 1631 |
| M3 (smoothed) | 24 | 0.092451531 | 1.144% | 3750 | 3599 |
| M3 (smoothed) | 32 | 0.093011709 | 0.545% | 6534 | 6335 |

Roof errors are relative to the historical `0.3024` and modern `0.3005924566` references, respectively.

MITC4 is a Mindlin-Reissner shell with 6 DOF per node; its plate response contains a small transverse-shear contribution that the Kirchhoff--Love Checkerboard/MidedgeTan models do not have.
