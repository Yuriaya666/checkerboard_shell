# Unified CMAME comparison (Round 4)

Checkerboard plate rows use the complete 3V production solve. MidedgeTan rows are rerun through one adapter on the same uniform logical grids. Roof errors are given against both references.

| benchmark | n | Checkerboard response (error) | MidedgeTan response (error) | CB raw/active | MT raw/active |
|---|---:|---:|---:|---:|---:|
| plate | 8 | 0.00046821809 (5.5475%) | 0.00046557919 (4.9526%) | 243/208 | 451/257 |
| plate | 16 | 0.00045003181 (1.4479%) | 0.00044914514 (1.2480%) | 867/800 | 1667/1025 |
| plate | 32 | 0.00044520816 (0.3605%) | 0.00044499585 (0.3126%) | 3267/3136 | 6403/4097 |
| plate | 64 | 0.00044390349 (0.0664%) | 0.00044395583 (0.0782%) | 12675/12416 | 25091/16385 |
| roof | 12 | 0.33703868 (11.455% / 12.125%) | 0.27488381 (9.099% / 8.553%) | 507/468 | 963/912 |
| roof | 16 | 0.32322295 (6.886% / 7.529%) | 0.28575551 (5.504% / 4.936%) | 867/816 | 1667/1600 |
| roof | 24 | 0.31347934 (3.664% / 4.287%) | 0.29574826 (2.200% / 1.612%) | 1875/1800 | 3651/3552 |
| roof | 32 | 0.31008252 (2.541% / 3.157%) | 0.29981084 (0.856% / 0.260%) | 3267/3168 | 6403/6272 |
| M3 | 6 | 0.053025732 (43.3011% / 42.2295%) | 0.053193383 (43.1218%) | 147/107 | 267/212 |
| M3 | 8 | 0.077795323 (16.8156% / 15.8473%) | 0.081163021 (13.2146%) | 243/191 | 451/380 |
| M3 | 12 | 0.090908217 (2.7944% / 2.2015%) | 0.097019751 (3.7405%) | 507/431 | 963/860 |
| M3 | 16 | 0.093073115 (0.4795% / 0.1105%) | 0.098206305 (5.0093%) | 867/767 | 1667/1532 |

For roof, the two errors are relative to the historical `0.3024` and the modern `0.3005924566` reference, respectively.

For M3, the two errors are the raw equator vertex (main, the same discrete observable as MidedgeTan) and the physical-midpoint moment reconstruction (in parentheses), respectively.
