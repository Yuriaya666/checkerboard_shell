# External baseline status

The historical source contains verified wrappers around official LibShell
commit `6072d9e4b3c9f99c8cfe538b8b5560954d87dfdf`.  The
`MidedgeAngleTanFormulation` binaries were found and rerun for plate,
Scordelis--Lo, and M3.  Their executable hashes and raw JSON rows are stored
under `00_driver_recovery/raw_data/`.

Round 2 also adds a geometry-input adapter around that same unmodified
`MidedgeAngleTanFormulation` for the surface-preserving distortion grids.  The
adapter reproduces all three uniform historical values (largest relative
difference `2.0e-10` from ASCII vertex transport/separate compilation), and
its source, build log, binary hash, and 30 raw rows are preserved.  This is
benchmark glue, not a replacement shell formulation.

No provenance-closed MITC4 or other verified four-node assumed-strain shell
exists in either repository or the verified dependency checkout.  Round 2
therefore does not introduce a homemade substitute.  A later external study
should use a documented Simo-type MITC4 implementation (commercial or a
well-validated open-source FE code) with the same geometry, loads, physical
BC, mesh family, and reference observables.
