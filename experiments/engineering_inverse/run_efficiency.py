"""Reproduce the strict efficiency study from repository root.

Run with PYTHONPATH=src:. and one OpenBLAS/OpenMP/MKL thread.  The aircraft
benchmark uses fresh processes per timed repeat.  This entry deliberately
does not modify either engineering model or the historical result files.
"""
from __future__ import annotations

import argparse
import os
import subprocess
import sys


def run():
    parser=argparse.ArgumentParser()
    parser.add_argument("--skip-aircraft",action="store_true")
    parser.add_argument("--skip-automotive",action="store_true")
    parser.add_argument("--repeats",type=int,default=3)
    parser.add_argument("--include-feasible-start",action="store_true",
                        help="also run supplemental three-method comparison from saved CB design")
    args=parser.parse_args()
    env=os.environ.copy()
    for variable in ("OPENBLAS_NUM_THREADS","OMP_NUM_THREADS","MKL_NUM_THREADS"):
        env[variable]="1"
    env["PYTHONPATH"]="src:."
    def stage(module,*arguments):
        subprocess.run([sys.executable,"-m",module,*arguments],check=True,env=env)
    stage("experiments.engineering_inverse.efficiency_environment")
    if not args.skip_aircraft:
        stage("experiments.engineering_inverse.efficiency_aircraft",
              "--repeats",str(args.repeats))
        stage("experiments.engineering_inverse.validate_aircraft_efficiency")
    if not args.skip_automotive:
        stage("experiments.engineering_inverse.validate_automotive_steps")
        stage("experiments.engineering_inverse.efficiency_forward")
        starts=("zero","cb_feasible") if args.include_feasible_start else ("zero",)
        for start in starts:
            for method in ("cb_consistent","cb_fd","asd_fd"):
                stage("experiments.engineering_inverse.efficiency_automotive",
                      method,"--start",start)
            stage("experiments.engineering_inverse.efficiency_automotive_reference",
                  "--start",start)
        stage("experiments.engineering_inverse.validate_automotive_stationarity")
        stage("experiments.engineering_inverse.automotive_branch_probe")
    stage("experiments.engineering_inverse.make_efficiency_figures")


if __name__=="__main__":run()
