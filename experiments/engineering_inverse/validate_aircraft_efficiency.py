"""Audit all 31 optimized layouts for exact-strategy equivalence."""
from __future__ import annotations

import json
import numpy as np

from experiments.engineering_inverse.aircraft_skin import OUT


def run():
    report={}
    for mesh in ("12x24","16x32","20x40"):
        records={m:json.loads((OUT/f"aircraft_efficiency_{mesh}_{m}_0.json").read_text())
                 for m in ("cb_direct","cb_shared","cb_shared_sparse",
                           "cb_shared_sparse_legacy",
                           "mitc4_direct","mitc4_shared")}
        checks={}
        for left,right in (("cb_direct","cb_shared"),
                           ("cb_shared","cb_shared_sparse"),
                           ("cb_shared_sparse_legacy","cb_shared_sparse"),
                           ("mitc4_direct","mitc4_shared")):
            maxima=dict(stroke_m=0.,reaction_N=0.,rms_error_mm=0.,
                        max_error_mm=0.,objective=0.)
            for a,b in zip(records[left]["layouts"],records[right]["layouts"]):
                if a["stations"]!=b["stations"]:raise AssertionError("layout order")
                for ca,cb in zip(a["cases"],b["cases"]):
                    maxima["stroke_m"]=max(maxima["stroke_m"],
                        float(np.max(abs(np.asarray(ca["strokes_m"])-cb["strokes_m"]))))
                    maxima["reaction_N"]=max(maxima["reaction_N"],
                        float(np.max(abs(np.asarray(ca["reactions_N"])-cb["reactions_N"]))))
                    for key in ("rms_error_mm","max_error_mm","objective"):
                        maxima[key]=max(maxima[key],abs(ca[key]-cb[key]))
            checks[f"{left}_vs_{right}"]=maxima
            if maxima["max_error_mm"]>1e-3 or maxima["reaction_N"]>1e-2:
                raise AssertionError((mesh,left,right,maxima))
        report[mesh]=checks
    (OUT/"aircraft_efficiency_equivalence.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps(report,indent=2))
    return report


if __name__=="__main__":run()
