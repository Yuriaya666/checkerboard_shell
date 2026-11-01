"""Small refinement and sensitivity cross-check on the Pq diagnostic.

Run `python -m experiments.cmame_positioning.diagnostic` first. Its n=8/12
rows are reused; all additional solves below are saved incrementally.
"""
from __future__ import annotations

import csv
import json
from .diagnostic import (
    RESULTS, cb_solve, grid, loaded_nodes, midpoint_map,
    mitc4_solve, physical_fields,
)


def save(rows):
    RESULTS.mkdir(exist_ok=True)
    path=RESULTS/"refinement.json"
    path.write_text(json.dumps(rows,indent=2)+"\n")
    flat=[]
    for row in rows:
        flat.append({k:v for k,v in row.items() if k!="payload"}|row["payload"])
    with (RESULTS/"refinement.csv").open("w",newline="") as file:
        keys=list(dict.fromkeys(k for row in flat for k in row))
        writer=csv.DictWriter(file,fieldnames=keys)
        writer.writeheader();writer.writerows(flat)


def main():
    low=json.loads((RESULTS/"diagnostic.json").read_text())
    rows=[]
    for r in low:
        for method,key in (("CB-Full","CB"),("ShellMITC4(Pq)","MITC4_Pq")):
            rows.append({"n":r["n"],"parameter":r["parameter"],
                         "method":method,"payload":r[key],
                         "provenance":"diagnostic.json"})
    save(rows)

    for n in (16,24,32):
        P,faces,_,_=midpoint_map(n)
        for parameter in (-.1,0.,.1):
            q=grid(n,parameter)
            y,_,_=physical_fields(q,P,faces,n)
            support,loaded,force=loaded_nodes(y,n)
            payload=cb_solve(q,P,faces,n,support,loaded,force)
            rows.append({"n":n,"parameter":parameter,"method":"CB-Full",
                         "payload":payload,"provenance":"fresh run"})
            save(rows)
            print(n,parameter,"CB-Full",payload["tip_z_m"],flush=True)

    external_cases={
        "ShellMITC4": {16:(0.,),24:(0.,),32:(0.,),48:(0.,),
                        64:(-.1,0.,.1),80:(0.,)},
        "ASDShellQ4": {8:(0.,),16:(0.,),32:(0.,),48:(0.,),
                        64:(-.1,0.,.1),80:(0.,)},
    }
    for element,cases in external_cases.items():
        for n,parameters in cases.items():
            P,faces,_,_=midpoint_map(n)
            for parameter in parameters:
                q=grid(n,parameter)
                y,_,_=physical_fields(q,P,faces,n)
                support,loaded,force=loaded_nodes(y,n)
                payload=mitc4_solve(y,faces,support,loaded,force,element=element)
                rows.append({"n":n,"parameter":parameter,
                             "method":element+"(Pq)",
                             "payload":payload,"provenance":"fresh run"})
                save(rows)
                print(n,parameter,element,payload["tip_z_m"],flush=True)

    summary=[]
    for method in ("CB-Full","ShellMITC4(Pq)","ASDShellQ4(Pq)"):
        for n in sorted({r["n"] for r in rows if r["method"]==method}):
            vals={r["parameter"]:r["payload"]["tip_z_m"] for r in rows
                  if r["method"]==method and r["n"]==n}
            summary.append({"method":method,"n":n,"response_at_zero_m":vals.get(0.),
                            "central_design_sensitivity_m_per_unit":
                            (vals[.1]-vals[-.1])/.2 if -.1 in vals and .1 in vals else None,
                            "reference_geometry":"identical physical midpoint nodes Pq"})
    (RESULTS/"summary.json").write_text(json.dumps(summary,indent=2)+"\n")
    with (RESULTS/"summary.csv").open("w",newline="") as file:
        writer=csv.DictWriter(file,fieldnames=list(summary[0]))
        writer.writeheader();writer.writerows(summary)


if __name__=="__main__":main()
