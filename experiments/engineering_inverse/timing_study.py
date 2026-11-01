"""Three repeated full linear-design loops, with initialization included."""

import json
import platform
import time

import numpy as np

from experiments.engineering_inverse.aircraft_skin import OUT,run as aircraft_run
from experiments.engineering_inverse.automotive_panel import run as automotive_run


def run(repeats=3):
    rows=[]
    preserved={name:(OUT/name).read_bytes() for name in
               ("aircraft_12x24.json","automotive_12x8.json")}
    for case in ("aircraft","automotive"):
        for repeat in range(repeats):
            start=time.perf_counter()
            result=(aircraft_run(12,24,True) if case=="aircraft" else
                    automotive_run(12,8,True,False))
            rows.append(dict(case=case,repeat=repeat,
                observed_wall_seconds=time.perf_counter()-start,
                checkerboard_seconds=(result["checkerboard"]["total_seconds"]),
                mitc4_seconds=(result["mitc4"]["search_seconds"] if case=="aircraft"
                                 else result["mitc4"]["total_seconds"])))
    for name,contents in preserved.items():
        (OUT/name).write_bytes(contents)
    summary={case:{method:float(np.median([r[method] for r in rows
                                          if r["case"]==case]))
                   for method in ("checkerboard_seconds","mitc4_seconds")}
             for case in ("aircraft","automotive")}
    data=dict(cpu=platform.processor(),threads="OPENBLAS_NUM_THREADS=1 recommended",
              repeats=repeats,rows=rows,median=summary,
              timing_scope="aircraft: complete 31-layout linear screen; automotive: linear design including response basis; excludes nonlinear audit")
    (OUT/"timing_repeats.json").write_text(json.dumps(data,indent=2)+"\n")
    return data


if __name__=="__main__":
    print(json.dumps(run()["median"],indent=2))
