"""Verify optional-thickness MidedgeTan binary reproduces the old adapter."""
from __future__ import annotations

import json

from experiments.cmame_round4.common import midedgetan_row
from .run_scan import DATA, T0, midedge


def main() -> None:
    rows = []
    for benchmark, n in (("plate", 16), ("roof", 16), ("M3", 12)):
        old = midedgetan_row(benchmark, n)
        new = midedge(benchmark, n, {"plate": 0.01, **T0}[benchmark])
        error = abs(float(new["observable"])-float(old["observable"]))
        rows.append({"benchmark": benchmark, "n": n,
                     "old_observable": old["observable"],
                     "new_observable": new["observable"],
                     "absolute_difference": error,
                     "same_state_counts": (old["raw_scalar_states"] == new["raw_scalar_states"]
                                           and old["active_states"] == new["active_states"])})
        if error > 1e-10 or not rows[-1]["same_state_counts"]:
            raise AssertionError(rows[-1])
    DATA.mkdir(exist_ok=True)
    (DATA/"midedgetan_gate.json").write_text(json.dumps(rows, indent=2)+"\n")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__": main()
