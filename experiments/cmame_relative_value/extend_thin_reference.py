"""Extend only the least-converged M3 thickness reference to n=128."""
from __future__ import annotations

import csv
from .run_scan import DATA, one, save_rows


def main() -> None:
    path = DATA/"reference_convergence.csv"
    rows = list(csv.DictReader(path.open()))
    for method in ("MidedgeTan", "MITC4"):
        if any(r["benchmark"] == "M3" and float(r["t_over_t0"]) == 0.0625
               and int(r["n"]) == 128 and r["method"] == method for r in rows):
            continue
        row = one("M3", 128, .0625, method)
        print(method, row["response"], row["total_wall_s"], flush=True)
        rows.append(row)
        save_rows(path, rows)


if __name__ == "__main__": main()
