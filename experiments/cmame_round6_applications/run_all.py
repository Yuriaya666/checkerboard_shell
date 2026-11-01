"""Run both application studies and regenerate figures."""

from __future__ import annotations

import json

from experiments.cmame_round6_applications import (
    aircraft_skin, automotive_panel, make_figures,
)


def main() -> dict:
    aircraft = aircraft_skin.run()
    automotive = automotive_panel.run()
    make_figures.main()
    return {
        "aircraft": aircraft["checkerboard"],
        "automotive": automotive["checkerboard"],
    }


if __name__ == "__main__":
    print(json.dumps(main(), indent=2))
