#!/usr/bin/env python3
# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Generate per-stage axis SVGs and composite overview SVGs for all preset geometries.

Usage (from the repository root)::

    python tools/generate_geometry_drawings.py

Output directory: ``docs/source/_static/geometries/<geometry>/``

Each geometry gets:
- One SVG per stage (``<stage>.svg``)
- One composite overview (``<geometry>_all.svg``)

Requires matplotlib (install with ``pip install ad_hoc_diffractometer[doc]``).
"""

from __future__ import annotations

import os

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt

from ad_hoc_diffractometer.drawing import draw_geometry_axes
from ad_hoc_diffractometer.drawing import draw_stage_axis
from ad_hoc_diffractometer.factories import list_geometries


def main() -> None:
    """Generate all geometry drawing SVGs."""
    base = os.path.join("docs", "source", "_static", "geometries")

    total = 0
    for factory in list_geometries().values():
        g = factory()
        gdir = os.path.join(base, g.name)
        os.makedirs(gdir, exist_ok=True)

        all_stages = list(g.sample_stages) + list(g.detector_stages)
        for stage in all_stages:
            fig = draw_stage_axis(stage, g.basis)
            path = os.path.join(gdir, f"{stage.name}.svg")
            fig.savefig(path, bbox_inches="tight")
            plt.close(fig)
            total += 1

        fig2 = draw_geometry_axes(g)
        path2 = os.path.join(gdir, f"{g.name}_all.svg")
        fig2.savefig(path2, bbox_inches="tight")
        plt.close(fig2)
        total += 1

        print(f"  {g.name}: {len(all_stages)} stages + 1 composite")

    print(f"\nGenerated {total} SVGs in {base}/")


if __name__ == "__main__":
    main()
