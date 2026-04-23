#!/usr/bin/env python3
# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Generate interactive Plotly HTML figures and fallback SVG images for all
preset geometries.

Usage (from the repository root)::

    python tools/generate_geometry_drawings.py

Output directory: ``docs/source/_static/geometries/<geometry>/``

Each geometry gets:
- One self-contained HTML file (``<geometry>.html``) containing an
  interactive :class:`~ad_hoc_diffractometer.drawing.GeometryAxisFigure`
  with all stages arranged in their parent-child tree.
- One static SVG fallback (``<geometry>.svg``) for environments where
  WebGL is unavailable (remote X sessions, PDF builds, etc.).

Requires plotly (``pip install plotly``) and kaleido (``pip install
kaleido``) for SVG export.
"""

from __future__ import annotations

import os

from ad_hoc_diffractometer.drawing import GeometryAxisFigure
from ad_hoc_diffractometer.factories import list_geometries


def main() -> None:
    """Generate all geometry Plotly HTML and SVG files."""
    base = os.path.join("docs", "source", "_static", "geometries")

    total_html = 0
    total_svg = 0
    for name, _factory in list_geometries().items():
        gdir = os.path.join(base, name)
        os.makedirs(gdir, exist_ok=True)

        fig = GeometryAxisFigure(name)

        # Interactive HTML
        html_path = os.path.join(gdir, f"{name}.html")
        fig.write_html(
            html_path,
            full_html=True,
            include_plotlyjs="cdn",
            config={"responsive": True},
        )
        total_html += 1

        # Static SVG fallback
        svg_path = os.path.join(gdir, f"{name}.svg")
        fig.write_image(svg_path)
        total_svg += 1

        print(f"  {name}: {html_path}  +  {svg_path}")

    print(f"\nGenerated {total_html} HTML + {total_svg} SVG files in {base}/")


if __name__ == "__main__":
    main()
