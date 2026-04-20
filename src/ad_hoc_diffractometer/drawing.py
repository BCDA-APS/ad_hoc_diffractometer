# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
drawing.py — geometry visualisation helpers.

Provides functions to generate:

1. **Graphviz DOT source** for stage coupling diagrams
   (:func:`geometry_dot`) — shows parent-child relationships, roles, and
   axis labels.

2. **Matplotlib per-stage axis diagrams** (:func:`draw_stage_axis`) —
   3D view of the lab-frame basis vectors with the stage's rotation axis
   highlighted.

3. **Matplotlib composite diagram** (:func:`draw_geometry_axes`) — all
   stages of a geometry in a single figure.

The DOT source is designed to be pasted into a Sphinx ``.. graphviz::``
directive in the geometry documentation pages.

Matplotlib is an optional dependency (listed in ``pyproject.toml``
``[project.optional-dependencies] doc``).  Functions that require it
raise :class:`ImportError` with a helpful message if it is not installed.
"""

from __future__ import annotations

import numpy as np

from .axes import axis_label


def _physical_label(axis_vec: np.ndarray, basis: dict | None) -> str:
    """Return a human-readable axis label with physical direction if possible.

    For standard signed basis vectors, returns e.g. ``"-lateral"`` or
    ``"+vertical"``.  For tilted axes (kappa), returns the numeric label.
    """
    if basis is None:
        return axis_label(axis_vec)

    atol = 1e-8
    for direction, bvec in basis.items():
        bvec = np.asarray(bvec, dtype=float)
        if np.allclose(axis_vec, bvec, atol=atol):
            return f"+{direction}"
        if np.allclose(axis_vec, -bvec, atol=atol):
            return f"-{direction}"
    # Tilted axis (kappa) — fall back to numeric label
    return axis_label(axis_vec)


def _handedness(axis_vec: np.ndarray, basis: dict | None) -> str:
    """Return 'RH' or 'LH' for the handedness of a stage axis.

    A positive basis vector means right-handed rotation; a negated basis
    vector means left-handed.  Non-standard axes (kappa) are always
    right-handed (the tilt encodes the direction).
    """
    if basis is None:
        return ""
    atol = 1e-8
    for bvec in basis.values():
        bvec = np.asarray(bvec, dtype=float)
        if np.allclose(axis_vec, bvec, atol=atol):
            return "RH"
        if np.allclose(axis_vec, -bvec, atol=atol):
            return "LH"
    return "RH"  # tilted axes are right-handed by convention


# ---------------------------------------------------------------------------
# Graphviz DOT generation
# ---------------------------------------------------------------------------


def geometry_dot(geometry) -> str:
    """Generate Graphviz DOT source for a stage coupling diagram.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        A fully constructed geometry instance.

    Returns
    -------
    str
        DOT source string suitable for a Sphinx ``.. graphviz::`` directive
        or ``graphviz.Source()`` rendering.

    Examples
    --------
    >>> from ad_hoc_diffractometer.presets import fourcv
    >>> from ad_hoc_diffractometer.drawing import geometry_dot
    >>> print(geometry_dot(fourcv()))
    digraph fourcv {
        ...
    }
    """
    basis = getattr(geometry, "basis", None)
    lines = [f"digraph {geometry.name} {{"]
    lines.append("    rankdir=BT;")
    lines.append(f'    label="{geometry.name}";')
    lines.append("    labelloc=t;")
    lines.append("    fontsize=14;")
    lines.append("    node [shape=box, style=filled, fontsize=11];")
    lines.append("")

    # Collect all stages with their properties
    all_stages = list(geometry.sample_stages) + list(geometry.detector_stages)

    for stage in all_stages:
        phys = _physical_label(stage.axis, basis)
        hand = _handedness(stage.axis, basis)
        role = stage.role
        color = "#a8d8ea" if role == "sample" else "#f8a5a5"

        label_parts = [stage.name, f"axis: {phys}"]
        if hand:
            label_parts.append(hand)

        node_label = "\\n".join(label_parts)
        lines.append(f'    {stage.name} [label="{node_label}", fillcolor="{color}"];')

    lines.append("")

    # Force root nodes (parent=None) to the same rank (bottom in BT layout)
    roots = [s.name for s in all_stages if s.parent is None]
    if len(roots) > 1:
        lines.append(f"    {{ rank=same; {'; '.join(roots)}; }}")
        lines.append("")

    # Edges: child → parent (arrow points down toward the base in BT layout)
    for stage in all_stages:
        if stage.parent is not None:
            lines.append(f"    {stage.name} -> {stage.parent};")

    # Add a legend (vertical layout)
    lines.append("")
    lines.append("    // Legend")
    lines.append("    subgraph cluster_legend {")
    lines.append('        label="Legend";')
    lines.append("        fontsize=8;")
    lines.append("        style=dashed;")
    lines.append("        color=gray;")
    lines.append(
        '        sample_legend [label="sample", '
        'fillcolor="#a8d8ea", shape=box, style=filled, fontsize=7];'
    )
    lines.append(
        '        detector_legend [label="detector", '
        'fillcolor="#f8a5a5", shape=box, style=filled, fontsize=7];'
    )
    lines.append("        sample_legend -> detector_legend [style=invis];")
    lines.append("    }")

    lines.append("}")
    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Matplotlib stage axis diagrams
# ---------------------------------------------------------------------------

# View angles (in display coordinates where ver=+z, lon=+y, lat=+x).
# The viewer looks roughly along the beam (lon), slightly offset so all
# three basis vectors are visible.
_VIEW_ELEV = -7.5
_VIEW_AZIM = 110

# Per-direction arc start offsets (radians) so the arc is not edge-on.
_ARC_START_OFFSETS = {
    "vertical": np.radians(90),
    "lateral": np.radians(130),
    "longitudinal": np.radians(90),
}

_BASIS_COLORS = {"vertical": "green", "longitudinal": "blue", "lateral": "red"}


def _require_matplotlib():
    """Import and return matplotlib.pyplot, or raise ImportError."""
    try:
        import matplotlib.pyplot as plt

        return plt
    except ImportError:
        raise ImportError(
            "matplotlib is required for drawing functions.  "
            "Install it with:  pip install ad_hoc_diffractometer[doc]"
        ) from None


def _display_rotation(basis: dict) -> np.ndarray:
    """Build a 3x3 rotation matrix that maps basis vectors to display coords.

    Display coordinates: lat→+x, lon→+y, ver→+z.  This ensures the same
    physical view (ver up, lon away, lat left) regardless of which
    Cartesian axes the basis maps to.

    Parameters
    ----------
    basis : dict
        Basis dictionary with ``"vertical"``, ``"longitudinal"``,
        ``"lateral"`` keys.

    Returns
    -------
    np.ndarray, shape (3, 3)
        Rotation matrix R such that ``R @ v`` maps a vector from the
        Cartesian frame into display coordinates.
    """
    lat = np.asarray(basis["lateral"], dtype=float)
    lon = np.asarray(basis["longitudinal"], dtype=float)
    ver = np.asarray(basis["vertical"], dtype=float)
    return np.array([lat, lon, ver])


def _stage_draw_direction(stage_axis: np.ndarray, basis: dict) -> tuple:
    """Determine the drawing direction for a stage axis.

    Negative-signed axes (e.g. -lateral) are drawn in the positive basis
    direction for visibility; the sign is conveyed by the subtitle and
    arc labels.

    Returns
    -------
    tuple of (draw_direction, axis_type, is_negated)
        draw_direction : np.ndarray — unit vector for drawing
        axis_type : str or None — "vertical", "lateral", "longitudinal", or None
        is_negated : bool — True if the stage axis is the negation of a basis vector
    """
    atol = 1e-8
    for direction, bvec in basis.items():
        bvec = np.asarray(bvec, dtype=float)
        if np.allclose(stage_axis, bvec, atol=atol):
            return bvec / np.linalg.norm(bvec), direction, False
        if np.allclose(stage_axis, -bvec, atol=atol):
            return bvec / np.linalg.norm(bvec), direction, True
    # Non-standard axis (kappa) — draw as-is
    return stage_axis / np.linalg.norm(stage_axis), None, False


def _subtitle(stage_axis: np.ndarray, basis: dict) -> str:
    """Build the subtitle string for a stage diagram."""
    phys = _physical_label(stage_axis, basis)
    cartesian = axis_label(stage_axis)
    if phys != cartesian:
        return f"{phys}  ({cartesian})"
    return f"[{stage_axis[0]:.3f}, {stage_axis[1]:.3f}, {stage_axis[2]:.3f}]"


def _draw_stage_on_axes(
    ax,
    stage,
    basis: dict,
    R: np.ndarray,
    *,
    geometry_name: str = "",
    basis_length: float = 1.3,
    axis_length: float = 0.7,
    arc_radius: float = 0.18,
    arc_center_t: float = 0.45,
    arc_sweep_deg: float = 180.0,
    label_fontsize: float = 9.0,
    arc_label_fontsize: float = 8.0,
) -> None:
    """Draw a single stage diagram on a matplotlib 3D Axes.

    All data is rotated into display coordinates via R before plotting.
    """

    def to_disp(v):
        return R @ np.asarray(v, dtype=float)

    stage_axis = np.asarray(stage.axis, dtype=float)
    draw_dir, axis_type, _is_negated = _stage_draw_direction(stage_axis, basis)

    # Perpendicular vectors for the rotation arc
    if abs(np.dot(draw_dir, np.array([1.0, 0, 0]))) < 0.9:
        perp1 = np.cross(draw_dir, np.array([1.0, 0, 0]))
    else:
        perp1 = np.cross(draw_dir, np.array([0, 1.0, 0]))
    perp1 /= np.linalg.norm(perp1)
    perp2 = np.cross(draw_dir, perp1)

    # Arc start angle with per-direction offset
    arc_start = np.radians(90)
    if axis_type in _ARC_START_OFFSETS:
        arc_start += _ARC_START_OFFSETS[axis_type]

    arc_center = draw_dir * arc_center_t
    arc_angles = np.linspace(arc_start, arc_start + np.radians(arc_sweep_deg), 60)
    arc_pts = np.array(
        [
            arc_center + arc_radius * (np.cos(a) * perp1 + np.sin(a) * perp2)
            for a in arc_angles
        ]
    )

    # --- Draw order: arc first (behind), then basis vectors, then stage axis ---

    # Arc (black)
    arc_disp = np.array([to_disp(p) for p in arc_pts])
    ax.plot(
        arc_disp[:, 0],
        arc_disp[:, 1],
        arc_disp[:, 2],
        color="black",
        linewidth=1.2,
    )

    # +/- labels
    arc_center_d = to_disp(arc_center)
    plus_dir = arc_disp[-1] - arc_center_d
    plus_dir /= np.linalg.norm(plus_dir)
    minus_dir = arc_disp[0] - arc_center_d
    minus_dir /= np.linalg.norm(minus_dir)
    ax.text(
        *(arc_disp[-1] + plus_dir * 0.12),
        "+",
        fontsize=arc_label_fontsize,
        fontweight="bold",
        ha="center",
        va="center",
        color="black",
    )
    ax.text(
        *(arc_disp[0] + minus_dir * 0.12),
        "\u2212",
        fontsize=arc_label_fontsize,
        fontweight="bold",
        ha="center",
        va="center",
        color="black",
    )

    # Basis vectors
    for direction, vec in basis.items():
        vec = np.asarray(vec, dtype=float)
        dv = to_disp(vec)
        color = _BASIS_COLORS.get(direction, "gray")
        bv = dv * basis_length
        ax.quiver(
            0,
            0,
            0,
            bv[0],
            bv[1],
            bv[2],
            color=color,
            alpha=0.4,
            arrow_length_ratio=0.08,
            linewidth=1.5,
        )
        tv = dv * (basis_length + 0.15)
        ax.text(
            tv[0],
            tv[1],
            tv[2],
            direction[:3],
            color=color,
            fontsize=label_fontsize,
            ha="center",
        )

    # Stage axis (grey line, on top)
    tip_d = to_disp(draw_dir * axis_length)
    ax.plot(
        [0, tip_d[0]],
        [0, tip_d[1]],
        [0, tip_d[2]],
        color="0.55",
        linewidth=3.0,
        solid_capstyle="round",
    )


def draw_stage_axis(
    stage,
    basis: dict,
    *,
    geometry_name: str = "",
    figsize: tuple[float, float] = (3.5, 2.8),
) -> matplotlib.figure.Figure:  # noqa: F821
    """Draw a 3D diagram for a single stage showing its rotation axis.

    The view is oriented so that the vertical basis vector points up,
    the longitudinal vector points away (into the screen), and the
    lateral vector points to the left.  This convention is independent
    of the Cartesian-axis mapping used by the geometry's basis.

    Parameters
    ----------
    stage : Stage
        The stage to draw.
    basis : dict
        Basis dictionary mapping physical direction names to unit vectors.
    geometry_name : str
        Geometry name for the title (e.g. ``"fourcv"``).
    figsize : tuple
        Figure size in inches.

    Returns
    -------
    matplotlib.figure.Figure
        The figure object.  Call ``fig.savefig(path)`` to write an SVG.
    """
    plt = _require_matplotlib()
    R = _display_rotation(basis)

    fig = plt.figure(figsize=figsize, frameon=False)
    ax = fig.add_subplot(111, projection="3d")

    _draw_stage_on_axes(ax, stage, basis, R, geometry_name=geometry_name)

    # Formatting — tighter limits to reduce white space
    lim = 1.4
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_axis_off()
    fig.patch.set_visible(False)
    ax.view_init(elev=_VIEW_ELEV, azim=_VIEW_AZIM)
    ax.set_aspect("equal")
    ax.set_position([0.0, -0.15, 1.0, 1.15])

    # Title and subtitle
    stage_axis = np.asarray(stage.axis, dtype=float)
    title = f"{geometry_name} \u2014 {stage.name}" if geometry_name else stage.name
    sub = _subtitle(stage_axis, basis)
    fig.text(
        0.5,
        0.97,
        title,
        fontsize=12,
        fontweight="bold",
        ha="center",
        va="top",
    )
    fig.text(0.5, 0.91, sub, fontsize=10, ha="center", va="top", color="0.4")

    return fig


def draw_geometry_axes(
    geometry,
    *,
    figsize_per_stage: tuple[float, float] = (3.5, 2.8),
    max_cols: int = 4,
) -> matplotlib.figure.Figure:  # noqa: F821
    """Draw a composite figure with all stages of a geometry.

    Each subplot uses the same physical view convention as
    :func:`draw_stage_axis`: vertical up, longitudinal away, lateral left.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The geometry to draw.
    figsize_per_stage : tuple
        Size per subplot in inches.
    max_cols : int
        Maximum columns in the subplot grid.

    Returns
    -------
    matplotlib.figure.Figure
        The composite figure.
    """
    plt = _require_matplotlib()
    basis = geometry.basis
    R = _display_rotation(basis)

    all_stages = list(geometry.sample_stages) + list(geometry.detector_stages)
    n = len(all_stages)
    cols = min(n, max_cols)
    rows = (n + cols - 1) // cols

    fig = plt.figure(
        figsize=(figsize_per_stage[0] * cols, figsize_per_stage[1] * rows),
        frameon=False,
    )

    for i, stage in enumerate(all_stages):
        ax = fig.add_subplot(rows, cols, i + 1, projection="3d")

        _draw_stage_on_axes(
            ax,
            stage,
            basis,
            R,
            geometry_name=geometry.name,
            label_fontsize=7,
            arc_label_fontsize=7,
        )

        # Per-subplot formatting — tighter limits
        lim = 1.4
        ax.set_xlim(-lim, lim)
        ax.set_ylim(-lim, lim)
        ax.set_zlim(-lim, lim)
        ax.set_axis_off()
        fig.patch.set_visible(False)
        ax.view_init(elev=_VIEW_ELEV, azim=_VIEW_AZIM)
        ax.set_aspect("equal")

        stage_axis = np.asarray(stage.axis, dtype=float)
        phys = _physical_label(stage_axis, basis)
        hand = _handedness(stage_axis, basis)
        ax.set_title(f"{stage.name}\n({phys}, {hand})", fontsize=9)

    fig.suptitle(
        geometry.name,
        fontsize=14,
        fontweight="bold",
        y=1.02,
    )
    fig.subplots_adjust(hspace=0.3, wspace=0.1)
    return fig
