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

4. **Interactive Plotly figures** (:class:`StageAxisFigure`,
   :class:`GeometryAxisFigure`) — 3D interactive diagrams for a single
   stage or an entire geometry, arranged in a directed-graph layout.
   These are the preferred output for documentation; pre-generate them
   with ``tools/generate_geometry_drawings.py`` and embed the resulting
   HTML files as ``<iframe>`` elements in the Sphinx pages.

The DOT source is designed to be pasted into a Sphinx ``.. graphviz::``
directive in the geometry documentation pages.

Matplotlib is an optional dependency (listed in ``pyproject.toml``
``[project.optional-dependencies] doc``).  Functions that require it
raise :class:`ImportError` with a helpful message if it is not installed.

Plotly is an optional dependency.  :class:`StageAxisFigure` and
:class:`GeometryAxisFigure` raise :class:`ImportError` with a helpful
message if it is not installed.
"""

from __future__ import annotations

import numpy as np

from .axes import axis_label


def _physical_label(axis_vec: np.ndarray, basis: dict | None) -> str:
    """Return a human-readable axis label with physical direction if possible.

    For standard signed basis vectors, returns e.g. ``"-transverse"`` or
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
    "transverse": np.radians(130),
    "longitudinal": np.radians(90),
}

_BASIS_COLORS = {"vertical": "green", "longitudinal": "blue", "transverse": "red"}


def _require_matplotlib():
    """Import and return matplotlib.pyplot, or raise ImportError."""
    try:
        import matplotlib.pyplot as plt

        return plt  # pragma: no cover
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
        ``"transverse"`` keys.

    Returns
    -------
    np.ndarray, shape (3, 3)
        Rotation matrix R such that ``R @ v`` maps a vector from the
        Cartesian frame into display coordinates.
    """
    lat = np.asarray(basis["transverse"], dtype=float)
    lon = np.asarray(basis["longitudinal"], dtype=float)
    ver = np.asarray(basis["vertical"], dtype=float)
    return np.array([lat, lon, ver])


def _stage_draw_direction(stage_axis: np.ndarray, basis: dict) -> tuple:
    """Determine the drawing direction for a stage axis.

    Negative-signed axes (e.g. -transverse) are drawn in the positive basis
    direction for visibility; the sign is conveyed by the subtitle and
    arc labels.

    Returns
    -------
    tuple of (draw_direction, axis_type, is_negated)
        draw_direction : np.ndarray — unit vector for drawing
        axis_type : str or None — "vertical", "transverse", "longitudinal", or None
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


def _draw_stage_on_axes(  # pragma: no cover
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


def draw_stage_axis(  # pragma: no cover
    stage,
    basis: dict,
    *,
    geometry_name: str = "",
    figsize: tuple[float, float] = (3.5, 2.8),
) -> matplotlib.figure.Figure:  # noqa: F821
    """Draw a 3D diagram for a single stage showing its rotation axis.

    The view is oriented so that the vertical basis vector points up,
    the longitudinal vector points away (into the screen), and the
    transverse vector points to the left.  This convention is independent
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


def draw_geometry_axes(  # pragma: no cover
    geometry,
    *,
    figsize_per_stage: tuple[float, float] = (3.5, 2.8),
    max_cols: int = 4,
) -> matplotlib.figure.Figure:  # noqa: F821
    """Draw a composite figure with all stages of a geometry.

    Each subplot uses the same physical view convention as
    :func:`draw_stage_axis`: vertical up, longitudinal away, transverse left.

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


# ---------------------------------------------------------------------------
# Plotly interactive figures
# ---------------------------------------------------------------------------

_ROLE_COLORS = {"sample": "aliceblue", "detector": "seashell", "other": "honeydew"}

# GeometryAxisFigure layout constants
_CELL_W = 400  # px per column
_CELL_H = 400  # px per row
_TITLE_PX = 70  # px reserved above each scene for the HTML title annotation
_GAP = 0.02  # fractional gap between adjacent scenes (paper coords)
_EYE_SCALE = 2.0  # camera zoom-out relative to StageAxisFigure default


def _require_plotly():  # pragma: no cover
    """Import and return plotly.graph_objects, or raise ImportError."""
    try:
        import plotly.graph_objects as go

        return go
    except ImportError:
        raise ImportError(
            "plotly is required for StageAxisFigure and GeometryAxisFigure.  "
            "Install it with:  pip install plotly"
        ) from None


def _tree_layout(stages):  # pragma: no cover
    """Compute bottom-to-top tree layout positions for a stage dict.

    Returns
    -------
    depth : dict  stage_name -> int
        0 for roots (bottom), increasing toward leaves (top).
    x_pos : dict  stage_name -> float
        Leaf nodes get integer slots 0, 1, 2, …
        Parent nodes are centered over their children (may be half-integer).
    children : dict  stage_name -> list[str]
    roots : list[str]
    """
    children = {n: [] for n in stages}
    roots = []
    for n, s in stages.items():
        if s.parent:
            children[s.parent].append(n)
        else:
            roots.append(n)

    # depth = longest path from any root
    memo = {}

    def get_depth(name):
        if name not in memo:
            p = stages[name].parent
            memo[name] = 0 if p is None else 1 + get_depth(p)
        return memo[name]

    depth = {n: get_depth(n) for n in stages}

    # leaf order: depth-first traversal from each root
    def subtree_leaves(name):
        if not children[name]:
            return [name]
        result = []
        for c in children[name]:
            result.extend(subtree_leaves(c))
        return result

    leaves = []
    for r in roots:
        leaves.extend(subtree_leaves(r))

    # x position: leaf nodes get integer slots; parents center over children
    x_pos = {}

    def assign_x(name):
        kids = children[name]
        if not kids:
            x_pos[name] = float(leaves.index(name))
        else:
            for k in kids:
                assign_x(k)
            x_pos[name] = (
                min(x_pos[k] for k in kids) + max(x_pos[k] for k in kids)
            ) / 2.0

    for r in roots:
        assign_x(r)

    return depth, x_pos, children, roots


def _scene_domain(col, row, n_cols, n_rows, title_frac):  # pragma: no cover
    """Paper-coordinate domain for a scene at grid position (col, row).

    row=0 is the bottom (roots); row=n_rows-1 is the top (deepest leaves).
    title_frac is the fraction of each row height reserved for the title.
    Returns (x_domain, y_domain) as two [lo, hi] lists.
    """
    col_w = 1.0 / n_cols
    row_h = 1.0 / n_rows
    x0 = col * col_w + _GAP / 2
    x1 = (col + 1) * col_w - _GAP / 2
    y0 = row * row_h + _GAP / 2
    y1 = (row + 1) * row_h - _GAP / 2 - title_frac * row_h
    return [x0, x1], [y0, y1]


class StageAxisFigure:  # pragma: no cover
    """Interactive Plotly figure for a single stage rotation-axis diagram.

    Draws the geometry's basis vectors (grey, semi-transparent) and the
    stage's rotation axis (red) with a direction arc, in a 3-D interactive
    Plotly scene.

    Parameters
    ----------
    geometry_name : str
        Name of a preset geometry (e.g. ``'kappa6c'``, ``'zaxis'``).
    axis_name : str
        Name of a stage in that geometry (e.g. ``'komega'``, ``'Z'``).
    axis_labels : bool
        If ``True``, label the basis vectors (V, L, T) in the figure.
    **kwargs
        Forwarded to :class:`plotly.graph_objects.Figure`.

    Examples
    --------
    >>> from ad_hoc_diffractometer.drawing import StageAxisFigure
    >>> fig = StageAxisFigure("fourcv", "omega")
    >>> fig.write_html("omega.html")
    """

    def __init__(
        self,
        geometry_name: str,
        axis_name: str,
        axis_labels: bool = False,
        **kwargs,
    ):
        go = _require_plotly()
        import ad_hoc_diffractometer.presets as _presets

        self._go = go
        self._fig = go.Figure(**kwargs)

        factory = getattr(_presets, geometry_name)
        self.geometry = factory()
        self.stage = self.geometry._stages[axis_name]
        self.R = None  # set in _draw_basis_vectors

        self._draw_basis_vectors(axis_labels=axis_labels)
        self._draw_axis()
        self._set_titles()
        self._set_scene()

    # ------------------------------------------------------------------
    # Public Figure-compatible interface
    # ------------------------------------------------------------------

    @property
    def data(self):
        """Plotly traces (delegated to the internal Figure)."""
        return self._fig.data

    @property
    def layout(self):
        """Plotly layout (delegated to the internal Figure)."""
        return self._fig.layout

    def add_trace(self, trace):
        """Add a trace to the internal figure."""
        self._fig.add_trace(trace)

    def update_layout(self, *args, **kwargs):
        """Update layout on the internal figure."""
        self._fig.update_layout(*args, **kwargs)

    def write_html(self, path, **kwargs):
        """Write the figure to an HTML file.

        Parameters
        ----------
        path : str or pathlib.Path
        **kwargs
            Forwarded to :meth:`plotly.graph_objects.Figure.write_html`.
        """
        self._fig.write_html(path, **kwargs)

    def write_image(self, path, **kwargs):
        """Write the figure to a static image file (SVG, PNG, etc.).

        Parameters
        ----------
        path : str or pathlib.Path
        **kwargs
            Forwarded to :meth:`plotly.graph_objects.Figure.write_image`.
        """
        self._fig.write_image(path, **kwargs)

    def show(self, **kwargs):
        """Display the figure (delegates to the internal Figure)."""
        self._fig.show(**kwargs)

    # ------------------------------------------------------------------
    # Drawing helpers
    # ------------------------------------------------------------------

    def _draw_basis_vectors(self, axis_labels: bool = False) -> None:
        """Draw basis vectors as lines with cone tips."""
        go = self._go
        ver = np.asarray(self.geometry.basis["vertical"], dtype=float)
        lon = np.asarray(self.geometry.basis["longitudinal"], dtype=float)
        lat = np.asarray(self.geometry.basis["transverse"], dtype=float)

        # R maps basis vectors to display coordinates:
        #   lat  -> +x; lon -> -y; ver -> +z (up on screen)
        self.R = np.array([lat, -lon, ver])

        ver_d = self.R @ ver
        lon_d = self.R @ lon
        lat_d = self.R @ lat

        basis_length = 0.5

        for name, vec, color in [
            ("V", ver_d, "grey"),
            ("L", lon_d, "grey"),
            ("T", lat_d, "grey"),
        ]:
            v = vec * basis_length
            self._fig.add_trace(
                go.Scatter3d(
                    x=[0, v[0]],
                    y=[0, v[1]],
                    z=[0, v[2]],
                    mode="lines",
                    line=dict(color=color, width=6),
                    opacity=0.4,
                    showlegend=False,
                )
            )
            self._fig.add_trace(
                go.Cone(
                    x=[v[0]],
                    y=[v[1]],
                    z=[v[2]],
                    u=[vec[0] * 0.01],
                    v=[vec[1] * 0.01],
                    w=[vec[2] * 0.01],
                    colorscale=[[0, color], [1, color]],
                    showscale=False,
                    sizemode="absolute",
                    sizeref=0.1,
                    opacity=0.4,
                )
            )
            if axis_labels:
                tv = vec * (basis_length + 0.15)
                self._fig.add_trace(
                    go.Scatter3d(
                        x=[tv[0]],
                        y=[tv[1]],
                        z=[tv[2]],
                        mode="text",
                        text=[name],
                        textfont=dict(size=14, color=color),
                        showlegend=False,
                    )
                )

    def _draw_axis(self) -> None:
        """Draw the axis vector and rotation arc."""
        axis_norm = self._draw_axis_vector()

        stage_axis = np.asarray(self.stage.axis, dtype=float)
        stage_norm = stage_axis / np.linalg.norm(stage_axis)
        lat_norm = np.asarray(self.geometry.basis["transverse"], dtype=float)
        lat_norm = lat_norm / np.linalg.norm(lat_norm)
        is_transverse = np.isclose(abs(np.dot(stage_norm, lat_norm)), 1.0, atol=1e-6)
        direction = -1 if is_transverse else 1

        self._draw_axis_arc(axis_norm, direction=direction)

    def _draw_axis_vector(self) -> np.ndarray:
        """Draw the rotation axis as a red line; return positive unit vector."""
        go = self._go
        stage_axis = np.asarray(self.stage.axis, dtype=float)
        axis_norm = stage_axis / np.linalg.norm(stage_axis)

        for component in axis_norm:
            if abs(component) > 1e-9:
                if component < 0:
                    axis_norm = -axis_norm
                break

        axis_d = self.R @ axis_norm
        v = axis_d * 0.4
        self._fig.add_trace(
            go.Scatter3d(
                x=[0, v[0]],
                y=[0, v[1]],
                z=[0, v[2]],
                mode="lines",
                line=dict(color="red", width=20),
                opacity=0.6,
                showlegend=False,
            )
        )
        return axis_norm

    def _draw_axis_arc(self, axis_norm: np.ndarray, direction: int = 1) -> None:
        """Draw the rotation arc and arrowhead."""
        go = self._go
        if direction == 0:
            return

        axis_d = self.R @ axis_norm

        eye = np.array([0.92, 1.18, -0.10])
        view_dir = -eye / np.linalg.norm(eye)
        cam_up = np.array([0.0, 0.0, 1.0])
        screen_right = np.cross(cam_up, view_dir)
        screen_right /= np.linalg.norm(screen_right)
        screen_up = np.cross(view_dir, screen_right)

        perp1_d = screen_right - np.dot(screen_right, axis_d) * axis_d
        if np.linalg.norm(perp1_d) < 1e-6:
            perp1_d = screen_up - np.dot(screen_up, axis_d) * axis_d
        perp1_d /= np.linalg.norm(perp1_d)
        perp2_d = np.cross(axis_d, perp1_d)

        r0 = np.array([np.dot(perp1_d, screen_right), np.dot(perp1_d, screen_up)])
        r1 = np.array([np.dot(perp2_d, screen_right), np.dot(perp2_d, screen_up)])
        cross_2d = r0[0] * r1[1] - r0[1] * r1[0]
        if cross_2d > 0:
            perp2_d = -perp2_d

        if direction < 0:
            perp2_d = -perp2_d

        arc_radius = 0.1
        arc_center_d = axis_d * 0.45
        arc_angles = np.linspace(-np.pi / 3, np.pi / 3, 60)
        arc_pts_d = np.array(
            [
                arc_center_d + arc_radius * (np.cos(a) * perp1_d + np.sin(a) * perp2_d)
                for a in arc_angles
            ]
        )
        arc_scale = 0.6

        self._fig.add_trace(
            go.Scatter3d(
                x=arc_pts_d[:, 0] * arc_scale,
                y=arc_pts_d[:, 1] * arc_scale,
                z=arc_pts_d[:, 2] * arc_scale,
                mode="lines",
                line=dict(color="black", width=3),
                opacity=0.6,
                showlegend=False,
            )
        )

        tangent = arc_pts_d[-1] - arc_pts_d[-3]
        tangent = tangent / np.linalg.norm(tangent)
        tip = arc_pts_d[-1] * arc_scale
        self._fig.add_trace(
            go.Cone(
                x=[tip[0]],
                y=[tip[1]],
                z=[tip[2]],
                u=[tangent[0] * 0.1],
                v=[tangent[1] * 0.1],
                w=[tangent[2] * 0.1],
                colorscale=[[0, "black"], [1, "black"]],
                showscale=False,
                sizemode="absolute",
                sizeref=0.06,
                opacity=0.6,
            )
        )

    def _set_titles(self) -> None:
        """Add title and subtitle to the figure."""
        stage_axis = np.asarray(self.stage.axis, dtype=float)
        phys = _physical_label(stage_axis, self.geometry.basis)
        cartesian = axis_label(stage_axis)
        if phys != cartesian:
            axis_desc = f"{phys}  ({cartesian})"
        else:
            axis_desc = (
                f"[{stage_axis[0]:.3f}, {stage_axis[1]:.3f}, {stage_axis[2]:.3f}]"
            )
        role = self.stage.role
        subtitle = f"{axis_desc} \u2014 {role}"
        title_text = (
            f"<b>{self.geometry.name} \u2014 {self.stage.name}</b><br>"
            f'<span style="color:gray;font-size:12px">{subtitle}</span>'
        )
        self._fig.update_layout(
            title=dict(
                text=title_text,
                x=0.5,
                xanchor="center",
                font=dict(size=14),
            ),
        )

    def _set_scene(self) -> None:
        """Configure the 3-D scene camera and axes."""
        self._fig.update_layout(
            scene=dict(
                xaxis=dict(visible=False),
                yaxis=dict(visible=False),
                zaxis=dict(visible=False),
                aspectmode="cube",
                camera=dict(
                    center=dict(x=0, y=0, z=0),
                    eye=dict(x=0.92, y=1.18, z=-0.10),
                    up=dict(x=0, y=0, z=1),
                ),
            ),
            margin=dict(l=0, r=0, t=50, b=0),
            width=400,
            height=450,
        )


class GeometryAxisFigure:  # pragma: no cover
    """Interactive Plotly figure showing all stages of a geometry.

    Stages are arranged as a directed graph following parent relationships.
    The graph flows bottom-to-top: root stages (``parent=None``) are at the
    bottom; each child is one row above its parent.  Leaf nodes are assigned
    integer column slots; parents are centered over their children.

    Role is indicated by scene background color:

    - sample   → aliceblue
    - detector → seashell
    - other    → honeydew

    Parent-child relationships are drawn as arrows on the figure.

    Parameters
    ----------
    geometry_name : str
        Name of a preset geometry (e.g. ``'zaxis'``, ``'kappa6c'``).
    axis_labels : bool
        If ``True``, label the basis vectors on each sub-plot.
    **kwargs
        Forwarded to :class:`plotly.graph_objects.Figure`.

    Examples
    --------
    >>> from ad_hoc_diffractometer.drawing import GeometryAxisFigure
    >>> fig = GeometryAxisFigure("fourcv")
    >>> fig.write_html("fourcv.html")
    """

    def __init__(
        self,
        geometry_name: str,
        axis_labels: bool = True,
        **kwargs,
    ):
        go = _require_plotly()
        import ad_hoc_diffractometer.presets as _presets

        self._go = go
        self._fig = go.Figure(**kwargs)

        factory = getattr(_presets, geometry_name)
        geometry = factory()
        stages = geometry._stages

        depth, x_pos, children, roots = _tree_layout(stages)
        max_depth = max(depth.values())
        n_rows = max_depth + 1
        leaf_xs = [x_pos[n] for n in stages if not children[n]]
        n_cols = int(max(leaf_xs)) + 1
        title_frac = _TITLE_PX / _CELL_H

        self._fig.update_layout(
            width=_CELL_W * n_cols,
            height=_CELL_H * n_rows,
            margin=dict(l=10, r=10, t=10, b=10),
            showlegend=False,
            xaxis=dict(visible=False),
            yaxis=dict(visible=False),
        )

        annotations = []
        scene_domains = {}
        scene_idx = 1

        for name, stage in stages.items():
            col = x_pos[name]
            row = depth[name]
            x_domain, y_domain = _scene_domain(col, row, n_cols, n_rows, title_frac)
            scene_domains[name] = (x_domain, y_domain)
            scene_key = "scene" if scene_idx == 1 else f"scene{scene_idx}"
            scene_idx += 1

            sf = StageAxisFigure(geometry_name, name, axis_labels=axis_labels)

            for trace in sf.data:
                trace.update(scene=scene_key)
                self._fig.add_trace(trace)

            role = getattr(stage, "role", "other")
            scene_layout = sf.layout.scene.to_plotly_json()
            scene_layout["domain"] = {"x": x_domain, "y": y_domain}
            scene_layout["bgcolor"] = _ROLE_COLORS.get(role, _ROLE_COLORS["other"])
            eye = scene_layout.get("camera", {}).get("eye", {})
            scene_layout["camera"]["eye"] = {k: v * _EYE_SCALE for k, v in eye.items()}
            self._fig.update_layout({scene_key: scene_layout})

            x_mid = (x_domain[0] + x_domain[1]) / 2
            y_top = y_domain[1]
            annotations.append(
                dict(
                    text=sf.layout.title.text,
                    x=x_mid,
                    y=y_top,
                    xref="paper",
                    yref="paper",
                    xanchor="center",
                    yanchor="bottom",
                    showarrow=False,
                    font=dict(size=12),
                )
            )

        # Parent-child edges
        fig_w = _CELL_W * n_cols
        fig_h = _CELL_H * n_rows
        margin_l = margin_t = 10
        plot_w = fig_w - 2 * margin_l
        plot_h = fig_h - 2 * margin_t

        def paper_to_px(xp, yp):
            return margin_l + xp * plot_w, margin_t + (1.0 - yp) * plot_h

        for name, stage in stages.items():
            if stage.parent is None:
                continue
            parent_name = stage.parent
            px, py = scene_domains[parent_name]
            cx, cy = scene_domains[name]

            x_start = (px[0] + px[1]) / 2
            y_start = py[1]
            x_end = (cx[0] + cx[1]) / 2
            y_end = cy[0]

            self._fig.add_shape(
                type="line",
                x0=x_start,
                y0=y_start,
                x1=x_end,
                y1=y_end,
                xref="paper",
                yref="paper",
                line=dict(color="dimgray", width=1.5),
            )

            sx_px, sy_px = paper_to_px(x_start, y_start)
            ex_px, ey_px = paper_to_px(x_end, y_end)
            dx = sx_px - ex_px
            dy = sy_px - ey_px
            dist = np.hypot(dx, dy)
            shaft_px = 12.0
            ax_px = dx / dist * shaft_px if dist > 0 else 0.0
            ay_px = dy / dist * shaft_px if dist > 0 else shaft_px
            annotations.append(
                dict(
                    x=x_end,
                    y=y_end,
                    xref="paper",
                    yref="paper",
                    ax=ax_px,
                    ay=ay_px,
                    axref="pixel",
                    ayref="pixel",
                    showarrow=True,
                    arrowhead=2,
                    arrowsize=1.2,
                    arrowwidth=1.5,
                    arrowcolor="dimgray",
                    text="",
                )
            )

        self._fig.update_layout(annotations=annotations)

    # ------------------------------------------------------------------
    # Public Figure-compatible interface
    # ------------------------------------------------------------------

    @property
    def data(self):
        """Plotly traces (delegated to the internal Figure)."""
        return self._fig.data

    @property
    def layout(self):
        """Plotly layout (delegated to the internal Figure)."""
        return self._fig.layout

    def write_html(self, path, **kwargs):
        """Write the figure to an HTML file.

        Parameters
        ----------
        path : str or pathlib.Path
        **kwargs
            Forwarded to :meth:`plotly.graph_objects.Figure.write_html`.
        """
        self._fig.write_html(path, **kwargs)

    def write_image(self, path, **kwargs):
        """Write the figure to a static image file (SVG, PNG, etc.).

        Parameters
        ----------
        path : str or pathlib.Path
        **kwargs
            Forwarded to :meth:`plotly.graph_objects.Figure.write_image`.
        """
        self._fig.write_image(path, **kwargs)

    def show(self, **kwargs):
        """Display the figure (delegates to the internal Figure)."""
        self._fig.show(**kwargs)
