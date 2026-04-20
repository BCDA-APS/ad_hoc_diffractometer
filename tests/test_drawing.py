# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.drawing.

Covers:
  - geometry_dot() produces valid DOT source for all preset geometries
  - DOT source contains expected node names, edges, and role colours
  - _physical_label() and _handedness() helper functions
  - draw_stage_axis() returns a matplotlib Figure (if matplotlib available)
  - draw_geometry_axes() returns a matplotlib Figure (if matplotlib available)
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

from ad_hoc_diffractometer.drawing import _handedness
from ad_hoc_diffractometer.drawing import _physical_label
from ad_hoc_diffractometer.drawing import geometry_dot
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.factories import BASIS_YOU
from ad_hoc_diffractometer.presets import fivec
from ad_hoc_diffractometer.presets import fourch
from ad_hoc_diffractometer.presets import fourcv
from ad_hoc_diffractometer.presets import kappa4ch
from ad_hoc_diffractometer.presets import kappa4cv
from ad_hoc_diffractometer.presets import kappa6c
from ad_hoc_diffractometer.presets import psic
from ad_hoc_diffractometer.presets import s2d2
from ad_hoc_diffractometer.presets import sixc
from ad_hoc_diffractometer.presets import zaxis

# ---------------------------------------------------------------------------
# _physical_label
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis_vec, basis, expected, context",
    [
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            BASIS_BL,
            "+lateral",
            does_not_raise(),
            id="positive-lateral-BL",
        ),
        pytest.param(
            np.array([-1.0, 0.0, 0.0]),
            BASIS_BL,
            "-lateral",
            does_not_raise(),
            id="negative-lateral-BL",
        ),
        pytest.param(
            np.array([0.0, 0.0, 1.0]),
            BASIS_BL,
            "+vertical",
            does_not_raise(),
            id="positive-vertical-BL",
        ),
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            BASIS_YOU,
            "+vertical",
            does_not_raise(),
            id="positive-vertical-YOU",
        ),
        pytest.param(
            np.array([0.5, 0.0, 0.5]),
            BASIS_BL,
            "[0.5, 0, 0.5]",  # falls back to axis_label for non-standard vectors
            does_not_raise(),
            id="non-standard-axis",
        ),
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            None,
            "+x",
            does_not_raise(),
            id="no-basis",
        ),
    ],
)
def test_physical_label(axis_vec, basis, expected, context):
    """_physical_label returns expected direction labels."""
    with context:
        result = _physical_label(axis_vec, basis)
        assert result == expected


# ---------------------------------------------------------------------------
# _handedness
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis_vec, basis, expected, context",
    [
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            BASIS_BL,
            "RH",
            does_not_raise(),
            id="positive-RH",
        ),
        pytest.param(
            np.array([-1.0, 0.0, 0.0]),
            BASIS_BL,
            "LH",
            does_not_raise(),
            id="negative-LH",
        ),
        pytest.param(
            np.array([0.5, 0.0, 0.5]),
            BASIS_BL,
            "RH",
            does_not_raise(),
            id="tilted-defaults-RH",
        ),
        pytest.param(
            np.array([1.0, 0.0, 0.0]),
            None,
            "",
            does_not_raise(),
            id="no-basis-empty",
        ),
    ],
)
def test_handedness(axis_vec, basis, expected, context):
    """_handedness returns RH, LH, or empty string."""
    with context:
        result = _handedness(axis_vec, basis)
        assert result == expected


# ---------------------------------------------------------------------------
# geometry_dot — all presets
# ---------------------------------------------------------------------------


ALL_PRESETS = [
    pytest.param(fivec, id="fivec"),
    pytest.param(fourch, id="fourch"),
    pytest.param(fourcv, id="fourcv"),
    pytest.param(kappa4ch, id="kappa4ch"),
    pytest.param(kappa4cv, id="kappa4cv"),
    pytest.param(kappa6c, id="kappa6c"),
    pytest.param(psic, id="psic"),
    pytest.param(s2d2, id="s2d2"),
    pytest.param(sixc, id="sixc"),
    pytest.param(zaxis, id="zaxis"),
]


@pytest.mark.parametrize("factory", ALL_PRESETS)
def test_geometry_dot_produces_valid_dot(factory):
    """geometry_dot returns a string starting with 'digraph' and ending with '}'."""
    g = factory()
    dot = geometry_dot(g)
    assert dot.startswith(f"digraph {g.name}")
    assert dot.rstrip().endswith("}")


@pytest.mark.parametrize("factory", ALL_PRESETS)
def test_geometry_dot_contains_all_stage_names(factory):
    """Every stage name appears as a node in the DOT source."""
    g = factory()
    dot = geometry_dot(g)
    all_stages = list(g.sample_stages) + list(g.detector_stages)
    for stage in all_stages:
        assert stage.name in dot, f"Stage {stage.name!r} missing from DOT"


@pytest.mark.parametrize("factory", ALL_PRESETS)
def test_geometry_dot_contains_edges(factory):
    """Child-to-parent edges appear in the DOT source (BT layout)."""
    g = factory()
    dot = geometry_dot(g)
    all_stages = list(g.sample_stages) + list(g.detector_stages)
    for stage in all_stages:
        if stage.parent is not None:
            edge = f"{stage.name} -> {stage.parent}"
            assert edge in dot, f"Edge {edge!r} missing from DOT"


@pytest.mark.parametrize("factory", ALL_PRESETS)
def test_geometry_dot_contains_role_colours(factory):
    """Sample nodes are blue, detector nodes are red."""
    g = factory()
    dot = geometry_dot(g)
    assert "#a8d8ea" in dot  # sample colour
    assert "#f8a5a5" in dot  # detector colour


def test_geometry_dot_sixc_shared_base():
    """sixc has alpha as a shared base for both sample and detector stacks."""
    g = sixc()
    dot = geometry_dot(g)
    assert "omega -> alpha" in dot
    assert "delta -> alpha" in dot


def test_geometry_dot_fivec_shared_base():
    """fivec has mu as a shared base for both sample and detector stacks."""
    g = fivec()
    dot = geometry_dot(g)
    assert "omega -> mu" in dot
    assert "ttheta -> mu" in dot


# ---------------------------------------------------------------------------
# Matplotlib tests (require matplotlib)
# ---------------------------------------------------------------------------


def test_require_matplotlib_missing():
    """_require_matplotlib raises ImportError when matplotlib is not available."""
    import builtins
    from unittest.mock import patch

    from ad_hoc_diffractometer.drawing import _require_matplotlib

    real_import = builtins.__import__

    def mock_import(name, *args, **kwargs):
        if name == "matplotlib.pyplot" or name == "matplotlib":
            raise ImportError("mocked")
        return real_import(name, *args, **kwargs)

    with patch("builtins.__import__", side_effect=mock_import):
        with pytest.raises(ImportError, match=re.escape("matplotlib is required")):
            _require_matplotlib()


def test_geometry_dot_no_basis():
    """geometry_dot works when geometry.basis is None (no handedness label)."""
    from ad_hoc_diffractometer.drawing import geometry_dot

    g = fourcv()
    g.basis = None
    dot = geometry_dot(g)
    assert "omega" in dot
    # No handedness label when basis is None
    assert "RH" not in dot
    assert "LH" not in dot


def test_display_rotation_bl1967():
    """_display_rotation maps BL1967 basis to lat=+x, lon=+y, ver=+z."""
    from ad_hoc_diffractometer.drawing import _display_rotation

    R = _display_rotation(BASIS_BL)
    # BL1967: lat=+x, lon=+y, ver=+z → R should be identity
    np.testing.assert_allclose(R, np.eye(3), atol=1e-12)


def test_display_rotation_you():
    """_display_rotation maps YOU basis so ver/lon/lat display correctly."""
    from ad_hoc_diffractometer.drawing import _display_rotation

    R = _display_rotation(BASIS_YOU)
    # YOU: ver=+x, lon=+y, lat=+z
    # Display: lat→row0, lon→row1, ver→row2
    # So R @ [1,0,0] (ver in YOU) should give [0,0,1] (ver in display)
    ver_disp = R @ np.array([1.0, 0, 0])
    np.testing.assert_allclose(ver_disp, [0, 0, 1], atol=1e-12)
    lat_disp = R @ np.array([0, 0, 1.0])
    np.testing.assert_allclose(lat_disp, [1, 0, 0], atol=1e-12)


def test_stage_draw_direction_positive():
    """_stage_draw_direction returns original direction for positive basis axis."""
    from ad_hoc_diffractometer.drawing import _stage_draw_direction

    axis = np.array([0.0, 1.0, 0.0])  # +lon in BL1967
    draw_dir, axis_type, is_neg = _stage_draw_direction(axis, BASIS_BL)
    np.testing.assert_allclose(draw_dir, [0, 1, 0], atol=1e-12)
    assert axis_type == "longitudinal"
    assert is_neg is False


def test_stage_draw_direction_negative():
    """_stage_draw_direction flips negative axis to positive for drawing."""
    from ad_hoc_diffractometer.drawing import _stage_draw_direction

    axis = np.array([-1.0, 0.0, 0.0])  # -lat in BL1967
    draw_dir, axis_type, is_neg = _stage_draw_direction(axis, BASIS_BL)
    np.testing.assert_allclose(draw_dir, [1, 0, 0], atol=1e-12)
    assert axis_type == "lateral"
    assert is_neg is True


def test_stage_draw_direction_kappa():
    """_stage_draw_direction returns as-is for non-standard (kappa) axis."""
    from ad_hoc_diffractometer.drawing import _stage_draw_direction

    axis = np.array([0.766, 0.0, 0.6428])
    draw_dir, axis_type, is_neg = _stage_draw_direction(axis, BASIS_BL)
    assert axis_type is None
    assert is_neg is False


def test_subtitle_standard_axis():
    """_subtitle shows physical and cartesian labels for standard axes."""
    from ad_hoc_diffractometer.drawing import _subtitle

    sub = _subtitle(np.array([0.0, 1.0, 0.0]), BASIS_BL)
    assert "+longitudinal" in sub
    assert "+y" in sub


def test_subtitle_kappa_axis():
    """_subtitle shows direction cosines for non-standard axes."""
    from ad_hoc_diffractometer.drawing import _subtitle

    sub = _subtitle(np.array([0.766, 0.0, 0.6428]), BASIS_BL)
    assert "0.766" in sub
    assert "0.643" in sub


def test_draw_stage_axis_returns_figure():
    """draw_stage_axis returns a matplotlib Figure."""
    import matplotlib

    matplotlib.use("Agg")
    from ad_hoc_diffractometer.drawing import draw_stage_axis

    g = fourcv()
    fig = draw_stage_axis(g.sample_stages[0], g.basis, geometry_name="fourcv")
    assert type(fig).__name__ == "Figure"
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_draw_stage_axis_kappa():
    """draw_stage_axis works for a tilted kappa axis (no arc offset)."""
    import matplotlib

    matplotlib.use("Agg")
    from ad_hoc_diffractometer.drawing import draw_stage_axis

    g = kappa4cv()
    # kappa stage is index 1 (komega, kappa, kphi, ttheta)
    fig = draw_stage_axis(g.sample_stages[1], g.basis, geometry_name="kappa4cv")
    assert type(fig).__name__ == "Figure"
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_draw_stage_axis_you_basis():
    """draw_stage_axis works with YOU basis (psic)."""
    import matplotlib

    matplotlib.use("Agg")
    from ad_hoc_diffractometer.drawing import draw_stage_axis

    g = psic()
    fig = draw_stage_axis(g.sample_stages[0], g.basis, geometry_name="psic")
    assert type(fig).__name__ == "Figure"
    import matplotlib.pyplot as plt

    plt.close(fig)


def test_draw_geometry_axes_returns_figure():
    """draw_geometry_axes returns a matplotlib Figure with correct subplot count."""
    import matplotlib

    matplotlib.use("Agg")
    from ad_hoc_diffractometer.drawing import draw_geometry_axes

    g = psic()
    fig = draw_geometry_axes(g)
    assert type(fig).__name__ == "Figure"
    # psic has 6 stages → 6 subplots
    axes_3d = [a for a in fig.get_axes() if hasattr(a, "view_init")]
    assert len(axes_3d) == 6
    import matplotlib.pyplot as plt

    plt.close(fig)
