"""
Unit tests for ad_hoc_diffractometer.

Covers:
  - axes: parse_axis(), axis_label(), axis_from_physical()
  - rotation_matrix()
  - Stage class
  - AdHocDiffractometer class (construction, validation, ordering,
    set_angle, sample_rotation_matrix, detector_rotation_matrix)
  - psic(), fourc_v(), sixc() factories
  - lattice_vectors()
  - reciprocal_vectors()
  - b_matrix()

Convention for context parameters:
  Each parametrize set includes a 'context' entry that is either
    does_not_raise()          -- the call must succeed
    pytest.raises(Exc, ...)   -- the call must raise Exc
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

from ad_hoc_diffractometer import XHAT
from ad_hoc_diffractometer import YHAT
from ad_hoc_diffractometer import ZHAT
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import Stage
from ad_hoc_diffractometer import axis_from_physical
from ad_hoc_diffractometer import axis_label
from ad_hoc_diffractometer import b_matrix
from ad_hoc_diffractometer import fivec
from ad_hoc_diffractometer import fourc_h
from ad_hoc_diffractometer import fourc_v
from ad_hoc_diffractometer import kappa4c
from ad_hoc_diffractometer import kappa4c_h
from ad_hoc_diffractometer import kappa6c
from ad_hoc_diffractometer import kappa_axis
from ad_hoc_diffractometer import lattice_vectors
from ad_hoc_diffractometer import list_geometries
from ad_hoc_diffractometer import parse_axis
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import reciprocal_vectors
from ad_hoc_diffractometer import rotation_matrix
from ad_hoc_diffractometer import s2d2
from ad_hoc_diffractometer import sixc
from ad_hoc_diffractometer import zaxis

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def assert_matrix_close(A, B, atol=1e-10):
    """Assert two arrays are element-wise close."""
    np.testing.assert_allclose(A, B, atol=atol)


def Rx(deg):
    """Reference right-handed rotation about +x."""
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([[1, 0, 0], [0, c, -s], [0, s, c]])


def Ry(deg):
    """Reference right-handed rotation about +y."""
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, 0, s], [0, 1, 0], [-s, 0, c]])


def Rz(deg):
    """Reference right-handed rotation about +z."""
    t = np.deg2rad(deg)
    c, s = np.cos(t), np.sin(t)
    return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1]])


# ---------------------------------------------------------------------------
# Tests for axes: parse_axis(), axis_label(), axis_from_physical()
# ---------------------------------------------------------------------------

_STANDARD_BASIS = {
    "vertical": XHAT,
    "longitudinal": YHAT,
    "lateral": ZHAT,
}


@pytest.mark.parametrize(
    "label, expected, context",
    [
        # Signed Cartesian labels
        pytest.param("+x", +XHAT, does_not_raise(), id="plus-x"),
        pytest.param("-x", -XHAT, does_not_raise(), id="minus-x"),
        pytest.param("+y", +YHAT, does_not_raise(), id="plus-y"),
        pytest.param("-y", -YHAT, does_not_raise(), id="minus-y"),
        pytest.param("+z", +ZHAT, does_not_raise(), id="plus-z"),
        pytest.param("-z", -ZHAT, does_not_raise(), id="minus-z"),
        # Case insensitive
        pytest.param("+X", +XHAT, does_not_raise(), id="plus-X-uppercase"),
        pytest.param("-Z", -ZHAT, does_not_raise(), id="minus-Z-uppercase"),
        # Physical direction names (require basis)
        pytest.param(
            "vertical",
            +XHAT,
            does_not_raise(),
            id="physical-vertical-no-sign",
        ),
        pytest.param(
            "+vertical",
            +XHAT,
            does_not_raise(),
            id="physical-plus-vertical",
        ),
        pytest.param(
            "-vertical",
            -XHAT,
            does_not_raise(),
            id="physical-minus-vertical",
        ),
        pytest.param(
            "-lateral",
            -ZHAT,
            does_not_raise(),
            id="physical-minus-lateral",
        ),
        pytest.param(
            "+longitudinal",
            +YHAT,
            does_not_raise(),
            id="physical-plus-longitudinal",
        ),
        # Failures
        pytest.param(
            "vertical",
            None,
            pytest.raises(ValueError, match=re.escape("basis dict is required")),
            id="invalid-physical-name-no-basis",
        ),
        pytest.param(
            "+bogus",
            None,
            pytest.raises(ValueError, match=re.escape("basis dict is required")),
            id="invalid-unknown-label-no-basis",
        ),
    ],
)
def test_parse_axis(label, expected, context):
    # Cartesian labels (+x, -z, etc.) don't need a basis.
    # Physical direction names do.  Failure cases explicitly test no-basis behaviour.
    is_cartesian = label.strip().lower().lstrip("+-") in ("x", "y", "z")
    is_failure = expected is None
    if is_cartesian or is_failure:
        basis = None
    else:
        basis = _STANDARD_BASIS
    with context:
        result = parse_axis(label, basis=basis)
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "label, basis, expected, context",
    [
        pytest.param(
            "+x", None, +XHAT, does_not_raise(), id="cartesian-no-basis-needed"
        ),
        pytest.param(
            "vertical",
            _STANDARD_BASIS,
            +XHAT,
            does_not_raise(),
            id="physical-with-basis",
        ),
        pytest.param(
            "-lateral",
            _STANDARD_BASIS,
            -ZHAT,
            does_not_raise(),
            id="physical-minus-with-basis",
        ),
        pytest.param(
            "+bogus",
            _STANDARD_BASIS,
            None,
            pytest.raises(ValueError, match=re.escape("not found in basis dict")),
            id="invalid-direction-not-in-basis",
        ),
    ],
)
def test_parse_axis_with_basis(label, basis, expected, context):
    with context:
        result = parse_axis(label, basis=basis)
        np.testing.assert_array_equal(result, expected)


@pytest.mark.parametrize(
    "vector, expected_label, context",
    [
        pytest.param(+XHAT, "+x", does_not_raise(), id="plus-xhat"),
        pytest.param(-XHAT, "-x", does_not_raise(), id="minus-xhat"),
        pytest.param(+YHAT, "+y", does_not_raise(), id="plus-yhat"),
        pytest.param(-YHAT, "-y", does_not_raise(), id="minus-yhat"),
        pytest.param(+ZHAT, "+z", does_not_raise(), id="plus-zhat"),
        pytest.param(-ZHAT, "-z", does_not_raise(), id="minus-zhat"),
        # Non-standard vector falls back to numeric formatting
        pytest.param(
            np.array([0.5, 0.5, 0.0]),
            "[0.5, 0.5, 0]",
            does_not_raise(),
            id="non-standard-vector-numeric",
        ),
    ],
)
def test_axis_label(vector, expected_label, context):
    with context:
        assert axis_label(vector) == expected_label


@pytest.mark.parametrize(
    "direction, sign, expected, context",
    [
        pytest.param("vertical", "+", +XHAT, does_not_raise(), id="plus-vertical"),
        pytest.param("vertical", "-", -XHAT, does_not_raise(), id="minus-vertical"),
        pytest.param("lateral", "-", -ZHAT, does_not_raise(), id="minus-lateral"),
        pytest.param(
            "longitudinal", "+", +YHAT, does_not_raise(), id="plus-longitudinal"
        ),
        pytest.param(
            "bogus",
            "+",
            None,
            pytest.raises(ValueError, match=re.escape("not found in basis dict")),
            id="invalid-unknown-direction",
        ),
    ],
)
def test_axis_from_physical(direction, sign, expected, context):
    with context:
        result = axis_from_physical(direction, sign, _STANDARD_BASIS)
        np.testing.assert_array_equal(result, expected)


# ---------------------------------------------------------------------------
# Tests for rotation_matrix()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "axis, angle_deg, expected, context",
    [
        pytest.param(
            XHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-x-is-identity"
        ),
        pytest.param(
            YHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-y-is-identity"
        ),
        pytest.param(
            ZHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-z-is-identity"
        ),
        pytest.param(XHAT, 90.0, Rx(90), does_not_raise(), id="90deg-about-x"),
        pytest.param(YHAT, 90.0, Ry(90), does_not_raise(), id="90deg-about-y"),
        pytest.param(ZHAT, 90.0, Rz(90), does_not_raise(), id="90deg-about-z"),
        pytest.param(XHAT, 180.0, Rx(180), does_not_raise(), id="180deg-about-x"),
        pytest.param(YHAT, 180.0, Ry(180), does_not_raise(), id="180deg-about-y"),
        pytest.param(ZHAT, 180.0, Rz(180), does_not_raise(), id="180deg-about-z"),
        pytest.param(
            XHAT, -45.0, Rx(-45), does_not_raise(), id="negative-angle-about-x"
        ),
        pytest.param(
            ZHAT, -90.0, Rz(-90), does_not_raise(), id="negative-angle-about-z"
        ),
        pytest.param(
            -ZHAT,
            30.0,
            Rz(-30),
            does_not_raise(),
            id="negated-axis-equals-negated-angle",
        ),
        pytest.param(
            -XHAT,
            45.0,
            Rx(-45),
            does_not_raise(),
            id="negated-x-axis-equals-negated-angle",
        ),
        pytest.param(
            2 * XHAT, 90.0, Rx(90), does_not_raise(), id="non-unit-axis-is-normalised"
        ),
        pytest.param(YHAT, 360.0, np.eye(3), does_not_raise(), id="360deg-is-identity"),
    ],
)
def test_rotation_matrix(axis, angle_deg, expected, context):
    with context:
        R = rotation_matrix(axis, angle_deg)
        assert R.shape == (3, 3)
        assert_matrix_close(R, expected)


@pytest.mark.parametrize(
    "axis, angle_deg, context",
    [
        pytest.param(XHAT, 37.0, does_not_raise(), id="x-axis-37deg"),
        pytest.param(YHAT, -73.0, does_not_raise(), id="y-axis-neg73deg"),
        pytest.param(ZHAT, 120.0, does_not_raise(), id="z-axis-120deg"),
        pytest.param(-ZHAT, 55.0, does_not_raise(), id="neg-z-axis-55deg"),
        pytest.param(
            np.array([1, 1, 1]), 60.0, does_not_raise(), id="diagonal-axis-60deg"
        ),
    ],
)
def test_rotation_matrix_orthogonal(axis, angle_deg, context):
    with context:
        R = rotation_matrix(axis, angle_deg)
        assert_matrix_close(R @ R.T, np.eye(3))
        assert abs(np.linalg.det(R) - 1.0) < 1e-10


# ---------------------------------------------------------------------------
# Tests for Stage
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, axis, parent, role, angle, context",
    [
        pytest.param(
            "mu", XHAT, None, "sample", 0.0, does_not_raise(), id="sample-floor-stage"
        ),
        pytest.param(
            "eta",
            -ZHAT,
            "mu",
            "sample",
            45.0,
            does_not_raise(),
            id="sample-stacked-stage",
        ),
        pytest.param(
            "delta",
            -ZHAT,
            "nu",
            "detector",
            0.0,
            does_not_raise(),
            id="detector-stacked-stage",
        ),
        pytest.param(
            "phi",
            -ZHAT,
            None,
            "sample",
            0.0,
            does_not_raise(),
            id="sample-floor-stage-neg-z",
        ),
    ],
)
def test_stage_construction(name, axis, parent, role, angle, context):
    with context:
        s = Stage(name, axis, parent=parent, role=role, angle=angle)
        assert s.name == name
        np.testing.assert_array_equal(s.axis, np.asarray(axis, dtype=float))
        assert s.parent == parent
        assert s.role == role
        assert s.angle == angle


@pytest.mark.parametrize(
    "axis, angle_deg, expected_R, context",
    [
        pytest.param(
            XHAT, 0.0, np.eye(3), does_not_raise(), id="zero-angle-is-identity"
        ),
        pytest.param(XHAT, 90.0, Rx(90), does_not_raise(), id="90deg-about-x"),
        pytest.param(-ZHAT, 30.0, Rz(-30), does_not_raise(), id="left-handed-z-30deg"),
        pytest.param(YHAT, -45.0, Ry(-45), does_not_raise(), id="neg45deg-about-y"),
    ],
)
def test_stage_rotation_matrix(axis, angle_deg, expected_R, context):
    with context:
        s = Stage("test", axis, angle=angle_deg)
        assert_matrix_close(s.rotation_matrix(), expected_R)


# ---------------------------------------------------------------------------
# Tests for basis vector validation
# ---------------------------------------------------------------------------

_VALID_STAGES = [Stage("a", XHAT, parent=None, role="sample")]


@pytest.mark.parametrize(
    "basis, context",
    [
        pytest.param(
            {"lateral": XHAT, "longitudinal": YHAT, "vertical": ZHAT},
            does_not_raise(),
            id="valid-orthogonal-three-vectors",
        ),
        pytest.param(
            {"longitudinal": YHAT, "lateral": XHAT, "vertical": ZHAT},
            does_not_raise(),
            id="valid-orthogonal-different-dict-order",
        ),
        pytest.param(
            {"lateral": XHAT, "longitudinal": YHAT, "vertical": -ZHAT},
            does_not_raise(),
            id="valid-orthogonal-negated-vector",
        ),
        pytest.param(
            {"lateral": -XHAT, "longitudinal": YHAT, "vertical": ZHAT},
            does_not_raise(),
            id="valid-orthogonal-negated-first-vector",
        ),
        # --- failures ---
        pytest.param(
            {
                "lateral": np.array([0.0, 0.0, 1.0]),
                "longitudinal": np.array([0.0, 1.0, 0.0]),
                "vertical": np.array([0.0, 0.0, 1.0]),
            },
            pytest.raises(ValueError, match=re.escape("are not orthogonal")),
            id="invalid-two-vectors-identical",
        ),
        pytest.param(
            {
                "lateral": XHAT,
                "longitudinal": np.array([0.5, 0.5, 0.0]),
                "vertical": ZHAT,
            },
            pytest.raises(ValueError, match=re.escape("are not orthogonal")),
            id="invalid-non-orthogonal-pair",
        ),
        pytest.param(
            {
                "lateral": XHAT,
                "longitudinal": np.array([0.0, 1.0, 1.0]) / np.sqrt(2),
                "vertical": ZHAT,
            },
            pytest.raises(ValueError, match=re.escape("are not orthogonal")),
            id="invalid-third-pair-not-orthogonal",
        ),
        pytest.param(
            {"lateral": XHAT, "longitudinal": YHAT},
            pytest.raises(ValueError, match=re.escape("exactly 3 vectors")),
            id="invalid-only-two-basis-vectors",
        ),
        pytest.param(
            {"a": XHAT, "b": YHAT, "c": ZHAT, "d": XHAT},
            pytest.raises(ValueError, match=re.escape("exactly 3 vectors")),
            id="invalid-four-basis-vectors",
        ),
        pytest.param(
            {
                "lateral": np.array([0.0, 0.0, 0.0]),
                "longitudinal": YHAT,
                "vertical": ZHAT,
            },
            pytest.raises(ValueError, match=re.escape("non-zero")),
            id="invalid-zero-basis-vector",
        ),
        pytest.param(
            {"lateral": np.array([1.0, 0.0]), "longitudinal": YHAT, "vertical": ZHAT},
            pytest.raises(ValueError, match=re.escape("3-dimensional")),
            id="invalid-2d-basis-vector",
        ),
        pytest.param(
            {
                "lateral": np.array([1.0, 0.0, 0.0, 0.0]),
                "longitudinal": YHAT,
                "vertical": ZHAT,
            },
            pytest.raises(ValueError, match=re.escape("3-dimensional")),
            id="invalid-4d-basis-vector",
        ),
    ],
)
def test_basis_validation(basis, context):
    with context:
        AdHocDiffractometer("test", _VALID_STAGES, basis=basis)


# ---------------------------------------------------------------------------
# Tests for AdHocDiffractometer construction and validation
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "stages, context",
    [
        pytest.param(
            [
                Stage("a", XHAT, parent=None, role="sample"),
                Stage("b", YHAT, parent=None, role="detector"),
            ],
            does_not_raise(),
            id="valid-two-independent-floor-stages",
        ),
        pytest.param(
            [
                Stage("a", XHAT, parent=None, role="sample"),
                Stage("b", YHAT, parent="a", role="sample"),
            ],
            does_not_raise(),
            id="valid-simple-chain",
        ),
        pytest.param(
            [
                Stage("a", XHAT, parent=None, role="sample"),
                Stage("a", YHAT, parent=None, role="sample"),
            ],
            pytest.raises(ValueError, match=re.escape("Stage names must be unique")),
            id="invalid-duplicate-stage-names",
        ),
        pytest.param(
            [Stage("a", XHAT, parent="missing", role="sample")],
            pytest.raises(
                ValueError, match=re.escape("which is not in the stage list")
            ),
            id="invalid-parent-not-in-list",
        ),
        pytest.param(
            [
                Stage("a", XHAT, parent="b", role="sample"),
                Stage("b", YHAT, parent="a", role="sample"),
            ],
            pytest.raises(ValueError, match=re.escape("Cycle detected")),
            id="invalid-cycle-a-b-a",
        ),
    ],
)
def test_geometry_construction(stages, context):
    with context:
        AdHocDiffractometer("test", stages)


# ---------------------------------------------------------------------------
# Tests for stacking order
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "stages, role, expected_order, context",
    [
        pytest.param(
            [
                Stage("c", XHAT, parent="b", role="sample"),
                Stage("a", XHAT, parent=None, role="sample"),
                Stage("b", XHAT, parent="a", role="sample"),
            ],
            "sample",
            ["a", "b", "c"],
            does_not_raise(),
            id="three-stage-chain-submitted-out-of-order",
        ),
        pytest.param(
            [
                Stage("nu", XHAT, parent=None, role="detector"),
                Stage("delta", XHAT, parent="nu", role="detector"),
            ],
            "detector",
            ["nu", "delta"],
            does_not_raise(),
            id="two-detector-stages-floor-first",
        ),
        pytest.param(
            [
                Stage("mu", XHAT, parent=None, role="sample"),
                Stage("eta", -ZHAT, parent="mu", role="sample"),
                Stage("chi", YHAT, parent="eta", role="sample"),
                Stage("phi", -ZHAT, parent="chi", role="sample"),
            ],
            "sample",
            ["mu", "eta", "chi", "phi"],
            does_not_raise(),
            id="psic-sample-stack-order",
        ),
    ],
)
def test_stacking_order(stages, role, expected_order, context):
    with context:
        g = AdHocDiffractometer("test", stages)
        ordered = g.sample_stages if role == "sample" else g.detector_stages
        assert [s.name for s in ordered] == expected_order


# ---------------------------------------------------------------------------
# Tests for set_angle and stage()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "stage_name, angle, context",
    [
        pytest.param("mu", 30.0, does_not_raise(), id="set-mu-to-30"),
        pytest.param("eta", 0.0, does_not_raise(), id="set-eta-to-zero"),
        pytest.param("phi", -90.0, does_not_raise(), id="set-phi-to-neg90"),
        pytest.param(
            "bogus", 10.0, pytest.raises(KeyError), id="unknown-stage-raises-KeyError"
        ),
    ],
)
def test_set_angle(stage_name, angle, context):
    with context:
        g = psic()
        g.set_angle(stage_name, angle)
        assert g.stage(stage_name).angle == angle


# ---------------------------------------------------------------------------
# Tests for sample_rotation_matrix and detector_rotation_matrix
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "angles, expected_Z, context",
    [
        pytest.param(
            {"mu": 0, "eta": 0, "chi": 0, "phi": 0},
            np.eye(3),
            does_not_raise(),
            id="all-zero-is-identity",
        ),
        pytest.param(
            {"mu": 90, "eta": 0, "chi": 0, "phi": 0},
            Rx(90),
            does_not_raise(),
            id="mu-90-only-rotates-about-x",
        ),
        pytest.param(
            {"mu": 0, "eta": 30, "chi": 0, "phi": 0},
            Rz(-30),
            does_not_raise(),
            id="eta-30-only-left-handed-z",
        ),
        pytest.param(
            {"mu": 0, "eta": 0, "chi": 45, "phi": 0},
            Ry(45),
            does_not_raise(),
            id="chi-45-only-rotates-about-y",
        ),
        pytest.param(
            {"mu": 0, "eta": 0, "chi": 0, "phi": 60},
            Rz(-60),
            does_not_raise(),
            id="phi-60-only-left-handed-z",
        ),
    ],
)
def test_psic_sample_rotation(angles, expected_Z, context):
    with context:
        g = psic()
        for name, val in angles.items():
            g.set_angle(name, val)
        assert_matrix_close(g.sample_rotation_matrix(), expected_Z)


@pytest.mark.parametrize(
    "angles, expected_D, context",
    [
        pytest.param(
            {"nu": 0, "delta": 0},
            np.eye(3),
            does_not_raise(),
            id="all-zero-is-identity",
        ),
        pytest.param(
            {"nu": 90, "delta": 0},
            Rx(90),
            does_not_raise(),
            id="nu-90-only-rotates-about-x",
        ),
        pytest.param(
            {"nu": 0, "delta": 45},
            Rz(-45),
            does_not_raise(),
            id="delta-45-only-left-handed-z",
        ),
    ],
)
def test_psic_detector_rotation(angles, expected_D, context):
    with context:
        g = psic()
        for name, val in angles.items():
            g.set_angle(name, val)
        assert_matrix_close(g.detector_rotation_matrix(), expected_D)


# ---------------------------------------------------------------------------
# Tests for kappa_axis()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "alpha_deg, expected, context",
    [
        # At alpha=0: pure vertical (degenerate) — invalid
        pytest.param(
            0.0,
            None,
            pytest.raises(ValueError, match=re.escape("must be in (0, 90)")),
            id="invalid-alpha-zero",
        ),
        # At alpha=90: pure lateral — invalid
        pytest.param(
            90.0,
            None,
            pytest.raises(ValueError, match=re.escape("must be in (0, 90)")),
            id="invalid-alpha-90",
        ),
        # Negative alpha — invalid
        pytest.param(
            -10.0,
            None,
            pytest.raises(ValueError, match=re.escape("must be in (0, 90)")),
            id="invalid-alpha-negative",
        ),
        # Default 50 deg: cos(50)*xhat + sin(50)*zhat
        pytest.param(
            50.0,
            np.cos(np.deg2rad(50)) * XHAT + np.sin(np.deg2rad(50)) * ZHAT,
            does_not_raise(),
            id="alpha-50-default",
        ),
        # Near-vertical: small alpha
        pytest.param(
            1.0,
            np.cos(np.deg2rad(1)) * XHAT + np.sin(np.deg2rad(1)) * ZHAT,
            does_not_raise(),
            id="alpha-1-near-vertical",
        ),
        # Near-lateral: large alpha
        pytest.param(
            89.0,
            np.cos(np.deg2rad(89)) * XHAT + np.sin(np.deg2rad(89)) * ZHAT,
            does_not_raise(),
            id="alpha-89-near-lateral",
        ),
        # 45 deg: equal components
        pytest.param(
            45.0,
            np.array([np.sqrt(2) / 2, 0.0, np.sqrt(2) / 2]),
            does_not_raise(),
            id="alpha-45-equal-components",
        ),
    ],
)
def test_kappa_axis(alpha_deg, expected, context):
    with context:
        result = kappa_axis(alpha_deg)
        np.testing.assert_allclose(result, expected, atol=1e-12)


@pytest.mark.parametrize(
    "alpha_deg, context",
    [
        pytest.param(50.0, does_not_raise(), id="kappa-axis-is-unit-vector-50"),
        pytest.param(30.0, does_not_raise(), id="kappa-axis-is-unit-vector-30"),
        pytest.param(70.0, does_not_raise(), id="kappa-axis-is-unit-vector-70"),
    ],
)
def test_kappa_axis_is_unit_vector(alpha_deg, context):
    with context:
        ax = kappa_axis(alpha_deg)
        np.testing.assert_allclose(np.linalg.norm(ax), 1.0, atol=1e-12)


@pytest.mark.parametrize(
    "context",
    [pytest.param(does_not_raise(), id="kappa-axis-missing-basis-key")],
)
def test_kappa_axis_bad_basis(context):
    with pytest.raises(ValueError, match=re.escape("missing: 'lateral'")):
        kappa_axis(50.0, basis={"vertical": XHAT, "longitudinal": YHAT})


# ---------------------------------------------------------------------------
# Tests for list_geometries()
# ---------------------------------------------------------------------------


def test_list_geometries_returns_all_factories():
    geoms = list_geometries()
    expected = {
        "psic",
        "fourc_v",
        "fourc_h",
        "sixc",
        "kappa4c",
        "kappa4c_h",
        "kappa6c",
        "zaxis",
        "s2d2",
        "fivec",
    }
    assert set(geoms.keys()) == expected


def test_list_geometries_values_are_callable():
    for name, func in list_geometries().items():
        assert callable(func), f"{name} is not callable"


def test_list_geometries_instantiate_all():
    """Every registered factory must instantiate without error."""
    for name, func in list_geometries().items():
        g = func()
        assert g.name == name, f"Factory {name!r} returned name {g.name!r}"


def test_list_geometries_returns_copy():
    """Mutating the returned dict must not affect the registry."""
    geoms = list_geometries()
    geoms.clear()
    assert len(list_geometries()) > 0


# ---------------------------------------------------------------------------
# Tests for geometry factories
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_name, sample_names, detector_names, context",
    [
        pytest.param(
            psic,
            "psic",
            ["mu", "eta", "chi", "phi"],
            ["nu", "delta"],
            does_not_raise(),
            id="psic-stage-lists",
        ),
        pytest.param(
            fourc_v,
            "fourc_v",
            ["omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fourc_v-stage-lists",
        ),
        pytest.param(
            fourc_h,
            "fourc_h",
            ["omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fourc_h-stage-lists",
        ),
        pytest.param(
            sixc,
            "sixc",
            ["alpha", "omega", "chi", "phi"],
            ["delta", "gamma"],
            does_not_raise(),
            id="sixc-stage-lists",
        ),
        pytest.param(
            kappa4c,
            "kappa4c",
            ["komega", "kappa", "kphi"],
            ["two_theta"],
            does_not_raise(),
            id="kappa4c-stage-lists",
        ),
        pytest.param(
            kappa4c_h,
            "kappa4c_h",
            ["komega", "kappa", "kphi"],
            ["two_theta"],
            does_not_raise(),
            id="kappa4c_h-stage-lists",
        ),
        pytest.param(
            kappa6c,
            "kappa6c",
            ["mu", "komega", "kappa", "kphi"],
            ["nu", "delta"],
            does_not_raise(),
            id="kappa6c-stage-lists",
        ),
        pytest.param(
            zaxis,
            "zaxis",
            ["alpha", "Z"],
            ["delta", "gamma"],
            does_not_raise(),
            id="zaxis-stage-lists",
        ),
        pytest.param(
            s2d2,
            "s2d2",
            ["mu", "Z"],
            ["nu", "delta"],
            does_not_raise(),
            id="s2d2-stage-lists",
        ),
        pytest.param(
            fivec,
            "fivec",
            ["mu", "omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fivec-stage-lists",
        ),
    ],
)
def test_geometry_factories(
    factory, expected_name, sample_names, detector_names, context
):
    with context:
        g = factory()
        assert g.name == expected_name
        assert [s.name for s in g.sample_stages] == sample_names
        assert [s.name for s in g.detector_stages] == detector_names


@pytest.mark.parametrize(
    "factory, stage_name, expected_parent, context",
    [
        # psic
        pytest.param(psic, "mu", None, does_not_raise(), id="psic-mu-on-floor"),
        pytest.param(psic, "nu", None, does_not_raise(), id="psic-nu-on-floor"),
        pytest.param(psic, "eta", "mu", does_not_raise(), id="psic-eta-on-mu"),
        pytest.param(psic, "chi", "eta", does_not_raise(), id="psic-chi-on-eta"),
        pytest.param(psic, "phi", "chi", does_not_raise(), id="psic-phi-on-chi"),
        pytest.param(psic, "delta", "nu", does_not_raise(), id="psic-delta-on-nu"),
        # fourc_v: omega and two_theta both on floor (decoupled)
        pytest.param(
            fourc_v, "omega", None, does_not_raise(), id="fourc_v-omega-on-floor"
        ),
        pytest.param(
            fourc_v,
            "two_theta",
            None,
            does_not_raise(),
            id="fourc_v-two_theta-decoupled",
        ),
        pytest.param(
            fourc_v, "chi", "omega", does_not_raise(), id="fourc_v-chi-on-omega"
        ),
        pytest.param(fourc_v, "phi", "chi", does_not_raise(), id="fourc_v-phi-on-chi"),
        # fourc_h
        pytest.param(
            fourc_h,
            "two_theta",
            None,
            does_not_raise(),
            id="fourc_h-two_theta-decoupled",
        ),
        # sixc: shared alpha base
        pytest.param(sixc, "alpha", None, does_not_raise(), id="sixc-alpha-on-floor"),
        pytest.param(
            sixc, "omega", "alpha", does_not_raise(), id="sixc-omega-on-alpha"
        ),
        pytest.param(
            sixc, "delta", "alpha", does_not_raise(), id="sixc-delta-on-alpha"
        ),
        pytest.param(
            sixc, "gamma", "delta", does_not_raise(), id="sixc-gamma-on-delta"
        ),
        # kappa4c: komega and two_theta both on floor (decoupled)
        pytest.param(
            kappa4c, "komega", None, does_not_raise(), id="kappa4c-komega-on-floor"
        ),
        pytest.param(
            kappa4c,
            "two_theta",
            None,
            does_not_raise(),
            id="kappa4c-two_theta-decoupled",
        ),
        pytest.param(
            kappa4c, "kappa", "komega", does_not_raise(), id="kappa4c-kappa-on-komega"
        ),
        pytest.param(
            kappa4c, "kphi", "kappa", does_not_raise(), id="kappa4c-kphi-on-kappa"
        ),
        # kappa6c
        pytest.param(kappa6c, "mu", None, does_not_raise(), id="kappa6c-mu-on-floor"),
        pytest.param(kappa6c, "nu", None, does_not_raise(), id="kappa6c-nu-on-floor"),
        pytest.param(
            kappa6c, "komega", "mu", does_not_raise(), id="kappa6c-komega-on-mu"
        ),
        pytest.param(
            kappa6c, "kappa", "komega", does_not_raise(), id="kappa6c-kappa-on-komega"
        ),
        pytest.param(
            kappa6c, "kphi", "kappa", does_not_raise(), id="kappa6c-kphi-on-kappa"
        ),
        # zaxis: shared alpha base
        pytest.param(zaxis, "alpha", None, does_not_raise(), id="zaxis-alpha-on-floor"),
        pytest.param(zaxis, "Z", "alpha", does_not_raise(), id="zaxis-Z-on-alpha"),
        pytest.param(
            zaxis, "delta", "alpha", does_not_raise(), id="zaxis-delta-on-alpha"
        ),
        pytest.param(
            zaxis, "gamma", "delta", does_not_raise(), id="zaxis-gamma-on-delta"
        ),
        # s2d2: mu/nu both on floor, Z on mu, delta on nu
        pytest.param(s2d2, "mu", None, does_not_raise(), id="s2d2-mu-on-floor"),
        pytest.param(s2d2, "nu", None, does_not_raise(), id="s2d2-nu-on-floor"),
        pytest.param(s2d2, "Z", "mu", does_not_raise(), id="s2d2-Z-on-mu"),
        pytest.param(s2d2, "delta", "nu", does_not_raise(), id="s2d2-delta-on-nu"),
        # fivec: shared mu base
        pytest.param(fivec, "mu", None, does_not_raise(), id="fivec-mu-on-floor"),
        pytest.param(fivec, "omega", "mu", does_not_raise(), id="fivec-omega-on-mu"),
        pytest.param(
            fivec, "two_theta", "mu", does_not_raise(), id="fivec-two_theta-on-mu"
        ),
        pytest.param(fivec, "chi", "omega", does_not_raise(), id="fivec-chi-on-omega"),
        pytest.param(fivec, "phi", "chi", does_not_raise(), id="fivec-phi-on-chi"),
    ],
)
def test_geometry_parent_chain(factory, stage_name, expected_parent, context):
    with context:
        g = factory()
        assert g.stage(stage_name).parent == expected_parent


@pytest.mark.parametrize(
    "factory, stage_name, expected_axis, context",
    [
        pytest.param(
            psic, "mu", +XHAT, does_not_raise(), id="psic-mu-vertical-right-handed"
        ),
        pytest.param(
            psic, "eta", -ZHAT, does_not_raise(), id="psic-eta-lateral-left-handed"
        ),
        pytest.param(
            psic,
            "chi",
            +YHAT,
            does_not_raise(),
            id="psic-chi-longitudinal-right-handed",
        ),
        pytest.param(
            psic, "phi", -ZHAT, does_not_raise(), id="psic-phi-lateral-left-handed"
        ),
        pytest.param(
            psic, "nu", +XHAT, does_not_raise(), id="psic-nu-vertical-right-handed"
        ),
        pytest.param(
            psic, "delta", -ZHAT, does_not_raise(), id="psic-delta-lateral-left-handed"
        ),
        # kappa axis is a unit vector tilted 50 deg from vertical
        pytest.param(
            kappa4c,
            "kappa",
            np.cos(np.deg2rad(50)) * np.array([0, 0, 1])
            + np.sin(np.deg2rad(50)) * np.array([1, 0, 0]),
            does_not_raise(),
            id="kappa4c-kappa-axis-tilted",
        ),
    ],
)
def test_geometry_axes(factory, stage_name, expected_axis, context):
    with context:
        g = factory()
        np.testing.assert_allclose(g.stage(stage_name).axis, expected_axis, atol=1e-12)


# ---------------------------------------------------------------------------
# Tests for lattice_vectors()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(5.0, 5.0, 5.0, 90.0, 90.0, 90.0, does_not_raise(), id="cubic"),
        pytest.param(
            2.0, 3.0, 4.0, 90.0, 90.0, 90.0, does_not_raise(), id="orthorhombic"
        ),
        pytest.param(
            3.0, 3.0, 5.0, 90.0, 90.0, 120.0, does_not_raise(), id="hexagonal"
        ),
        pytest.param(
            4.0, 5.0, 6.0, 90.0, 100.0, 90.0, does_not_raise(), id="monoclinic"
        ),
    ],
)
def test_lattice_vectors_shapes(a, b, c, alpha, beta, gamma, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        assert a1.shape == (3,)
        assert a2.shape == (3,)
        assert a3.shape == (3,)


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(5.0, 5.0, 5.0, 90.0, 90.0, 90.0, does_not_raise(), id="cubic"),
        pytest.param(
            2.0, 3.0, 4.0, 90.0, 90.0, 90.0, does_not_raise(), id="orthorhombic"
        ),
        pytest.param(
            3.0, 3.0, 5.0, 90.0, 90.0, 120.0, does_not_raise(), id="hexagonal"
        ),
    ],
)
def test_lattice_vectors_magnitudes(a, b, c, alpha, beta, gamma, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        np.testing.assert_allclose(np.linalg.norm(a1), a, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(a2), b, atol=1e-12)
        np.testing.assert_allclose(np.linalg.norm(a3), c, atol=1e-12)


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(
            2.0,
            3.0,
            4.0,
            90.0,
            90.0,
            90.0,
            does_not_raise(),
            id="orthorhombic-aligned-with-axes",
        ),
    ],
)
def test_lattice_vectors_orthorhombic(a, b, c, alpha, beta, gamma, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        np.testing.assert_allclose(a1, [a, 0, 0], atol=1e-12)
        np.testing.assert_allclose(a2, [0, b, 0], atol=1e-12)
        np.testing.assert_allclose(a3, [0, 0, c], atol=1e-12)


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(
            3.0, 3.0, 5.0, 90.0, 90.0, 120.0, does_not_raise(), id="hexagonal-gamma-120"
        ),
        pytest.param(
            4.0, 5.0, 6.0, 90.0, 100.0, 90.0, does_not_raise(), id="monoclinic-beta-100"
        ),
    ],
)
def test_lattice_vectors_angles(a, b, c, alpha, beta, gamma, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        cos_alpha = np.dot(a2, a3) / (np.linalg.norm(a2) * np.linalg.norm(a3))
        np.testing.assert_allclose(np.rad2deg(np.arccos(cos_alpha)), alpha, atol=1e-10)
        cos_beta = np.dot(a1, a3) / (np.linalg.norm(a1) * np.linalg.norm(a3))
        np.testing.assert_allclose(np.rad2deg(np.arccos(cos_beta)), beta, atol=1e-10)
        cos_gamma = np.dot(a1, a2) / (np.linalg.norm(a1) * np.linalg.norm(a2))
        np.testing.assert_allclose(np.rad2deg(np.arccos(cos_gamma)), gamma, atol=1e-10)


# ---------------------------------------------------------------------------
# Tests for reciprocal_vectors()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(5.0, 5.0, 5.0, 90.0, 90.0, 90.0, does_not_raise(), id="cubic"),
        pytest.param(
            2.0, 3.0, 4.0, 90.0, 90.0, 90.0, does_not_raise(), id="orthorhombic"
        ),
        pytest.param(
            3.0, 3.0, 5.0, 90.0, 90.0, 120.0, does_not_raise(), id="hexagonal"
        ),
        pytest.param(
            4.0, 5.0, 6.0, 90.0, 100.0, 90.0, does_not_raise(), id="monoclinic"
        ),
    ],
)
def test_reciprocal_orthogonality(a, b, c, alpha, beta, gamma, context):
    """b_i . a_j = 2*pi * delta_ij for all i, j."""
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        twopi = 2 * np.pi
        np.testing.assert_allclose(np.dot(b1, a1), twopi, atol=1e-10)
        np.testing.assert_allclose(np.dot(b2, a2), twopi, atol=1e-10)
        np.testing.assert_allclose(np.dot(b3, a3), twopi, atol=1e-10)
        np.testing.assert_allclose(np.dot(b1, a2), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b1, a3), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b2, a1), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b2, a3), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b3, a1), 0.0, atol=1e-10)
        np.testing.assert_allclose(np.dot(b3, a2), 0.0, atol=1e-10)


@pytest.mark.parametrize(
    "a, b, c, context",
    [
        pytest.param(3.0, 3.0, 3.0, does_not_raise(), id="cubic-a3"),
        pytest.param(5.0, 5.0, 5.0, does_not_raise(), id="cubic-a5"),
    ],
)
def test_reciprocal_vectors_cubic(a, b, c, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, 90, 90, 90)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        twopi = 2 * np.pi
        np.testing.assert_allclose(b1, [twopi / a, 0, 0], atol=1e-12)
        np.testing.assert_allclose(b2, [0, twopi / b, 0], atol=1e-12)
        np.testing.assert_allclose(b3, [0, 0, twopi / c], atol=1e-12)


# ---------------------------------------------------------------------------
# Tests for b_matrix()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(5.0, 5.0, 5.0, 90.0, 90.0, 90.0, does_not_raise(), id="cubic"),
        pytest.param(
            2.0, 3.0, 4.0, 90.0, 90.0, 90.0, does_not_raise(), id="orthorhombic"
        ),
        pytest.param(
            3.0, 3.0, 5.0, 90.0, 90.0, 120.0, does_not_raise(), id="hexagonal"
        ),
    ],
)
def test_b_matrix_shape(a, b, c, alpha, beta, gamma, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        B = b_matrix(b1, b2, b3)
        assert B.shape == (3, 3)


@pytest.mark.parametrize(
    "a, b, c, context",
    [
        pytest.param(5.0, 5.0, 5.0, does_not_raise(), id="cubic-a5"),
        pytest.param(3.0, 3.0, 3.0, does_not_raise(), id="cubic-a3"),
    ],
)
def test_b_matrix_cubic(a, b, c, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, 90, 90, 90)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        B = b_matrix(b1, b2, b3)
        np.testing.assert_allclose(np.diag(B), [1 / a, 1 / b, 1 / c], atol=1e-12)
        off = B - np.diag(np.diag(B))
        np.testing.assert_allclose(off, np.zeros((3, 3)), atol=1e-12)


@pytest.mark.parametrize(
    "a, b, c, context",
    [
        pytest.param(2.0, 3.0, 4.0, does_not_raise(), id="orthorhombic-2-3-4"),
    ],
)
def test_b_matrix_orthorhombic(a, b, c, context):
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, 90, 90, 90)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        B = b_matrix(b1, b2, b3)
        np.testing.assert_allclose(np.diag(B), [1 / a, 1 / b, 1 / c], atol=1e-12)


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, h, context",
    [
        pytest.param(
            5.0,
            5.0,
            5.0,
            90,
            90,
            90,
            np.array([1, 0, 0]),
            does_not_raise(),
            id="cubic-h100",
        ),
        pytest.param(
            5.0,
            5.0,
            5.0,
            90,
            90,
            90,
            np.array([0, 1, 0]),
            does_not_raise(),
            id="cubic-h010",
        ),
        pytest.param(
            5.0,
            5.0,
            5.0,
            90,
            90,
            90,
            np.array([0, 0, 1]),
            does_not_raise(),
            id="cubic-h001",
        ),
        pytest.param(
            5.0,
            5.0,
            5.0,
            90,
            90,
            90,
            np.array([1, 2, 3]),
            does_not_raise(),
            id="cubic-h123",
        ),
        pytest.param(
            3.0,
            3.0,
            5.0,
            90,
            90,
            120,
            np.array([1, 0, 0]),
            does_not_raise(),
            id="hexagonal-h100",
        ),
        pytest.param(
            3.0,
            3.0,
            5.0,
            90,
            90,
            120,
            np.array([1, 1, 0]),
            does_not_raise(),
            id="hexagonal-h110",
        ),
        pytest.param(
            4.0,
            5.0,
            6.0,
            90,
            100,
            90,
            np.array([0, 0, 1]),
            does_not_raise(),
            id="monoclinic-h001",
        ),
    ],
)
def test_b_matrix_transforms_miller(a, b, c, alpha, beta, gamma, h, context):
    """
    2*pi * B.T @ h == h[0]*b1 + h[1]*b2 + h[2]*b3 for any Miller index h.

    The I16 convention is (b1, b2, b3) = 2*pi * B.T, meaning the j-th column
    of 2*pi*B.T is b_j.  Therefore:

        2*pi * B.T @ h = sum_j h[j] * b_j

    Note: B @ h (without transpose) gives the Cartesian coordinates of h in
    the crystal frame (hc = B @ h, Busing & Levy eq. 3), which is NOT the
    same as sum_j h[j]*b_j for non-orthogonal lattices.
    """
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        B = b_matrix(b1, b2, b3)
        expected = h[0] * b1 + h[1] * b2 + h[2] * b3
        np.testing.assert_allclose(2 * np.pi * B.T @ h, expected, atol=1e-10)


@pytest.mark.parametrize(
    "a, b, c, alpha, beta, gamma, context",
    [
        pytest.param(5.0, 5.0, 5.0, 90, 90, 90, does_not_raise(), id="cubic"),
        pytest.param(3.0, 3.0, 5.0, 90, 90, 120, does_not_raise(), id="hexagonal"),
        pytest.param(4.0, 5.0, 6.0, 90, 100, 90, does_not_raise(), id="monoclinic"),
    ],
)
def test_b_matrix_i16_convention(a, b, c, alpha, beta, gamma, context):
    """Verify (b1, b2, b3) = 2*pi * B.T (I16 diffractometer convention)."""
    with context:
        a1, a2, a3 = lattice_vectors(a, b, c, alpha, beta, gamma)
        b1, b2, b3 = reciprocal_vectors(a1, a2, a3)
        B = b_matrix(b1, b2, b3)
        rec_matrix = np.column_stack([b1, b2, b3])
        np.testing.assert_allclose(rec_matrix, 2 * np.pi * B.T, atol=1e-10)
