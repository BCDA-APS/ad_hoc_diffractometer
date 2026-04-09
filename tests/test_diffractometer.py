"""
Unit tests for ad_hoc_diffractometer.

Covers:
  - rotation_matrix()
  - Stage class
  - AdHocDiffractometer class (construction, validation, ordering,
    set_angle, sample_rotation_matrix, detector_rotation_matrix)
  - geometry_psic(), geometry_fourc(), geometry_sixc() factories
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
from ad_hoc_diffractometer import b_matrix
from ad_hoc_diffractometer import geometry_fourc
from ad_hoc_diffractometer import geometry_psic
from ad_hoc_diffractometer import geometry_sixc
from ad_hoc_diffractometer import lattice_vectors
from ad_hoc_diffractometer import reciprocal_vectors
from ad_hoc_diffractometer import rotation_matrix

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
        g = geometry_psic()
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
        g = geometry_psic()
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
        g = geometry_psic()
        for name, val in angles.items():
            g.set_angle(name, val)
        assert_matrix_close(g.detector_rotation_matrix(), expected_D)


# ---------------------------------------------------------------------------
# Tests for geometry factories
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, expected_name, sample_names, detector_names, context",
    [
        pytest.param(
            geometry_psic,
            "psic",
            ["mu", "eta", "chi", "phi"],
            ["nu", "delta"],
            does_not_raise(),
            id="psic-stage-lists",
        ),
        pytest.param(
            geometry_fourc,
            "fourc",
            ["omega", "chi", "phi"],
            ["two_theta"],
            does_not_raise(),
            id="fourc-stage-lists",
        ),
        pytest.param(
            geometry_sixc,
            "sixc",
            ["alpha", "omega", "chi", "phi"],
            ["delta", "gamma"],
            does_not_raise(),
            id="sixc-stage-lists",
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
        pytest.param(
            geometry_psic, "mu", None, does_not_raise(), id="psic-mu-on-floor"
        ),
        pytest.param(
            geometry_psic, "nu", None, does_not_raise(), id="psic-nu-on-floor"
        ),
        pytest.param(geometry_psic, "eta", "mu", does_not_raise(), id="psic-eta-on-mu"),
        pytest.param(
            geometry_psic, "chi", "eta", does_not_raise(), id="psic-chi-on-eta"
        ),
        pytest.param(
            geometry_psic, "phi", "chi", does_not_raise(), id="psic-phi-on-chi"
        ),
        pytest.param(
            geometry_psic, "delta", "nu", does_not_raise(), id="psic-delta-on-nu"
        ),
        pytest.param(
            geometry_fourc, "omega", None, does_not_raise(), id="fourc-omega-on-floor"
        ),
        pytest.param(
            geometry_fourc,
            "two_theta",
            None,
            does_not_raise(),
            id="fourc-two_theta-on-floor-decoupled",
        ),
        pytest.param(
            geometry_fourc, "chi", "omega", does_not_raise(), id="fourc-chi-on-omega"
        ),
        pytest.param(
            geometry_fourc, "phi", "chi", does_not_raise(), id="fourc-phi-on-chi"
        ),
        pytest.param(
            geometry_sixc, "alpha", None, does_not_raise(), id="sixc-alpha-on-floor"
        ),
        pytest.param(
            geometry_sixc, "omega", "alpha", does_not_raise(), id="sixc-omega-on-alpha"
        ),
        pytest.param(
            geometry_sixc, "delta", "alpha", does_not_raise(), id="sixc-delta-on-alpha"
        ),
        pytest.param(
            geometry_sixc, "gamma", "delta", does_not_raise(), id="sixc-gamma-on-delta"
        ),
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
            geometry_psic,
            "mu",
            +XHAT,
            does_not_raise(),
            id="psic-mu-vertical-right-handed",
        ),
        pytest.param(
            geometry_psic,
            "eta",
            -ZHAT,
            does_not_raise(),
            id="psic-eta-lateral-left-handed",
        ),
        pytest.param(
            geometry_psic,
            "chi",
            +YHAT,
            does_not_raise(),
            id="psic-chi-longitudinal-right-handed",
        ),
        pytest.param(
            geometry_psic,
            "phi",
            -ZHAT,
            does_not_raise(),
            id="psic-phi-lateral-left-handed",
        ),
        pytest.param(
            geometry_psic,
            "nu",
            +XHAT,
            does_not_raise(),
            id="psic-nu-vertical-right-handed",
        ),
        pytest.param(
            geometry_psic,
            "delta",
            -ZHAT,
            does_not_raise(),
            id="psic-delta-lateral-left-handed",
        ),
    ],
)
def test_geometry_axes(factory, stage_name, expected_axis, context):
    with context:
        g = factory()
        np.testing.assert_array_equal(g.stage(stage_name).axis, expected_axis)


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
