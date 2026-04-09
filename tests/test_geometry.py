"""
Unit tests for ad_hoc_diffractometer.geometry.

Covers:
  - AdHocDiffractometer construction and validation
  - Basis vector validation
  - Stacking order
  - set_angle() and stage()
  - sample_rotation_matrix() and detector_rotation_matrix()
  - check_limits()
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import Rx
from helpers import Ry
from helpers import Rz

from ad_hoc_diffractometer import XHAT
from ad_hoc_diffractometer import YHAT
from ad_hoc_diffractometer import ZHAT
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import Stage

_VALID_STAGES = [Stage("a", XHAT, parent=None, role="sample")]


# ---------------------------------------------------------------------------
# Basis vector validation
# ---------------------------------------------------------------------------


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
# Construction and validation
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
# Stacking order
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
# set_angle() and stage()
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
def test_set_angle(stage_name, angle, context, psic_geom):
    with context:
        psic_geom.set_angle(stage_name, angle)
        assert psic_geom.stage(stage_name).angle == angle


# ---------------------------------------------------------------------------
# sample_rotation_matrix() — psic
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
def test_psic_sample_rotation(angles, expected_Z, context, psic_geom):
    with context:
        for name, val in angles.items():
            psic_geom.set_angle(name, val)
        np.testing.assert_allclose(
            psic_geom.sample_rotation_matrix(), expected_Z, atol=1e-10
        )


# ---------------------------------------------------------------------------
# detector_rotation_matrix() — psic
# ---------------------------------------------------------------------------


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
def test_psic_detector_rotation(angles, expected_D, context, psic_geom):
    with context:
        for name, val in angles.items():
            psic_geom.set_angle(name, val)
        np.testing.assert_allclose(
            psic_geom.detector_rotation_matrix(), expected_D, atol=1e-10
        )


# ---------------------------------------------------------------------------
# check_limits()
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "stage_limits, angles, context",
    [
        pytest.param(
            {},
            {"mu": 0.0, "eta": 0.0, "chi": 0.0, "phi": 0.0},
            does_not_raise(),
            id="all-zero-default-limits",
        ),
        pytest.param(
            {},
            {"mu": 179.9, "eta": -179.9},
            does_not_raise(),
            id="near-boundary-inside",
        ),
        pytest.param(
            {},
            {"mu": 180.0},
            does_not_raise(),
            id="exactly-at-max-boundary",
        ),
        pytest.param(
            {},
            {"mu": -180.0},
            does_not_raise(),
            id="exactly-at-min-boundary",
        ),
        pytest.param(
            {},
            {"mu": 181.0},
            pytest.raises(ValueError, match=re.escape("outside their limits")),
            id="single-violation-above-max",
        ),
        pytest.param(
            {},
            {"mu": -181.0},
            pytest.raises(ValueError, match=re.escape("outside their limits")),
            id="single-violation-below-min",
        ),
        pytest.param(
            {"mu": (-90.0, 90.0)},
            {"mu": 91.0},
            pytest.raises(ValueError, match=re.escape("outside their limits")),
            id="custom-limit-violated",
        ),
        pytest.param(
            {"mu": (-90.0, 90.0)},
            {"mu": 90.0},
            does_not_raise(),
            id="custom-limit-at-boundary",
        ),
    ],
)
def test_check_limits(stage_limits, angles, context, psic_geom):
    for stage_name, lims in stage_limits.items():
        psic_geom.stage(stage_name).limits = lims
    with context:
        psic_geom.check_limits(**angles)


def test_check_limits_multiple_violations(psic_geom):
    """All violations are reported in a single ValueError."""
    psic_geom.stage("mu").limits = (-90.0, 90.0)
    psic_geom.stage("eta").limits = (-45.0, 45.0)
    with pytest.raises(ValueError) as exc_info:
        psic_geom.check_limits(mu=100.0, eta=50.0, chi=0.0)
    msg = str(exc_info.value)
    assert "'mu'" in msg
    assert "'eta'" in msg
    assert "'chi'" not in msg


def test_check_limits_unknown_stage_raises(psic_geom):
    """An unknown stage name must raise KeyError (not ValueError)."""
    with pytest.raises(KeyError):
        psic_geom.check_limits(bogus=10.0)


def test_check_limits_empty_call(psic_geom):
    """Calling check_limits with no arguments must not raise."""
    psic_geom.check_limits()


# ---------------------------------------------------------------------------
# wavelength attribute
# ---------------------------------------------------------------------------


def test_wavelength_default_is_none(psic_geom):
    """Default wavelength is None (unset)."""
    assert psic_geom.wavelength is None


@pytest.mark.parametrize(
    "value, context",
    [
        pytest.param(1.5406, does_not_raise(), id="cu-kalpha"),
        pytest.param(0.7107, does_not_raise(), id="mo-kalpha"),
        pytest.param(0.001, does_not_raise(), id="very-short"),
        pytest.param(100.0, does_not_raise(), id="very-long"),
        pytest.param(None, does_not_raise(), id="unset-none"),
        pytest.param(
            0.0,
            pytest.raises(ValueError, match=re.escape("must be > 0")),
            id="invalid-zero",
        ),
        pytest.param(
            -1.0,
            pytest.raises(ValueError, match=re.escape("must be > 0")),
            id="invalid-negative",
        ),
    ],
)
def test_wavelength_assignment(value, context, psic_geom):
    with context:
        psic_geom.wavelength = value
        assert psic_geom.wavelength == value


def test_wavelength_constructor():
    """wavelength can be supplied at construction time."""
    from ad_hoc_diffractometer import psic

    g = psic()
    # default
    assert g.wavelength is None
    # via setter after construction
    g.wavelength = 1.5406
    assert g.wavelength == pytest.approx(1.5406)


def test_wavelength_invalid_at_construction():
    """Supplying an invalid wavelength at construction raises ValueError."""
    from ad_hoc_diffractometer import XHAT
    from ad_hoc_diffractometer import Stage

    stages = [Stage("a", XHAT, role="sample")]
    with pytest.raises(ValueError, match=re.escape("must be > 0")):
        AdHocDiffractometer("test", stages, wavelength=-1.0)


def test_wavelength_reset_to_none(psic_geom):
    """Wavelength can be cleared back to None after being set."""
    psic_geom.wavelength = 1.5406
    psic_geom.wavelength = None
    assert psic_geom.wavelength is None


def test_wavelength_in_summary(psic_geom, capsys):
    """summary() reports wavelength when set and 'not set' when None."""
    psic_geom.summary()
    out = capsys.readouterr().out
    assert "not set" in out

    psic_geom.wavelength = 1.5406
    psic_geom.summary()
    out = capsys.readouterr().out
    assert "1.540600 Å" in out


# ---------------------------------------------------------------------------
# kappa_alpha_deg attribute
# ---------------------------------------------------------------------------


def test_kappa_alpha_deg_default_is_none(psic_geom):
    """Non-kappa geometries have kappa_alpha_deg = None."""
    assert psic_geom.kappa_alpha_deg is None


def test_kappa_alpha_deg_settable():
    """kappa_alpha_deg can be set on a bare AdHocDiffractometer."""
    from ad_hoc_diffractometer import XHAT
    from ad_hoc_diffractometer import Stage

    stages = [Stage("a", XHAT, role="sample")]
    g = AdHocDiffractometer("test", stages, kappa_alpha_deg=50.0)
    assert g.kappa_alpha_deg == pytest.approx(50.0)


def test_kappa_alpha_deg_none_for_non_kappa():
    """All non-kappa factories return None for kappa_alpha_deg."""
    from ad_hoc_diffractometer import fivec
    from ad_hoc_diffractometer import fourch
    from ad_hoc_diffractometer import fourcv
    from ad_hoc_diffractometer import psic
    from ad_hoc_diffractometer import s2d2
    from ad_hoc_diffractometer import sixc
    from ad_hoc_diffractometer import zaxis

    for factory in (psic, fourcv, fourch, sixc, zaxis, s2d2, fivec):
        g = factory()
        assert g.kappa_alpha_deg is None, f"{factory.__name__} should have None"
