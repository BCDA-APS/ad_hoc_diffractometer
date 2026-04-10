"""
Unit tests for ad_hoc_diffractometer.geometry.

Covers:
  - AdHocDiffractometer construction and validation
  - Basis vector validation
  - Stacking order
  - set_angle() and stage()
  - sample_rotation_matrix() and detector_rotation_matrix()
  - check_limits()
  - azimuthal_reference property: storage, validation, default None
  - psi() method: psi=0 when n in scattering plane, psi=90 when n perp,
    uses current angles when called with no args, error cases
  - wh(print=False): str output, content, graceful fallbacks, stdout control
  - pa(print=False): str output, content, graceful fallbacks, stdout control
  - _spec_motor_name static helper: SPEC column-name mapping
"""

import math as _math
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


# ---------------------------------------------------------------------------
# inverse()
# ---------------------------------------------------------------------------

_LAMBDA_CU_KA = 1.5406

# psic angles for a convenient non-zero reflection
_PSIC_ANGLES = {
    "mu": 0.0,
    "eta": 20.97,
    "chi": 90.0,
    "phi": 0.0,
    "nu": 0.0,
    "delta": 41.94,
}


def _psic_with_identity_UB(wavelength=_LAMBDA_CU_KA):
    """Return a psic geometry with wavelength set and UB = I."""
    from ad_hoc_diffractometer import psic
    from ad_hoc_diffractometer import ub_identity

    g = psic()
    g.wavelength = wavelength
    ub_identity(g.sample)  # UB = B = I  (default cubic a=1 Å)
    return g


def test_inverse_all_zero_gives_zero():
    """All motor angles zero → no scattering → hkl = (0, 0, 0)."""
    g = _psic_with_identity_UB()
    hkl = g.inverse({"mu": 0, "eta": 0, "chi": 0, "phi": 0, "nu": 0, "delta": 0})
    assert hkl == pytest.approx((0.0, 0.0, 0.0), abs=1e-12)


def test_inverse_returns_tuple_of_floats():
    """inverse() returns a tuple of three Python floats."""
    g = _psic_with_identity_UB()
    hkl = g.inverse({"mu": 0, "eta": 0, "chi": 0, "phi": 0, "nu": 0, "delta": 0})
    assert isinstance(hkl, tuple)
    assert len(hkl) == 3
    assert all(isinstance(x, float) for x in hkl)


def test_inverse_UB_identity_equals_Q_phi():
    """
    When UB = I, inverse() must return exactly Q_phi (in Å⁻¹).

    UB = I implies hkl = UB⁻¹ @ Q_phi = I @ Q_phi = Q_phi.
    """
    from ad_hoc_diffractometer import angles_to_phi_vector
    from ad_hoc_diffractometer import psic

    g = psic()
    g.wavelength = _LAMBDA_CU_KA
    g.sample.UB = np.eye(3)

    Q_phi = angles_to_phi_vector(g, **_PSIC_ANGLES)
    hkl = g.inverse(_PSIC_ANGLES)

    np.testing.assert_allclose(hkl, Q_phi, atol=1e-10)


def test_inverse_UB_scaled_identity():
    """
    When UB = s·I, hkl = Q_phi / s.

    Scaling UB by a scalar s means hkl = (1/s)·Q_phi.
    """
    from ad_hoc_diffractometer import angles_to_phi_vector
    from ad_hoc_diffractometer import psic

    s = 2.5
    g = psic()
    g.wavelength = _LAMBDA_CU_KA
    g.sample.UB = s * np.eye(3)

    Q_phi = angles_to_phi_vector(g, **_PSIC_ANGLES)
    hkl = g.inverse(_PSIC_ANGLES)

    np.testing.assert_allclose(hkl, np.array(Q_phi) / s, atol=1e-10)


def test_inverse_round_trip_ub_identity():
    """
    Round-trip: UB @ hkl == Q_phi for the result of inverse().

    Whatever hkl inverse() returns, plugging it back into UB @ hkl
    must recover Q_phi.
    """
    from ad_hoc_diffractometer import angles_to_phi_vector
    from ad_hoc_diffractometer import psic
    from ad_hoc_diffractometer import ub_identity

    g = psic()
    g.wavelength = _LAMBDA_CU_KA
    ub_identity(g.sample)

    Q_phi = angles_to_phi_vector(g, **_PSIC_ANGLES)
    hkl = g.inverse(_PSIC_ANGLES)

    np.testing.assert_allclose(g.sample.UB @ np.array(hkl), Q_phi, atol=1e-10)


def test_inverse_round_trip_ub_from_one_reflection():
    """
    Round-trip with a physically meaningful UB from ub_from_one_reflection.

    UB is set from a sapphire (006) reflection; calling inverse() on the
    same angles must satisfy UB @ hkl == Q_phi.
    """
    from ad_hoc_diffractometer import Lattice
    from ad_hoc_diffractometer import angles_to_phi_vector
    from ad_hoc_diffractometer import psic
    from ad_hoc_diffractometer import ub_from_one_reflection

    g = psic()
    g.wavelength = _LAMBDA_CU_KA
    g.add_sample("sap", Lattice(a=4.758, c=12.991))
    g.sample = "sap"
    g.add_reflection("r1", hkl=(0, 0, 6), angles=_PSIC_ANGLES)
    ub_from_one_reflection(
        g.sample,
        "r1",
        reference_hkl=(0, 0, 1),
        reference_stage=g.stage("phi"),
    )

    Q_phi = angles_to_phi_vector(g, **_PSIC_ANGLES)
    hkl = g.inverse(_PSIC_ANGLES)

    np.testing.assert_allclose(g.sample.UB @ np.array(hkl), Q_phi, atol=1e-10)


def test_inverse_partial_angles_uses_current():
    """Stages not in the dict keep their current angle."""
    from ad_hoc_diffractometer import psic

    g = psic()
    g.wavelength = _LAMBDA_CU_KA
    g.sample.UB = np.eye(3)

    # Pre-set all angles on the geometry
    for name, angle in _PSIC_ANGLES.items():
        g.set_angle(name, angle)

    # Call with empty dict — uses the pre-set values
    hkl_empty = g.inverse({})
    hkl_explicit = g.inverse(_PSIC_ANGLES)

    np.testing.assert_allclose(hkl_empty, hkl_explicit, atol=1e-12)


def test_inverse_restores_stage_angles():
    """Motor angles are restored to their original values after inverse()."""
    g = _psic_with_identity_UB()
    g.set_angle("eta", 99.9)
    g.set_angle("phi", 45.0)

    g.inverse({"mu": 0, "eta": 20.97, "chi": 90, "phi": 0, "nu": 0, "delta": 41.94})

    assert g.stage("eta").angle == pytest.approx(99.9)
    assert g.stage("phi").angle == pytest.approx(45.0)


def test_inverse_no_wavelength_raises():
    """Raises ValueError when wavelength is not set."""
    from ad_hoc_diffractometer import psic
    from ad_hoc_diffractometer import ub_identity

    g = psic()
    ub_identity(g.sample)
    assert g.wavelength is None

    with pytest.raises(ValueError, match=re.escape("wavelength must be set")):
        g.inverse({"mu": 0, "eta": 0, "chi": 0, "phi": 0, "nu": 0, "delta": 0})


def test_inverse_no_UB_raises():
    """Raises ValueError when the active sample has no UB matrix."""
    g = _psic_with_identity_UB()
    g.sample.UB = None  # clear it

    with pytest.raises(ValueError, match=re.escape("no UB matrix")):
        g.inverse({"mu": 0, "eta": 0, "chi": 0, "phi": 0, "nu": 0, "delta": 0})


def test_inverse_unknown_stage_raises():
    """Raises KeyError for a stage name not in the geometry."""
    g = _psic_with_identity_UB()
    with pytest.raises(KeyError):
        g.inverse({"no_such_stage": 0.0})


# ---------------------------------------------------------------------------
# azimuthal_reference property (#11)
# ---------------------------------------------------------------------------


def test_azimuthal_reference_default_is_none():
    """azimuthal_reference is None by default."""
    from ad_hoc_diffractometer import fourcv

    assert fourcv().azimuthal_reference is None


def test_azimuthal_reference_set_tuple():
    """Setting to a 3-tuple stores it as (float, float, float)."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.azimuthal_reference = (0, 0, 1)
    assert g.azimuthal_reference == (0.0, 0.0, 1.0)


def test_azimuthal_reference_set_list():
    """Setting to a list of 3 numbers works (converts to tuple of floats)."""
    from ad_hoc_diffractometer import psic

    g = psic()
    g.azimuthal_reference = [1, 1, 0]
    assert g.azimuthal_reference == (1.0, 1.0, 0.0)


def test_azimuthal_reference_clear_with_none():
    """Setting to None clears the reference."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.azimuthal_reference = (0, 0, 1)
    g.azimuthal_reference = None
    assert g.azimuthal_reference is None


def test_azimuthal_reference_constructor():
    """azimuthal_reference can be set at construction via keyword argument."""
    from ad_hoc_diffractometer import fourcv

    # Use fourcv factory result and check the property works after setting
    g = fourcv()
    g.azimuthal_reference = (0, 1, 0)
    assert g.azimuthal_reference == (0.0, 1.0, 0.0)


def test_azimuthal_reference_zero_vector_raises():
    """Setting to (0, 0, 0) raises ValueError."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    with pytest.raises(ValueError, match=re.escape("non-zero")):
        g.azimuthal_reference = (0, 0, 0)


def test_azimuthal_reference_bad_type_raises():
    """Setting to a non-sequence raises ValueError."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    with pytest.raises(ValueError, match=re.escape("length-3 sequence")):
        g.azimuthal_reference = 42


def test_azimuthal_reference_wrong_length_raises():
    """Setting to a 2-element tuple raises ValueError."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    with pytest.raises(ValueError):
        g.azimuthal_reference = (0, 1)


# ---------------------------------------------------------------------------
# psi() method (#11)
# ---------------------------------------------------------------------------
#
# Setup: fourcv, cubic a=1 (B=I), lambda=2*pi, UB=I.
# At chi=0, phi=0, Q = XHAT (lateral in BL convention).
# Scattering plane = span(YHAT_BL, XHAT_BL) = lateral-longitudinal plane.
#
# n=(0,1,0)=YHAT_BL: lies in scattering plane → psi=0.
# n=(0,0,1)=ZHAT_BL: perpendicular to scattering plane → psi=90.


def _fourcv_identity():
    """fourcv with B=I (a=1), lambda=2pi, UB=I, azimuthal_reference=(0,0,1)."""
    from ad_hoc_diffractometer import Lattice
    from ad_hoc_diffractometer import fourcv
    from ad_hoc_diffractometer import ub_identity

    g = fourcv()
    g.wavelength = 2 * _math.pi
    g.sample.lattice = Lattice(a=1.0)
    ub_identity(g.sample)
    g.azimuthal_reference = (0, 0, 1)  # ZHAT_BL = vertical ⊥ scattering plane at chi=0
    return g


def test_psi_returns_float():
    """psi() returns a float."""
    g = _fourcv_identity()
    g.set_angle("omega", 30)
    g.set_angle("two_theta", 60)
    result = g.psi()
    assert isinstance(result, float)


def test_psi_n_perpendicular_to_scattering_plane_is_90():
    """n=(0,0,1) ⊥ scattering plane at chi=0 → psi=90."""
    g = _fourcv_identity()
    g.azimuthal_reference = (0, 0, 1)  # ZHAT_BL ⊥ lateral-longitudinal plane
    angles = {"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0}
    psi = g.psi(angles)
    assert abs(psi - 90.0) < 1e-8


def test_psi_n_in_scattering_plane_is_0():
    """n=(0,1,0) = YHAT_BL lies in scattering plane at chi=0 → psi=0."""
    g = _fourcv_identity()
    g.azimuthal_reference = (0, 1, 0)  # YHAT_BL = longitudinal = in scattering plane
    angles = {"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0}
    psi = g.psi(angles)
    assert abs(psi - 0.0) < 1e-8


def test_psi_uses_current_angles_when_none_passed():
    """psi() with no args uses the current motor angles."""
    g = _fourcv_identity()
    g.set_angle("omega", 30.0)
    g.set_angle("chi", 0.0)
    g.set_angle("phi", 0.0)
    g.set_angle("two_theta", 60.0)
    g.azimuthal_reference = (0, 0, 1)
    psi_implicit = g.psi()
    psi_explicit = g.psi({"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0})
    assert abs(psi_implicit - psi_explicit) < 1e-10


def test_psi_range_within_180():
    """psi() always returns a value in (-180, 180]."""
    g = _fourcv_identity()
    for chi in range(0, 180, 30):
        for phi in range(0, 360, 45):
            angles = {
                "omega": 30.0,
                "chi": float(chi),
                "phi": float(phi),
                "two_theta": 60.0,
            }
            try:
                psi = g.psi(angles)
                assert -180.0 < psi <= 180.0 or abs(abs(psi) - 180.0) < 1e-8
            except ValueError:
                pass  # n || Q at some positions — expected


def test_psi_no_reference_raises():
    """psi() raises ValueError when azimuthal_reference is None."""
    g = _fourcv_identity()
    g.azimuthal_reference = None
    with pytest.raises(ValueError, match=re.escape("azimuthal_reference")):
        g.psi({"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0})


def test_psi_no_wavelength_raises():
    """psi() raises ValueError when wavelength is not set."""
    g = _fourcv_identity()
    g.wavelength = None
    with pytest.raises(ValueError, match=re.escape("wavelength")):
        g.psi({"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0})


def test_psi_no_ub_raises():
    """psi() raises ValueError when sample.UB is None."""
    g = _fourcv_identity()
    g.sample.UB = None
    with pytest.raises(ValueError, match=re.escape("UB matrix")):
        g.psi({"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0})


def test_psi_n_parallel_to_q_raises():
    """psi() raises ValueError when the reference is parallel to Q."""
    g = _fourcv_identity()
    # At chi=0 Q is along XHAT_BL=(1,0,0). Set n=(1,0,0) → parallel to Q.
    g.azimuthal_reference = (1, 0, 0)
    with pytest.raises(ValueError, match=re.escape("parallel to Q")):
        g.psi({"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0})


def test_psi_pa_shows_azimuthal_reference():
    """pa property output includes the azimuthal reference hkl when set."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.wavelength = 1.5406
    g.azimuthal_reference = (0, 0, 1)
    assert "Azimuthal Reference" in g.pa(print=False)
    assert "0  0  1" in g.pa(print=False)


def test_psi_pa_shows_not_set_when_none():
    """pa method shows 'not set' for azimuthal reference when not configured."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.wavelength = 1.5406
    assert "not set" in g.pa(print=False)


def test_psi_wh_shows_psi_when_available():
    """wh method includes a Psi line when psi can be computed."""
    g = _fourcv_identity()
    g.set_angle("omega", 30.0)
    g.set_angle("chi", 0.0)
    g.set_angle("phi", 0.0)
    g.set_angle("two_theta", 60.0)
    assert "Psi" in g.wh(print=False)


def test_psi_wh_shows_not_available_when_no_reference():
    """wh method shows 'not available' for Psi when reference is not set."""
    from ad_hoc_diffractometer import Lattice
    from ad_hoc_diffractometer import fourcv
    from ad_hoc_diffractometer import ub_identity

    g = fourcv()
    g.wavelength = 2 * _math.pi
    g.sample.lattice = Lattice(a=1.0)
    ub_identity(g.sample)
    out = g.wh(print=False)
    assert "Psi" in out
    assert "not available" in out


# ---------------------------------------------------------------------------
# _spec_motor_name static helper
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "internal,expected",
    [
        pytest.param("two_theta", "TwoTheta", id="two_theta"),
        pytest.param("omega", "Theta", id="omega"),
        pytest.param("chi", "Chi", id="chi"),
        pytest.param("phi", "Phi", id="phi"),
        pytest.param("mu", "Mu", id="mu"),
        pytest.param("eta", "Eta", id="eta"),
        pytest.param("nu", "Nu", id="nu"),
        pytest.param("delta", "Delta", id="delta"),
        pytest.param("unknown_stage", "unknown_stage", id="unknown"),
    ],
)
def test_spec_motor_name(internal, expected):
    """_spec_motor_name is a static method on AdHocDiffractometer."""
    assert AdHocDiffractometer._spec_motor_name(internal) == expected


# ---------------------------------------------------------------------------
# wh() and pa() methods
# ---------------------------------------------------------------------------


class TestWhMethod:
    """Tests for AdHocDiffractometer.wh(print=False) — capture mode."""

    def test_wh_is_callable(self):
        from ad_hoc_diffractometer import fourcv

        assert callable(fourcv().wh)

    def test_wh_returns_str(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert isinstance(g.wh(print=False), str)

    def test_wh_contains_hkl(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "H K L" in g.wh(print=False)

    def test_wh_contains_lambda(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "1.5406" in g.wh(print=False)

    def test_wh_contains_motor_table(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        out = g.wh(print=False)
        assert "TwoTheta" in out
        assert "Theta" in out

    def test_wh_graceful_no_ub(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "not available" in g.wh(print=False)

    def test_wh_graceful_no_wavelength(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        out = g.wh(print=False)
        assert "H K L" in out
        assert "Lambda" in out

    def test_wh_print_true_prints(self, capsys):
        """Default print=True writes to stdout."""
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        g.wh()
        captured = capsys.readouterr()
        assert "Lambda" in captured.out

    def test_wh_print_false_no_stdout(self, capsys):
        """print=False produces no stdout output."""
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        g.wh(print=False)
        assert capsys.readouterr().out == ""

    def test_wh_psic_stage_names(self):
        from ad_hoc_diffractometer import psic

        g = psic()
        g.wavelength = 1.5406
        out = g.wh(print=False)
        assert "Mu" in out
        assert "Eta" in out
        assert "Nu" in out
        assert "Delta" in out


class TestPaMethod:
    """Tests for AdHocDiffractometer.pa(print=False) — capture mode."""

    def test_pa_is_callable(self):
        from ad_hoc_diffractometer import fourcv

        assert callable(fourcv().pa)

    def test_pa_returns_str(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert isinstance(g.pa(print=False), str)

    def test_pa_contains_geometry_name(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "fourcv" in g.pa(print=False)

    def test_pa_contains_lattice_section(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "Lattice" in g.pa(print=False)

    def test_pa_contains_wavelength(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "1.5406" in g.pa(print=False)

    def test_pa_no_reflections_shows_not_set(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        assert "not set" in g.pa(print=False)

    def test_pa_graceful_no_wavelength(self):
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        out = g.pa(print=False)
        assert "fourcv" in out
        assert "Lattice" in out

    def test_pa_print_true_prints(self, capsys):
        """Default print=True writes to stdout."""
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        g.pa()
        captured = capsys.readouterr()
        assert "Geometry" in captured.out

    def test_pa_print_false_no_stdout(self, capsys):
        """print=False produces no stdout output."""
        from ad_hoc_diffractometer import fourcv

        g = fourcv()
        g.wavelength = 1.5406
        g.pa(print=False)
        assert capsys.readouterr().out == ""

    def test_pa_psic_geometry_name(self):
        from ad_hoc_diffractometer import psic

        g = psic()
        g.wavelength = 1.5406
        assert "psic" in g.pa(print=False)


# ---------------------------------------------------------------------------
# Geometry — additional branch coverage
# ---------------------------------------------------------------------------


def test_geometry_repr():
    """AdHocDiffractometer.__repr__ includes the geometry name and stage names."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    r = repr(g)
    assert "fourcv" in r
    assert "omega" in r


def test_geometry_stacking_order_with_scrambled_input():
    """_ordered_stages correctly sorts stages even when passed in reverse parent order."""
    from ad_hoc_diffractometer import AdHocDiffractometer

    phi = Stage("phi", ZHAT, role="sample", parent=None)
    chi = Stage("chi", XHAT, role="sample", parent="phi")
    omega = Stage("omega", YHAT, role="sample", parent="chi")
    det = Stage("delta", ZHAT, role="detector", parent=None)
    # Pass in scrambled order — requires multiple while-loop passes in _ordered_stages
    g = AdHocDiffractometer("test", [omega, chi, phi, det])
    assert [s.name for s in g.sample_stages] == ["phi", "chi", "omega"]


def test_geometry_cycle_in_parent_chain_raises():
    """Stage parent chain with a cycle raises ValueError."""
    from ad_hoc_diffractometer import AdHocDiffractometer

    a = Stage("a", XHAT, role="sample", parent="b")
    b = Stage("b", ZHAT, role="sample", parent="a")
    with pytest.raises(ValueError, match="Cycle"):
        AdHocDiffractometer("test", [a, b])


def test_geometry_samples_property_is_read_only():
    """g._samples is a read-only property; direct assignment raises AttributeError."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    with pytest.raises(AttributeError):
        g._samples = {}


@pytest.mark.parametrize(
    "angles, match",
    [
        pytest.param(
            {"omega": 0.0, "chi": 0.0, "phi": 0.0, "two_theta": 0.0},
            "Q = 0",
            id="psi-q-zero",
        ),
    ],
)
def test_psi_q_zero_raises(angles, match):
    """psi() raises when the motor angles produce Q=0."""
    from ad_hoc_diffractometer import fourcv
    from ad_hoc_diffractometer import ub_identity

    g = fourcv()
    g.wavelength = 2 * _math.pi
    ub_identity(g.sample)
    g.azimuthal_reference = (0, 0, 1)
    with pytest.raises(ValueError, match=match):
        g.psi(angles)


def test_psi_n_maps_to_zero_raises():
    """psi() raises when UB @ n_hkl is zero."""
    from ad_hoc_diffractometer import Lattice
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.wavelength = 2 * _math.pi
    g.sample.lattice = Lattice(a=1.0)
    g.sample.UB = np.zeros((3, 3))
    g.sample.U = np.zeros((3, 3))
    g.azimuthal_reference = (0, 0, 1)
    with pytest.raises(ValueError, match="zero in the phi frame"):
        g.psi({"omega": 30.0, "chi": 0.0, "phi": 0.0, "two_theta": 60.0})


def test_pa_reflection_wavelength_none_falls_back_to_geometry_wavelength():
    """pa() uses the geometry wavelength when a reflection has wavelength=None."""
    from ad_hoc_diffractometer import fourcv
    from ad_hoc_diffractometer.reflection import Reflection

    g = fourcv()
    g.wavelength = 1.5406
    r = Reflection(
        "r1",
        hkl=(0, 0, 6),
        angles={"omega": 20.97, "chi": 90.0, "phi": 0.0, "two_theta": 41.94},
        wavelength=None,
    )
    g.sample.reflections._data["r1"] = r
    g.sample.reflections.setor1("r1")
    assert "1.5406" in g.pa(print=False)


def test_geometry_summary_with_wavelength(capsys):
    """summary() prints the geometry name and wavelength."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.wavelength = 1.5406
    g.summary()
    out = capsys.readouterr().out
    assert "fourcv" in out
    assert "1.5406" in out


def test_geometry_summary_no_wavelength(capsys):
    """summary() reports 'not set' when wavelength is unset."""
    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.summary()
    assert "not set" in capsys.readouterr().out


# ---------------------------------------------------------------------------
# Logging — wh() emits DEBUG when HKL / psi cannot be computed
# ---------------------------------------------------------------------------


def test_wh_logs_debug_when_hkl_unavailable(caplog):
    """wh() emits a DEBUG log when HKL cannot be computed (no UB set)."""
    import logging

    from ad_hoc_diffractometer import fourcv

    g = fourcv()
    g.wavelength = 1.5406
    with caplog.at_level(logging.DEBUG, logger="ad_hoc_diffractometer.geometry"):
        g.wh(print=False)

    assert any("wh" in r.message.lower() for r in caplog.records)
    assert any(r.levelno == logging.DEBUG for r in caplog.records)


def test_wh_logs_debug_when_psi_unavailable(caplog):
    """wh() emits a DEBUG log when psi cannot be computed (no azimuthal reference)."""
    import logging

    from ad_hoc_diffractometer import fourcv
    from ad_hoc_diffractometer import ub_identity

    g = fourcv()
    g.wavelength = 1.5406
    ub_identity(g.sample)
    # No azimuthal_reference → psi() raises → logger.debug fires
    with caplog.at_level(logging.DEBUG, logger="ad_hoc_diffractometer.geometry"):
        g.wh(print=False)

    assert any("psi" in r.message.lower() for r in caplog.records)
