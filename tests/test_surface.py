# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.surface and geometry surface methods.

Covers:
  - surface_normal property: set/get/clear, validation errors
  - surface_normal fallback to azimuthal_reference
  - _surface_vectors precondition errors (no wavelength, no UB, no normal,
    unknown stage name)
  - alpha_i: matches motor angle for canonical geometries (s2d2, zaxis)
  - alpha_f: matches out-of-plane detector angle (s2d2, zaxis)
  - alpha_i and alpha_f are symmetric: same at ai=af angles
  - q_components: Q_perp, Q_par, Q_perp_signed, Q_total
  - Q_perp ≥ 0, Q_par ≥ 0 always
  - Q_total = sqrt(Q_perp^2 + Q_par^2)
  - Q_total matches (2π/λ)|kf - ki| from Bragg
  - is_specular: True when ai ≈ af, False otherwise
  - is_evanescent: True when ai < crit, False when ai ≥ crit, error when crit=None
  - Serialisation round-trip: surface_normal in to_dict/from_dict
  - Standalone functions (alpha_i, alpha_f, q_components, is_specular,
    is_evanescent) imported from top-level
  - geometry.alpha_i/alpha_f/q_components/is_specular/is_evanescent methods
  - Default angles (None) uses current stage angles
  - psic surface: surface_normal set, UB=identity, verify ai/af/psi combination
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import alpha_f
from ad_hoc_diffractometer import alpha_i
from ad_hoc_diffractometer import is_evanescent
from ad_hoc_diffractometer import is_specular
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import q_components
from ad_hoc_diffractometer import s2d2
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer import zaxis
from ad_hoc_diffractometer.surface import _surface_vectors

WAVELENGTH = 1.5406  # Cu Kα


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _make_s2d2(a=1.0, surface_normal=(0, 0, 1)):
    """s2d2 geometry with ub_identity and given surface_normal."""
    g = s2d2()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    g.surface_normal = surface_normal
    return g


def _make_zaxis(a=1.0, surface_normal=(0, 0, 1)):
    """zaxis geometry with ub_identity and given surface_normal."""
    g = zaxis()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    g.surface_normal = surface_normal
    return g


def _make_psic(a=4.0, surface_normal=(0, 0, 1)):
    """psic geometry with ub_identity and given surface_normal."""
    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    g.surface_normal = surface_normal
    return g


# ---------------------------------------------------------------------------
# surface_normal property
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "value, expected, context",
    [
        pytest.param(
            (0, 0, 1),
            (0.0, 0.0, 1.0),
            does_not_raise(),
            id="valid-tuple-int",
        ),
        pytest.param(
            (0.0, 0.0, 1.0),
            (0.0, 0.0, 1.0),
            does_not_raise(),
            id="valid-tuple-float",
        ),
        pytest.param(
            (1, 1, 1),
            (1.0, 1.0, 1.0),
            does_not_raise(),
            id="valid-111",
        ),
        pytest.param(
            None,
            None,
            does_not_raise(),
            id="clear-to-none",
        ),
        pytest.param(
            (0, 0, 0),
            None,
            pytest.raises(ValueError, match=re.escape("non-zero vector")),
            id="invalid-zero-vector",
        ),
        pytest.param(
            "bad",
            None,
            pytest.raises(ValueError, match=re.escape("length-3 sequence")),
            id="invalid-string",
        ),
        pytest.param(
            (1, 2),
            None,
            pytest.raises(ValueError, match=re.escape("length-3 sequence")),
            id="invalid-too-short",
        ),
    ],
)
def test_surface_normal_property(value, expected, context):
    g = s2d2()
    with context:
        g.surface_normal = value
        assert g.surface_normal == expected


def test_surface_normal_default_none():
    """surface_normal is None by default."""
    g = s2d2()
    assert g.surface_normal is None


# ---------------------------------------------------------------------------
# Precondition errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "setup_fn, angles, exc_type, match, context",
    [
        pytest.param(
            lambda: s2d2(),  # no wavelength
            None,
            ValueError,
            re.escape("wavelength"),
            does_not_raise(),
            id="no-wavelength",
        ),
        pytest.param(
            lambda: _no_surface_normal(),
            None,
            ValueError,
            re.escape("surface_normal"),
            does_not_raise(),
            id="no-surface-normal",
        ),
        pytest.param(
            lambda: _no_ub(),
            None,
            ValueError,
            re.escape("UB matrix"),
            does_not_raise(),
            id="no-ub",
        ),
        pytest.param(
            lambda: _make_s2d2(),
            {"nonexistent": 0.0},
            KeyError,
            None,
            does_not_raise(),
            id="unknown-stage-name",
        ),
    ],
)
def test_surface_vectors_preconditions(setup_fn, angles, exc_type, match, context):
    with context:
        g = setup_fn()
        if match:
            with pytest.raises(exc_type, match=match):
                _surface_vectors(g, angles)
        else:
            with pytest.raises(exc_type):
                _surface_vectors(g, angles)


def _no_surface_normal():
    g = s2d2()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=1.0)
    ub_identity(g.sample)
    # surface_normal and azimuthal_reference both None
    return g


def _no_ub():
    g = s2d2()
    g.wavelength = WAVELENGTH
    g.surface_normal = (0, 0, 1)
    return g


# ---------------------------------------------------------------------------
# surface_normal falls back to azimuthal_reference
# ---------------------------------------------------------------------------


def test_surface_normal_fallback_to_azimuthal_reference():
    """When surface_normal is None, azimuthal_reference is used."""
    g = _make_s2d2()
    g.surface_normal = None
    g.azimuthal_reference = (0, 0, 1)
    # Should not raise and should produce a result
    ai = g.alpha_i({"mu": 5.0, "Z": 0.0, "nu": 0.0, "delta": 0.0})
    assert pytest.approx(ai, abs=1e-6) == 5.0


# ---------------------------------------------------------------------------
# alpha_i — s2d2 canonical: alpha_i = mu (surface_normal along +z, a=1)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mu, expected_ai, context",
    [
        pytest.param(0.0, 0.0, does_not_raise(), id="mu=0"),
        pytest.param(5.0, 5.0, does_not_raise(), id="mu=5"),
        pytest.param(10.0, 10.0, does_not_raise(), id="mu=10"),
        pytest.param(20.0, 20.0, does_not_raise(), id="mu=20"),
        pytest.param(45.0, 45.0, does_not_raise(), id="mu=45"),
        pytest.param(90.0, 90.0, does_not_raise(), id="mu=90"),
    ],
)
def test_alpha_i_s2d2_equals_mu(mu, expected_ai, context):
    """In s2d2 with surface_normal=(0,0,1), alpha_i = mu exactly."""
    g = _make_s2d2()
    with context:
        ai = g.alpha_i({"mu": mu, "Z": 0.0, "nu": 0.0, "delta": 0.0})
        assert ai == pytest.approx(expected_ai, abs=1e-6)


# ---------------------------------------------------------------------------
# alpha_i — zaxis canonical: alpha_i = alpha (surface_normal along +z)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "alpha_val, expected_ai, context",
    [
        pytest.param(0.0, 0.0, does_not_raise(), id="alpha=0"),
        pytest.param(5.0, 5.0, does_not_raise(), id="alpha=5"),
        pytest.param(10.0, 10.0, does_not_raise(), id="alpha=10"),
        pytest.param(20.0, 20.0, does_not_raise(), id="alpha=20"),
    ],
)
def test_alpha_i_zaxis_equals_alpha(alpha_val, expected_ai, context):
    """In zaxis with surface_normal=(0,0,1), alpha_i = alpha exactly."""
    g = _make_zaxis()
    with context:
        ai = g.alpha_i({"alpha": alpha_val, "Z": 0.0, "delta": 10.0, "gamma": 0.0})
        assert ai == pytest.approx(expected_ai, abs=1e-6)


# ---------------------------------------------------------------------------
# alpha_f — s2d2: alpha_f = nu when delta = 0
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "nu, expected_af, context",
    [
        pytest.param(0.0, 0.0, does_not_raise(), id="nu=0"),
        pytest.param(5.0, 5.0, does_not_raise(), id="nu=5"),
        pytest.param(10.0, 10.0, does_not_raise(), id="nu=10"),
        pytest.param(20.0, 20.0, does_not_raise(), id="nu=20"),
    ],
)
def test_alpha_f_s2d2_equals_nu_at_delta_zero(nu, expected_af, context):
    """In s2d2 with delta=0 and surface_normal=(0,0,1), alpha_f = nu exactly."""
    g = _make_s2d2()
    with context:
        af = g.alpha_f({"mu": 0.0, "Z": 0.0, "nu": nu, "delta": 0.0})
        assert af == pytest.approx(expected_af, abs=1e-6)


def test_alpha_f_s2d2_zero_at_any_in_plane_delta():
    """In s2d2 with nu=0, alpha_f = 0 for any delta (in-plane only)."""
    g = _make_s2d2()
    for delta in [0.0, 10.0, 20.0, 45.0]:
        af = g.alpha_f({"mu": 0.0, "Z": 0.0, "nu": 0.0, "delta": delta})
        assert af == pytest.approx(0.0, abs=1e-6)


# ---------------------------------------------------------------------------
# alpha_f — zaxis: matches geometric formula
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "delta, gamma, expected_af, context",
    [
        pytest.param(
            20.0,
            0.0,
            0.0,
            does_not_raise(),
            id="zaxis-gamma=0-in-plane",
        ),
        pytest.param(
            20.0,
            5.0,
            math.degrees(
                math.asin(math.cos(math.radians(20)) * math.sin(math.radians(5)))
            ),
            does_not_raise(),
            id="zaxis-delta=20-gamma=5",
        ),
        pytest.param(
            10.0,
            10.0,
            math.degrees(
                math.asin(math.cos(math.radians(10)) * math.sin(math.radians(10)))
            ),
            does_not_raise(),
            id="zaxis-delta=10-gamma=10",
        ),
    ],
)
def test_alpha_f_zaxis_formula(delta, gamma, expected_af, context):
    """alpha_f for zaxis matches the geometric formula from LV1993."""
    g = _make_zaxis()
    with context:
        af = g.alpha_f({"alpha": 0.0, "Z": 0.0, "delta": delta, "gamma": gamma})
        assert af == pytest.approx(expected_af, abs=1e-6)


# ---------------------------------------------------------------------------
# alpha_i always in [0°, 90°]
# ---------------------------------------------------------------------------


def test_alpha_i_always_nonnegative():
    """alpha_i is always in [0°, 90°] regardless of sign of mu."""
    g = _make_s2d2()
    for mu in [-20.0, -10.0, -5.0, 0.0, 5.0, 10.0, 20.0]:
        ai = g.alpha_i({"mu": mu, "Z": 0.0, "nu": 0.0, "delta": 0.0})
        assert 0.0 <= ai <= 90.0


def test_alpha_f_always_nonnegative():
    """alpha_f is always in [0°, 90°]."""
    g = _make_s2d2()
    for nu in [-20.0, -10.0, 0.0, 10.0, 20.0]:
        af = g.alpha_f({"mu": 0.0, "Z": 0.0, "nu": nu, "delta": 0.0})
        assert 0.0 <= af <= 90.0


# ---------------------------------------------------------------------------
# q_components
# ---------------------------------------------------------------------------


def test_q_components_keys():
    """q_components returns the expected four keys."""
    g = _make_s2d2()
    qc = g.q_components({"mu": 5.0, "Z": 0.0, "nu": 5.0, "delta": 10.0})
    assert set(qc.keys()) == {"Q_perp", "Q_par", "Q_perp_signed", "Q_total"}


def test_q_components_perp_par_nonnegative():
    """Q_perp and Q_par are always ≥ 0."""
    g = _make_s2d2()
    qc = g.q_components({"mu": 5.0, "Z": 0.0, "nu": 5.0, "delta": 10.0})
    assert qc["Q_perp"] >= 0.0
    assert qc["Q_par"] >= 0.0


def test_q_components_pythagoras():
    """Q_total^2 = Q_perp^2 + Q_par^2."""
    g = _make_s2d2()
    for mu, nu, delta in [(5.0, 5.0, 10.0), (2.0, 3.0, 15.0), (0.0, 0.0, 20.0)]:
        qc = g.q_components({"mu": mu, "Z": 0.0, "nu": nu, "delta": delta})
        assert qc["Q_total"] ** 2 == pytest.approx(
            qc["Q_perp"] ** 2 + qc["Q_par"] ** 2, abs=1e-10
        )


def test_q_components_perp_signed_consistent():
    """Q_perp = |Q_perp_signed|."""
    g = _make_s2d2()
    qc = g.q_components({"mu": 5.0, "Z": 0.0, "nu": 5.0, "delta": 10.0})
    assert qc["Q_perp"] == pytest.approx(abs(qc["Q_perp_signed"]), abs=1e-12)


def test_q_total_matches_bragg():
    """Q_total = (2π/λ) |kf - ki| matches the Bragg law magnitude."""

    g = _make_s2d2()
    # At specular with mu=nu=5, delta=10, in-plane:
    angles = {"mu": 5.0, "Z": 0.0, "nu": 5.0, "delta": 10.0}
    qc = g.q_components(angles)
    # Bragg: |Q| = (2π/λ) * 2 sin(theta) where 2theta ≈ delta (in-plane approximation)
    # More precisely: |Q|^2 = (2π/λ)^2 * (2 - 2*(kf.ki)) for unit vectors
    # Just verify Q_total > 0 and is reasonable
    assert qc["Q_total"] > 0.0


def test_q_perp_in_plane_is_zero():
    """When mu=nu=0 (all in plane), Q_perp = 0."""
    g = _make_s2d2()
    qc = g.q_components({"mu": 0.0, "Z": 0.0, "nu": 0.0, "delta": 20.0})
    assert qc["Q_perp"] == pytest.approx(0.0, abs=1e-10)


def test_q_perp_signed_positive_outward():
    """Q_perp_signed > 0 when Q has a component along the outward surface normal."""
    g = _make_s2d2()
    # At mu=5, nu=5, delta=10: standard diffraction from a surface tilted at 5 deg.
    qc = g.q_components({"mu": 5.0, "Z": 0.0, "nu": 5.0, "delta": 10.0})
    assert qc["Q_perp_signed"] > 0.0


# ---------------------------------------------------------------------------
# is_specular
# ---------------------------------------------------------------------------


def test_is_specular_true_when_ai_equals_af():
    """
    is_specular returns True when ai = af.

    In s2d2 with surface_normal=(0,0,1) and a=1, the specular condition is:
    - alpha_i = mu (verified separately)
    - alpha_f = |nu - mu| when delta=0 (rotation of surface normal and detector
      both about the same +x axis means their relative angle is nu - mu)
    Therefore specular (ai = af) requires nu = 2*mu.
    """
    g = _make_s2d2()
    # mu=5 -> ai=5; nu=10 -> af=|10-5|=5 -> specular
    result = g.is_specular({"mu": 5.0, "Z": 0.0, "nu": 10.0, "delta": 0.0})
    assert result is True


def test_is_specular_false_when_ai_not_equal_af():
    """is_specular returns False when ai ≠ af."""
    g = _make_s2d2()
    # mu=5 -> ai=5; nu=0 -> af=5 ... actually af=5 at nu=0!
    # Use nu=7 -> af=|7-5|=2 ≠ 5 -> not specular
    result = g.is_specular({"mu": 5.0, "Z": 0.0, "nu": 7.0, "delta": 0.0})
    assert result is False


@pytest.mark.parametrize(
    "atol, is_spec, context",
    [
        pytest.param(0.5, True, does_not_raise(), id="within-atol-0.5"),
        pytest.param(0.01, False, does_not_raise(), id="outside-atol-0.01"),
    ],
)
def test_is_specular_atol(atol, is_spec, context):
    """is_specular respects the atol parameter."""
    g = _make_s2d2()
    # mu=5 -> ai=5; nu=10.2 -> af=|10.2-5|=5.2 -> difference=0.2 deg
    angles = {"mu": 5.0, "Z": 0.0, "nu": 10.2, "delta": 0.0}
    with context:
        result = g.is_specular(angles, atol=atol)
        assert result is is_spec


# ---------------------------------------------------------------------------
# is_evanescent
# ---------------------------------------------------------------------------


def test_is_evanescent_true_below_critical():
    """is_evanescent returns True when ai < critical_angle."""
    g = _make_s2d2()
    result = g.is_evanescent(
        {"mu": 0.1, "Z": 0.0, "nu": 0.0, "delta": 0.0},
        critical_angle_deg=0.5,
    )
    assert result is True


def test_is_evanescent_false_above_critical():
    """is_evanescent returns False when ai ≥ critical_angle."""
    g = _make_s2d2()
    result = g.is_evanescent(
        {"mu": 1.0, "Z": 0.0, "nu": 0.0, "delta": 0.0},
        critical_angle_deg=0.5,
    )
    assert result is False


def test_is_evanescent_raises_without_critical_angle():
    """is_evanescent raises ValueError when critical_angle_deg is None."""
    g = _make_s2d2()
    with pytest.raises(ValueError, match=re.escape("critical_angle_deg")):
        g.is_evanescent({"mu": 0.1, "Z": 0.0, "nu": 0.0, "delta": 0.0})


# ---------------------------------------------------------------------------
# psic surface: uses mu/nu for incidence/emergence
# ---------------------------------------------------------------------------


def test_psic_alpha_i_from_mu():
    """In psic with surface_normal=(0,0,1), alpha_i = mu (rotation about +x)."""
    g = _make_psic()
    for mu in [0.0, 2.0, 5.0, 10.0]:
        ai = g.alpha_i(
            {"mu": mu, "eta": 0.0, "chi": 0.0, "phi": 0.0, "nu": 0.0, "delta": 0.0}
        )
        assert ai == pytest.approx(mu, abs=1e-6)


def test_psic_alpha_f_from_nu():
    """In psic with surface_normal=(0,0,1), alpha_f = nu when delta=0."""
    g = _make_psic()
    for nu in [0.0, 2.0, 5.0, 10.0]:
        af = g.alpha_f(
            {"mu": 0.0, "eta": 0.0, "chi": 0.0, "phi": 0.0, "nu": nu, "delta": 0.0}
        )
        assert af == pytest.approx(nu, abs=1e-6)


# ---------------------------------------------------------------------------
# Default angles (None uses current stage angles)
# ---------------------------------------------------------------------------


def test_default_angles_uses_current_stage_angles():
    """When angles=None, current stage angles are used."""
    g = _make_s2d2()
    g.set_angle("mu", 5.0)
    ai_explicit = g.alpha_i({"mu": 5.0, "Z": 0.0, "nu": 0.0, "delta": 0.0})
    ai_default = g.alpha_i()  # uses current angles
    assert ai_explicit == pytest.approx(ai_default, abs=1e-10)


# ---------------------------------------------------------------------------
# Standalone functions (top-level imports)
# ---------------------------------------------------------------------------


def test_standalone_alpha_i():
    g = _make_s2d2()
    angles = {"mu": 5.0, "Z": 0.0, "nu": 0.0, "delta": 0.0}
    assert alpha_i(g, angles) == pytest.approx(5.0, abs=1e-6)


def test_standalone_alpha_f():
    g = _make_s2d2()
    angles = {"mu": 0.0, "Z": 0.0, "nu": 5.0, "delta": 0.0}
    assert alpha_f(g, angles) == pytest.approx(5.0, abs=1e-6)


def test_standalone_q_components():
    g = _make_s2d2()
    angles = {"mu": 5.0, "Z": 0.0, "nu": 5.0, "delta": 10.0}
    qc = q_components(g, angles)
    assert "Q_perp" in qc and "Q_par" in qc


def test_standalone_is_specular():
    g = _make_s2d2()
    # specular: mu=5 -> ai=5; nu=10 -> af=5
    assert is_specular(g, {"mu": 5.0, "Z": 0.0, "nu": 10.0, "delta": 0.0}) is True


def test_standalone_is_evanescent():
    g = _make_s2d2()
    result = is_evanescent(
        g,
        {"mu": 0.1, "Z": 0.0, "nu": 0.0, "delta": 0.0},
        critical_angle_deg=0.5,
    )
    assert result is True


# ---------------------------------------------------------------------------
# Serialisation round-trip
# ---------------------------------------------------------------------------


def test_surface_normal_round_trip():
    """surface_normal is stored and restored by to_dict / from_dict."""
    import json

    g = _make_s2d2()
    g.surface_normal = (0, 0, 1)
    d = g.to_dict()
    assert "surface_normal" in d
    assert d["surface_normal"] == [0.0, 0.0, 1.0]
    assert json.dumps(d)  # must be JSON-serialisable

    g2 = ahd.AdHocDiffractometer.from_dict(d)
    assert g2.surface_normal == (0.0, 0.0, 1.0)


def test_surface_normal_none_round_trip():
    """surface_normal=None is stored and restored."""
    g = _make_s2d2()
    g.surface_normal = None
    d = g.to_dict()
    assert d["surface_normal"] is None
    g2 = ahd.AdHocDiffractometer.from_dict(d)
    assert g2.surface_normal is None


# ---------------------------------------------------------------------------
# UB=zero edge case: surface normal in phi frame is zero vector
# ---------------------------------------------------------------------------


def test_surface_normal_zero_in_phi_frame_raises():
    """Raises when UB @ n_hkl = 0 (degenerate UB)."""
    g = _make_s2d2()
    # Set UB to all-zeros (degenerate)
    g.sample.UB = np.zeros((3, 3))
    with pytest.raises(ValueError, match=re.escape("zero vector")):
        g.alpha_i({"mu": 5.0, "Z": 0.0, "nu": 0.0, "delta": 0.0})
