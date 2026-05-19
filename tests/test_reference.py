# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Unit tests for ad_hoc_diffractometer.reference — reference pseudo-angles.

Covers:
  - incidence_angle: requires surface_normal; raises when None
  - exit_angle: requires surface_normal; raises when None
  - psi_angle: requires azimuthal_reference; raises when None
  - naz_angle: requires surface_normal; raises when None; vertical n̂ gives 0
  - ReferenceConstraint.is_implemented(): True when reference set, False when None
  - Serialisation round-trip with surface_normal and azimuthal_reference set
  - Smoke tests: reasonable output range for known geometry configurations
"""

import re

import pytest
from helpers import psic
from helpers import s2d2
from helpers import sixc
from helpers import zaxis

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer import ReferenceConstraint
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.reference import exit_angle
from ad_hoc_diffractometer.reference import incidence_angle
from ad_hoc_diffractometer.reference import naz_angle
from ad_hoc_diffractometer.reference import psi_angle

WAVELENGTH = 1.5406


def _setup_psic():
    """Return a psic geometry with wavelength, cubic lattice, and UB=B."""
    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ub_identity(g.sample)
    return g


# ---------------------------------------------------------------------------
# incidence_angle
# ---------------------------------------------------------------------------


def test_incidence_angle_raises_without_surface_normal():
    """incidence_angle raises ValueError when surface_normal is None."""
    g = _setup_psic()
    assert g.surface_normal is None
    with pytest.raises(ValueError, match=re.escape("surface_normal must be set")):
        incidence_angle(g)


def test_incidence_angle_with_surface_normal():
    """incidence_angle returns a float in [-90, 90] when surface_normal is set."""
    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    g.mode_name = "bisecting_vertical"
    sols = g.forward(1, 0, 0)
    assert len(sols) > 0
    for s in sols:
        ai = incidence_angle(g, angles=s)
        assert isinstance(ai, float)
        assert -90.0 <= ai <= 90.0


def test_incidence_angle_uses_current_angles_when_none():
    """incidence_angle uses current stage angles when angles=None."""
    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    ai = incidence_angle(g, angles=None)
    assert isinstance(ai, float)


# ---------------------------------------------------------------------------
# exit_angle
# ---------------------------------------------------------------------------


def test_exit_angle_raises_without_surface_normal():
    """exit_angle raises ValueError when surface_normal is None."""
    g = _setup_psic()
    with pytest.raises(ValueError, match=re.escape("surface_normal must be set")):
        exit_angle(g)


def test_exit_angle_with_surface_normal():
    """exit_angle returns a float in [-90, 90] when surface_normal is set."""
    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    g.mode_name = "bisecting_vertical"
    sols = g.forward(1, 0, 0)
    for s in sols:
        af = exit_angle(g, angles=s)
        assert isinstance(af, float)
        assert -90.0 <= af <= 90.0


def test_specular_condition_alpha_i_equals_alpha_f():
    """At bisecting with surface normal ⊥ to scattering plane, alpha_i ≈ alpha_f."""
    g = _setup_psic()
    # Surface normal along transverse axis — perpendicular to the scattering plane
    g.surface_normal = (0, 0, 1)
    g.mode_name = "bisecting_vertical"
    sols = g.forward(1, 0, 0)
    for s in sols:
        ai = incidence_angle(g, angles=s)
        af = exit_angle(g, angles=s)
        # At bisecting in vertical plane with transverse surface normal, ai ≈ af
        assert ai == pytest.approx(af, abs=1e-6)


# ---------------------------------------------------------------------------
# psi_angle
# ---------------------------------------------------------------------------


def test_psi_angle_raises_without_azimuthal_reference():
    """psi_angle raises ValueError when azimuthal_reference is None."""
    g = _setup_psic()
    assert g.azimuthal_reference is None
    with pytest.raises(ValueError, match=re.escape("azimuthal_reference must be set")):
        psi_angle(g)


def test_psi_angle_with_azimuthal_reference():
    """psi_angle returns a float in (-180, 180] when azimuthal_reference is set.

    Uses ``(0, 1, 0)`` instead of ``(1, 0, 0)``: under issue #280
    ub_identity the crystal a* axis is along the beam, so Q_phi(1,0,0)
    is parallel to the beam and psi is undefined.  ``(0, 1, 0)``
    produces a Q_phi off the beam axis.
    """
    g = _setup_psic()
    g.azimuthal_reference = (0, 0, 1)
    g.mode_name = "bisecting_vertical"
    sols = g.forward(0, 1, 0)
    for s in sols:
        psi = psi_angle(g, angles=s)
        assert isinstance(psi, float)
        assert -180.0 < psi <= 180.0


def test_psi_angle_uses_current_angles_when_none():
    """psi_angle uses current stage angles when angles=None.

    Uses ``(0, 1, 0)`` for the same reason as
    :func:`test_psi_angle_with_azimuthal_reference`.
    """
    g = _setup_psic()
    g.azimuthal_reference = (0, 0, 1)
    g.mode_name = "bisecting_vertical"
    sols = g.forward(0, 1, 0)
    s = sols[0]
    for name, value in s.items():
        g.set_angle(name, value)
    psi = psi_angle(g, angles=None)
    assert isinstance(psi, float)


# ---------------------------------------------------------------------------
# naz_angle
# ---------------------------------------------------------------------------


def test_naz_angle_raises_without_surface_normal():
    """naz_angle raises ValueError when surface_normal is None."""
    g = _setup_psic()
    with pytest.raises(ValueError, match=re.escape("surface_normal must be set")):
        naz_angle(g)


def test_naz_angle_with_surface_normal():
    """naz_angle returns a float in (-180, 180] when surface_normal is set."""
    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    naz = naz_angle(g)
    assert isinstance(naz, float)
    assert -180.0 < naz <= 180.0


def test_naz_angle_with_explicit_angles():
    """naz_angle accepts an explicit angles dict."""
    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    angles = {s.name: s.angle for s in g._stages.values()}
    naz = naz_angle(g, angles=angles)
    assert isinstance(naz, float)


def test_naz_angle_uses_current_angles_when_none():
    """naz_angle uses current stage angles when angles=None."""
    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    naz = naz_angle(g, angles=None)
    assert isinstance(naz, float)


def test_naz_angle_vertical_normal_returns_zero():
    """naz_angle returns 0 when surface normal is vertical (undefined by convention).

    Under issue #280 ub_identity, ``UB @ (0, 1, 0) = U[:, 1] · |b2*|`` is
    physically along ``+vertical`` (column 1 of U is the vertical basis
    vector).  Pre-#280 the same physical configuration was selected by
    ``surface_normal = (1, 0, 0)`` because the basis-relative ``U = I``
    placed the crystal a-axis along the basis-x direction (= vertical
    in psic-YOU).
    """
    g = _setup_psic()
    g.surface_normal = (0, 1, 0)
    naz = naz_angle(g)
    assert naz == 0.0


# ---------------------------------------------------------------------------
# ReferenceConstraint.is_implemented
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "name, value, ref_attr, ref_value, expected",
    [
        # surface-normal constraints: implemented when surface_normal is set
        pytest.param(
            "alpha_i", 0.0, "surface_normal", (0, 0, 1), True, id="alpha_i-with-sn"
        ),
        pytest.param("alpha_i", 0.0, "surface_normal", None, False, id="alpha_i-no-sn"),
        pytest.param(
            "beta_out", 0.0, "surface_normal", (0, 0, 1), True, id="beta_out-with-sn"
        ),
        pytest.param(
            "beta_out", 0.0, "surface_normal", None, False, id="beta_out-no-sn"
        ),
        pytest.param(
            "a_eq_b", True, "surface_normal", (0, 0, 1), True, id="a_eq_b-with-sn"
        ),
        pytest.param("a_eq_b", True, "surface_normal", None, False, id="a_eq_b-no-sn"),
        # psi: implemented when azimuthal_reference is set
        pytest.param(
            "psi",
            0.0,
            "azimuthal_reference",
            (0, 0, 1),
            True,
            id="psi-with-azref",
        ),
        pytest.param(
            "psi",
            0.0,
            "azimuthal_reference",
            None,
            False,
            id="psi-no-azref",
        ),
        # naz: not yet implemented regardless of reference vector
        pytest.param(
            "naz",
            0.0,
            "azimuthal_reference",
            (0, 0, 1),
            False,
            id="naz-not-implemented",
        ),
    ],
)
def test_reference_constraint_is_implemented(
    name, value, ref_attr, ref_value, expected
):
    """ReferenceConstraint.is_implemented() reflects solver availability."""
    g = _setup_psic()
    setattr(g, ref_attr, ref_value)
    rc = ReferenceConstraint(name, value)
    assert rc.is_implemented(g) is expected


@pytest.mark.parametrize(
    "name, value, ref_attr, ref_value, expected",
    [
        pytest.param("alpha_i", 0.0, "surface_normal", None, False, id="alpha_i-no-sn"),
        pytest.param(
            "alpha_i", 0.0, "surface_normal", (0, 0, 1), True, id="alpha_i-with-sn"
        ),
        pytest.param(
            "beta_out", 0.0, "surface_normal", None, False, id="beta_out-no-sn"
        ),
        pytest.param(
            "beta_out", 0.0, "surface_normal", (0, 0, 1), True, id="beta_out-with-sn"
        ),
        pytest.param("a_eq_b", True, "surface_normal", None, False, id="a_eq_b-no-sn"),
        pytest.param(
            "a_eq_b", True, "surface_normal", (0, 0, 1), True, id="a_eq_b-with-sn"
        ),
        pytest.param("psi", 0.0, "azimuthal_reference", None, False, id="psi-no-ar"),
        pytest.param(
            "psi", 0.0, "azimuthal_reference", (0, 0, 1), True, id="psi-with-ar"
        ),
        pytest.param("naz", 0.0, "azimuthal_reference", None, False, id="naz-no-ar"),
        pytest.param(
            "naz", 0.0, "azimuthal_reference", (0, 0, 1), True, id="naz-with-ar"
        ),
    ],
)
def test_reference_constraint_has_reference_vector(
    name, value, ref_attr, ref_value, expected
):
    """has_reference_vector() reflects whether the reference vector is set."""
    g = _setup_psic()
    setattr(g, ref_attr, ref_value)
    rc = ReferenceConstraint(name, value)
    assert rc.has_reference_vector(g) is expected


# ---------------------------------------------------------------------------
# Serialisation with reference vectors
# ---------------------------------------------------------------------------


def test_surface_normal_round_trip():
    """surface_normal survives to_dict / from_dict round-trip."""
    import json

    g = _setup_psic()
    g.surface_normal = (0, 0, 1)
    d = g.to_dict()
    assert json.dumps(d)
    assert d["surface_normal"] == [0.0, 0.0, 1.0]
    g2 = AdHocDiffractometer.from_dict(d)
    assert g2.surface_normal == (0.0, 0.0, 1.0)


def test_azimuthal_reference_round_trip():
    """azimuthal_reference survives to_dict / from_dict round-trip."""
    import json

    g = _setup_psic()
    g.azimuthal_reference = (1, 0, 0)
    d = g.to_dict()
    assert json.dumps(d)
    assert d["azimuthal_reference"] == [1.0, 0.0, 0.0]
    g2 = AdHocDiffractometer.from_dict(d)
    assert g2.azimuthal_reference == (1.0, 0.0, 0.0)


def test_surface_normal_none_serialisation():
    """surface_normal=None serialises and restores correctly."""
    g = _setup_psic()
    assert g.surface_normal is None
    d = g.to_dict()
    assert d["surface_normal"] is None
    g2 = AdHocDiffractometer.from_dict(d)
    assert g2.surface_normal is None


# ---------------------------------------------------------------------------
# Reference constraint forward() — still NotImplementedError despite is_implemented=True
# ---------------------------------------------------------------------------


def test_fixed_psi_implemented_with_azref():
    """fixed_psi is_implemented=True when azimuthal_reference is set.

    With azimuthal_reference set, the fixed_psi forward solver is available.
    It acts as a validation filter: it returns bisecting solutions only when
    the natural psi for (h,k,l) matches the stored target.
    """
    g = _setup_psic()
    g.azimuthal_reference = (0, 0, 1)
    g.mode_name = "fixed_psi_vertical"
    cs = g.modes["fixed_psi_vertical"]
    # is_implemented is True — solver available
    assert cs.is_implemented(g) is True
    # has_reference_vector is True — the prerequisite is met
    rc = cs.reference_constraint
    assert rc is not None
    assert rc.has_reference_vector(g) is True
    # forward() with correct psi target returns solutions
    # Natural psi for (1,0,0) with ref=(0,0,1) is 90.0 on psic with identity UB.
    # The mode's default psi target is 0.0, which does NOT match → empty list.
    solutions = g.forward(1, 0, 0)
    assert solutions == []  # natural psi=90 != target psi=0


def test_fixed_psi_not_implemented_without_azref():
    """fixed_psi is_implemented=False when azimuthal_reference is not set."""
    g = _setup_psic()
    assert g.azimuthal_reference is None
    g.mode_name = "fixed_psi_vertical"
    cs = g.modes["fixed_psi_vertical"]
    assert cs.is_implemented(g) is False
    with pytest.raises(NotImplementedError):
        g.forward(1, 0, 0)


# ---------------------------------------------------------------------------
# Issue #175 — surface diffraction forward solvers
# ---------------------------------------------------------------------------

WAVELENGTH = 1.5406


def _setup_surface(factory, surface_normal=(0, 0, 1)):
    """Return a geometry with wavelength, lattice, UB, and surface_normal set."""
    g = factory()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ub_identity(g.sample)
    g.surface_normal = surface_normal
    return g


@pytest.mark.parametrize(
    "factory, mode_name",
    [
        pytest.param(zaxis, "zaxis", id="zaxis-zaxis"),
        pytest.param(zaxis, "reflectivity", id="zaxis-reflectivity"),
        pytest.param(s2d2, "reflectivity", id="s2d2-reflectivity"),
        pytest.param(sixc, "fixed_alpha_zaxis", id="sixc-fixed_alpha_zaxis"),
        pytest.param(sixc, "fixed_beta_zaxis", id="sixc-fixed_beta_zaxis"),
        pytest.param(sixc, "alpha_eq_beta_zaxis", id="sixc-alpha_eq_beta_zaxis"),
    ],
)
def test_surface_mode_is_implemented_with_surface_normal(factory, mode_name):
    """Surface reference modes return is_implemented=True when surface_normal is set."""
    g = _setup_surface(factory)
    assert g.modes[mode_name].is_implemented(g) is True


@pytest.mark.parametrize(
    "factory, mode_name",
    [
        pytest.param(zaxis, "zaxis", id="zaxis-zaxis"),
        pytest.param(zaxis, "reflectivity", id="zaxis-reflectivity"),
        pytest.param(s2d2, "reflectivity", id="s2d2-reflectivity"),
        pytest.param(sixc, "fixed_alpha_zaxis", id="sixc-fixed_alpha_zaxis"),
        pytest.param(sixc, "fixed_beta_zaxis", id="sixc-fixed_beta_zaxis"),
        pytest.param(sixc, "alpha_eq_beta_zaxis", id="sixc-alpha_eq_beta_zaxis"),
    ],
)
def test_surface_mode_not_implemented_without_surface_normal(factory, mode_name):
    """Surface reference modes return is_implemented=False without surface_normal."""
    g = factory()
    assert g.modes[mode_name].is_implemented(g) is False


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l",
    [
        pytest.param(zaxis, "zaxis", 0, 1, 0, id="zaxis-zaxis"),
        pytest.param(zaxis, "reflectivity", 0, 0, 1, id="zaxis-reflectivity"),
        pytest.param(s2d2, "reflectivity", 0, 1, 0, id="s2d2-reflectivity"),
        pytest.param(sixc, "fixed_alpha_zaxis", 0, 1, 0, id="sixc-fixed_alpha_zaxis"),
        pytest.param(sixc, "fixed_beta_zaxis", 0, 1, 0, id="sixc-fixed_beta_zaxis"),
        pytest.param(
            sixc,
            "alpha_eq_beta_zaxis",
            0,
            1,
            0,
            id="sixc-alpha_eq_beta_zaxis",
        ),
    ],
)
def test_surface_mode_returns_solutions(factory, mode_name, h, k, l):  # noqa: E741
    """Surface modes return at least one solution when surface_normal is set."""
    g = _setup_surface(factory)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l",
    [
        pytest.param(zaxis, "zaxis", 0, 1, 0, id="zaxis-zaxis-alpha_i=0"),
        pytest.param(
            sixc,
            "fixed_alpha_zaxis",
            0,
            1,
            0,
            id="sixc-fixed_alpha-alpha_i=0",
        ),
    ],
)
def test_surface_alpha_i_fixed_constraint_satisfied(factory, mode_name, h, k, l):  # noqa: E741
    """alpha_i modes: incidence angle equals declared target (0°) in all solutions."""
    g = _setup_surface(factory)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        ai = incidence_angle(g, angles=sol)
        assert ai == pytest.approx(0.0, abs=1e-4)


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l",
    [
        pytest.param(
            sixc,
            "fixed_beta_zaxis",
            0,
            1,
            0,
            id="sixc-fixed_beta-beta_out=0",
        ),
    ],
)
def test_surface_beta_out_fixed_constraint_satisfied(factory, mode_name, h, k, l):  # noqa: E741
    """beta_out modes: exit angle equals declared target (0°) in all solutions."""
    g = _setup_surface(factory)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        bo = exit_angle(g, angles=sol)
        assert bo == pytest.approx(0.0, abs=1e-4)


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l",
    [
        pytest.param(zaxis, "reflectivity", 0, 0, 1, id="zaxis-reflectivity"),
        pytest.param(s2d2, "reflectivity", 0, 1, 0, id="s2d2-reflectivity"),
        pytest.param(sixc, "alpha_eq_beta_zaxis", 0, 1, 0, id="sixc-alpha_eq_beta"),
    ],
)
def test_surface_a_eq_b_constraint_satisfied(factory, mode_name, h, k, l):  # noqa: E741
    """a_eq_b modes: alpha_i ≈ beta_out in all solutions."""
    g = _setup_surface(factory)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        ai = incidence_angle(g, angles=sol)
        bo = exit_angle(g, angles=sol)
        assert ai == pytest.approx(bo, abs=1e-4)


def test_surface_mode_not_implemented_raises():
    """Surface modes without surface_normal raise NotImplementedError on forward()."""
    g = zaxis()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=4.0)
    ub_identity(g.sample)
    # No surface_normal set
    g.mode_name = "zaxis"
    with pytest.raises(NotImplementedError):
        g.forward(0, 1, 0)


# ---------------------------------------------------------------------------
# omega_pseudo (SPEC OMEGA = Q[6], issue #264)
# ---------------------------------------------------------------------------


from contextlib import nullcontext as does_not_raise  # noqa: E402

from ad_hoc_diffractometer.reference import omega_pseudo  # noqa: E402


def test_omega_pseudo_requires_wavelength():
    """omega_pseudo raises ValueError when wavelength is None."""
    g = psic()
    # Do NOT set wavelength
    with pytest.raises(ValueError, match=re.escape("wavelength")):
        omega_pseudo(g)


def test_omega_pseudo_requires_chi_stage():
    """omega_pseudo raises KeyError on a geometry with no 'chi' stage."""
    g = zaxis()
    g.wavelength = WAVELENGTH
    with pytest.raises(KeyError, match=re.escape("no sample stage named 'chi'")):
        omega_pseudo(g)


def test_omega_pseudo_does_not_require_surface_normal():
    """omega_pseudo works with no surface_normal or azimuthal_reference."""
    g = _setup_psic()
    assert g.surface_normal is None
    assert g.azimuthal_reference is None
    g.mode_name = "bisecting_vertical"
    sols = g.forward(1, 0, 0)
    for s in sols:
        om = omega_pseudo(g, angles=s)
        assert isinstance(om, float)


def test_omega_pseudo_uses_current_angles_when_none():
    """omega_pseudo uses current stage angles when angles=None."""
    g = _setup_psic()
    om = omega_pseudo(g, angles=None)
    assert isinstance(om, float)


def test_omega_pseudo_zero_at_bisecting_vertical():
    """At bisecting_vertical (mu=nu=0, eta=delta/2), OMEGA = 0."""
    g = _setup_psic()
    g.mode_name = "bisecting_vertical"
    sols = g.forward(1, 0, 0)
    assert len(sols) > 0
    for s in sols:
        om = omega_pseudo(g, angles=s)
        assert om == pytest.approx(0.0, abs=1e-6), (
            f"OMEGA should be 0 at bisecting; got {om} for {s}"
        )


def test_omega_pseudo_zero_at_bisecting_horizontal():
    """At bisecting_horizontal (eta=delta=0, mu=nu/2), OMEGA = 0."""
    g = _setup_psic()
    g.mode_name = "bisecting_horizontal"
    sols = g.forward(0, 0, 1)
    assert len(sols) > 0
    for s in sols:
        om = omega_pseudo(g, angles=s)
        assert om == pytest.approx(0.0, abs=1e-6), (
            f"OMEGA should be 0 at bisecting_horizontal; got {om} for {s}"
        )


def test_omega_pseudo_independent_of_phi():
    """OMEGA depends only on the outer sample stages and the detector;
    it is independent of phi."""
    g = _setup_psic()
    angles_a = {
        "mu": 5.0,
        "eta": 12.0,
        "chi": 30.0,
        "phi": 0.0,
        "nu": 0.0,
        "delta": 24.0,
    }
    angles_b = dict(angles_a)
    angles_b["phi"] = 73.0
    om_a = omega_pseudo(g, angles=angles_a)
    om_b = omega_pseudo(g, angles=angles_b)
    assert om_a == pytest.approx(om_b, abs=1e-9), (
        f"OMEGA must be independent of phi; got {om_a} vs {om_b}"
    )


def test_omega_pseudo_independent_of_chi():
    """OMEGA depends only on the outer sample stages and the detector;
    it is independent of chi (chi rotates Q and the chi-circle plane
    together, preserving their relative angle)."""
    g = _setup_psic()
    angles_a = {
        "mu": 5.0,
        "eta": 12.0,
        "chi": 30.0,
        "phi": 17.0,
        "nu": 0.0,
        "delta": 24.0,
    }
    angles_b = dict(angles_a)
    angles_b["chi"] = 91.0
    om_a = omega_pseudo(g, angles=angles_a)
    om_b = omega_pseudo(g, angles=angles_b)
    assert om_a == pytest.approx(om_b, abs=1e-9)


@pytest.mark.parametrize(
    "name, expected",
    [
        pytest.param("psi", True, id="psi"),
        pytest.param("alpha_i", True, id="alpha_i"),
        pytest.param("beta_out", True, id="beta_out"),
        pytest.param("a_eq_b", True, id="a_eq_b"),
        pytest.param("naz", True, id="naz"),
        pytest.param("omega", True, id="omega"),
        pytest.param("not_a_pseudo_angle", False, id="invalid"),
    ],
)
def test_reference_constraint_accepts_omega(name, expected):
    """ReferenceConstraint('omega', value) is now a valid constraint."""
    context = (
        does_not_raise()
        if expected
        else pytest.raises(ValueError, match=re.escape("ReferenceConstraint name"))
    )
    with context:
        if name == "a_eq_b":
            ReferenceConstraint(name, True)
        else:
            ReferenceConstraint(name, 0.0)


def test_reference_constraint_omega_is_implemented_on_psic():
    """ReferenceConstraint('omega', 0).is_implemented(psic) is True."""
    g = _setup_psic()
    rc = ReferenceConstraint("omega", 0.0)
    assert rc.is_implemented(g) is True
    assert rc.has_reference_vector(g) is True


def test_reference_constraint_omega_not_implemented_on_zaxis():
    """ReferenceConstraint('omega', 0).is_implemented(zaxis) is False — no chi."""
    g = zaxis()
    g.wavelength = WAVELENGTH
    rc = ReferenceConstraint("omega", 0.0)
    assert rc.is_implemented(g) is False
    # has_reference_vector still True (omega needs no reference vector)
    assert rc.has_reference_vector(g) is True


def test_reference_constraint_omega_serialization_round_trip():
    """ReferenceConstraint('omega', value) round-trips through to_dict / from_dict."""
    rc = ReferenceConstraint("omega", 12.5)
    d = rc.to_dict()
    assert d["type"] == "ReferenceConstraint"
    assert d["name"] == "omega"
    assert d["value"] == 12.5
    rc2 = ReferenceConstraint.from_dict(d)
    assert rc == rc2
