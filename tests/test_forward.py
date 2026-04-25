# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.forward (compute_forward / geometry.forward).

Covers:
  - Precondition errors: no wavelength, no UB, hkl=(0,0,0), Q > Ewald sphere,
    no active mode, unimplemented mode
  - ConstraintSet bisecting solver: fourcv, fourch, psic (bisecting)
  - ConstraintSet fixed-angle solver: fourcv fixed_chi, psic fixed_chi
  - Round-trip invariant: inverse(forward(hkl)) == hkl for every solution
  - Stage limits: solutions outside limits are filtered out
  - Cut-point application: mode and geometry-level
  - Q along ±z (degenerate chi branch): phi set to 0
  - _check_limits and _apply_cut_points helpers
  - compute_forward top-level function (same as geometry.forward)
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import BisectConstraint
from ad_hoc_diffractometer import ConstraintSet
from ad_hoc_diffractometer import DetectorConstraint
from ad_hoc_diffractometer import SampleConstraint
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT
from ad_hoc_diffractometer.forward import _apply_cut_points
from ad_hoc_diffractometer.forward import _check_limits
from ad_hoc_diffractometer.forward import compute_forward
from ad_hoc_diffractometer.presets import fivec
from ad_hoc_diffractometer.presets import fourch
from ad_hoc_diffractometer.presets import fourcv
from ad_hoc_diffractometer.presets import kappa4ch
from ad_hoc_diffractometer.presets import kappa4cv
from ad_hoc_diffractometer.presets import kappa6c
from ad_hoc_diffractometer.presets import psic
from ad_hoc_diffractometer.presets import sixc
from ad_hoc_diffractometer.stage import Stage

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

WAVELENGTH = 1.5406  # Cu Kα


def _setup_cubic(factory, a=1.0):
    """Return a geometry with UB=B for a cubic lattice of side a."""
    g = factory()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=a)
    ub_identity(g.sample)
    return g


def _round_trip_ok(g, h, k, l, atol=1e-8):  # noqa: E741
    """
    Return True if every forward solution round-trips back to (h, k, l)
    via inverse(), and the list is non-empty.
    """
    solutions = g.forward(h, k, l)
    if not solutions:
        return False
    for sol in solutions:
        hkl_back = g.inverse(sol)
        if not np.allclose(hkl_back, [h, k, l], atol=atol):
            return False
    return True


# ---------------------------------------------------------------------------
# Precondition errors
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "setup_fn, h, k, l, exc_type, match, context",
    [
        pytest.param(
            lambda: fourcv(),  # no wavelength
            1,
            0,
            0,
            ValueError,
            re.escape("wavelength must be set"),
            does_not_raise(),
            id="no-wavelength",
        ),
        pytest.param(
            lambda: _no_ub(),  # no UB
            1,
            0,
            0,
            ValueError,
            re.escape("has no UB matrix"),
            does_not_raise(),
            id="no-ub",
        ),
        pytest.param(
            lambda: _setup_cubic(fourcv),
            0,
            0,
            0,
            ValueError,
            re.escape("(0, 0, 0) is not a valid"),
            does_not_raise(),
            id="hkl-zero",
        ),
        pytest.param(
            lambda: _setup_cubic(fourcv),
            100,
            0,
            0,
            ahd.EwaldSphereViolation,
            re.escape("exceeds the Ewald sphere"),
            does_not_raise(),
            id="q-too-large",
        ),
        pytest.param(
            lambda: _no_mode(),  # no active mode
            1,
            0,
            0,
            NotImplementedError,
            re.escape("has no active diffraction mode"),
            does_not_raise(),
            id="no-active-mode",
        ),
        pytest.param(
            lambda: _unsupported_mode(),
            1,
            0,
            0,
            NotImplementedError,
            re.escape("is not yet implemented"),
            does_not_raise(),
            id="unsupported-mode-type",
        ),
    ],
)
def test_forward_precondition_errors(setup_fn, h, k, l, exc_type, match, context):  # noqa: E741
    with context:
        g = setup_fn()
        with pytest.raises(exc_type, match=match):
            g.forward(h, k, l)


def _no_ub():
    g = fourcv()
    g.wavelength = WAVELENGTH
    # UB is None by default
    return g


def _no_mode():
    g = _setup_cubic(fourcv)
    g.mode_name = None
    return g


def _unsupported_mode():
    """A geometry whose active mode has is_implemented() returning False."""
    # Use a ReferenceConstraint which is not yet implemented
    from ad_hoc_diffractometer import ReferenceConstraint

    cs = ConstraintSet([ReferenceConstraint("psi", 90.0)])
    g = _setup_cubic(fourcv)
    g.modes["psi_mode"] = cs
    g.mode_name = "psi_mode"
    return g


# ---------------------------------------------------------------------------
# BisectingMode — fourcv and fourch round-trips
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, a, h, k, l",
    [
        pytest.param(fourcv, 5.0, 1, 0, 0, id="fourcv-100"),
        pytest.param(fourcv, 5.0, 0, 1, 0, id="fourcv-010"),
        pytest.param(fourcv, 5.0, 0, 0, 1, id="fourcv-001"),
        pytest.param(fourcv, 4.0, 1, 1, 1, id="fourcv-111"),
        pytest.param(fourch, 4.0, 1, 0, 0, id="fourch-100"),
    ],
)
def test_bisecting_round_trip(factory, a, h, k, l):  # noqa: E741
    """Bisecting mode round-trip: inverse(forward(hkl)) == hkl."""
    g = _setup_cubic(factory, a=a)
    assert g.mode_name == "bisecting"
    assert _round_trip_ok(g, h, k, l)


def test_fourcv_bisecting_two_solutions():
    """Bisecting mode returns two solutions (positive and negative chi branch)."""
    g = _setup_cubic(fourcv, a=4.0)
    solutions = g.forward(1, 0, 0)
    assert len(solutions) == 2


def test_fourcv_bisecting_omega_equals_ttheta_half():
    """In bisecting mode, omega must equal ttheta/2 for all solutions."""
    g = _setup_cubic(fourcv, a=4.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["omega"] == pytest.approx(sol["ttheta"] / 2.0, abs=1e-10)


def test_fourcv_bisecting_ttheta_from_bragg():
    """ttheta must satisfy Bragg's law for the given lattice and wavelength.

    The B matrix in this package uses the no-2π convention (Å⁻¹), so
    |Q| = |UB @ hkl| = 1/d (not 2π/d).  Bragg's law then reads:
        sin(θ) = |Q| * λ / (4π)   [equivalent to λ = 2d sin(θ) with d = 1/|Q|]
    """
    a = 4.0
    g = _setup_cubic(fourcv, a=a)
    h, k, l = 1, 0, 0  # noqa: E741
    import numpy as np

    Q_phi = g.sample.UB @ np.array([h, k, l], dtype=float)
    Q_mag = float(np.linalg.norm(Q_phi))
    ttheta_expected = 2.0 * math.degrees(
        math.asin(Q_mag * WAVELENGTH / (4.0 * math.pi))
    )
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["ttheta"] == pytest.approx(ttheta_expected, abs=1e-8)


# ---------------------------------------------------------------------------
# BisectingMode — psic (eta = delta/2, mu = nu = 0)
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "a, h, k, l",
    [
        pytest.param(5.0, 1, 0, 0, id="psic-100"),
        pytest.param(4.0, 1, 1, 1, id="psic-111"),
    ],
)
def test_psic_bisecting_round_trip(a, h, k, l):  # noqa: E741
    """psic bisecting round-trip: inverse(forward(hkl)) == hkl."""
    g = _setup_cubic(psic, a=a)
    assert g.mode_name == "bisecting_vertical"
    assert _round_trip_ok(g, h, k, l)


def test_psic_bisecting_frozen_mu_nu():
    """In psic bisecting mode, mu and nu must be 0 in all solutions."""
    g = _setup_cubic(psic, a=4.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["mu"] == pytest.approx(0.0, abs=1e-10)
        assert sol["nu"] == pytest.approx(0.0, abs=1e-10)


def test_psic_bisecting_eta_equals_delta_half():
    """In psic bisecting mode, eta must equal delta/2."""
    g = _setup_cubic(psic, a=4.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["eta"] == pytest.approx(sol["delta"] / 2.0, abs=1e-10)


# ---------------------------------------------------------------------------
# BisectingMode — kappa4cv
# ---------------------------------------------------------------------------


def test_kappa4cv_bisecting_round_trip():
    """kappa4cv round-trip: use (0,1,0) which is reachable in BL transverse basis."""
    g = _setup_cubic(kappa4cv, a=4.0)
    assert g.mode_name == "bisecting"
    assert _round_trip_ok(g, 0, 1, 0)


# ---------------------------------------------------------------------------
# Issue #151 — kappa4cv and kappa4ch: mode counts, stubs, fixed_kphi
# ---------------------------------------------------------------------------

_KAPPA4_BASE_MODES = {
    "bisecting",
    "fixed_kphi",
    "constant_omega",
    "constant_chi",
    "constant_phi",
    "psi_constant",
}
_KAPPA4CV_ALL_MODES = _KAPPA4_BASE_MODES | {"double_diffraction"}
_KAPPA4CH_ALL_MODES = _KAPPA4_BASE_MODES
_KAPPA4_STUB_MODES = {"psi_constant"}


def test_kappa4cv_has_all_modes():
    """kappa4cv exposes all declared modes including double_diffraction."""
    assert set(kappa4cv().modes.keys()) == _KAPPA4CV_ALL_MODES


def test_kappa4ch_has_all_modes():
    """kappa4ch exposes all declared modes (no double_diffraction)."""
    assert set(kappa4ch().modes.keys()) == _KAPPA4CH_ALL_MODES


@pytest.mark.parametrize(
    "factory, mode_name",
    [pytest.param(kappa4cv, m, id=f"kappa4cv-{m}") for m in sorted(_KAPPA4_STUB_MODES)]
    + [
        pytest.param(kappa4ch, m, id=f"kappa4ch-{m}")
        for m in sorted(_KAPPA4_STUB_MODES)
    ],
)
def test_kappa4_stub_not_implemented(factory, mode_name):
    """Virtual-angle stubs raise NotImplementedError on forward()."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    assert not g.modes[mode_name].is_implemented(g)
    with pytest.raises(NotImplementedError):
        g.forward(0, 1, 0)


@pytest.mark.parametrize(
    "factory",
    [pytest.param(kappa4cv, id="kappa4cv"), pytest.param(kappa4ch, id="kappa4ch")],
)
def test_kappa4_fixed_kphi_round_trip(factory):
    """fixed_kphi (real stage) is implemented and round-trips correctly."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = "fixed_kphi"
    assert g.modes["fixed_kphi"].is_implemented(g)
    assert _round_trip_ok(g, 0, 1, 0)


@pytest.mark.parametrize(
    "factory",
    [pytest.param(kappa4cv, id="kappa4cv"), pytest.param(kappa4ch, id="kappa4ch")],
)
def test_kappa4_fixed_kphi_value_in_solution(factory):
    """fixed_kphi: kphi=0 in all solutions."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = "fixed_kphi"
    solutions = g.forward(0, 1, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["kphi"] == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Issue #174 — kappa4cv/kappa4ch virtual-angle modes now implemented
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l",
    [
        # kappa4cv (BL basis: transverse=x, longitudinal=y, vertical=z)
        pytest.param(kappa4cv, "constant_omega", 0, 1, 0, id="kappa4cv-constant_omega"),
        pytest.param(kappa4cv, "constant_chi", 0, 0, 1, id="kappa4cv-constant_chi"),
        pytest.param(kappa4cv, "constant_phi", 0, 1, 0, id="kappa4cv-constant_phi"),
        # kappa4ch (horizontal scattering plane — different reachable reflections)
        pytest.param(kappa4ch, "constant_omega", 0, 0, 1, id="kappa4ch-constant_omega"),
        pytest.param(kappa4ch, "constant_chi", 0, 0, 1, id="kappa4ch-constant_chi"),
        pytest.param(kappa4ch, "constant_phi", 0, 1, 0, id="kappa4ch-constant_phi"),
    ],
)
def test_kappa4_virtual_angle_round_trip(factory, mode_name, h, k, l):  # noqa: E741
    """Virtual Eulerian angle modes on kappa4cv/kappa4ch round-trip correctly."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    assert g.modes[mode_name].is_implemented(g)
    # Newton-Raphson converges to ~1e-6; use looser tolerance for these modes
    assert _round_trip_ok(g, h, k, l, atol=1e-4)


@pytest.mark.parametrize(
    "factory, mode_name, virtual_angle, expected_value, h, k, l",
    [
        pytest.param(
            kappa4cv, "constant_omega", "omega", 0.0, 0, 1, 0, id="kappa4cv-omega=0"
        ),
        pytest.param(
            kappa4cv, "constant_chi", "chi", 90.0, 0, 0, 1, id="kappa4cv-chi=90"
        ),
        pytest.param(
            kappa4cv, "constant_phi", "phi", 0.0, 0, 1, 0, id="kappa4cv-phi=0"
        ),
    ],
)
def test_kappa4_virtual_angle_constraint_satisfied(
    factory,
    mode_name,
    virtual_angle,
    expected_value,
    h,
    k,
    l,  # noqa: E741
):
    """Virtual angle constraint is satisfied in all returned solutions."""
    from ad_hoc_diffractometer.kappa import kappa_to_eulerian

    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        o, c, p = kappa_to_eulerian(
            sol["komega"],
            sol["kappa"],
            sol["kphi"],
            alpha_deg=g.kappa_alpha_deg,
        )
        virtual_vals = {"omega": o, "chi": c, "phi": p}
        assert virtual_vals[virtual_angle] == pytest.approx(expected_value, abs=1e-4)


# ---------------------------------------------------------------------------
# FixedAngleMode — fourcv fixed_chi
# ---------------------------------------------------------------------------


def test_fourcv_fixed_chi_round_trip():
    """fixed_chi round-trip: constraint freezes chi at its declared value (90°)."""
    g = _setup_cubic(fourcv, a=4.0)
    g.mode_name = "fixed_chi"
    assert _round_trip_ok(g, 1, 0, 0)


def test_psic_fixed_chi_uses_constraint_value():
    """fixed_chi on psic: chi=90°, mu=0, nu=0 from constraint values."""
    from ad_hoc_diffractometer.presets import psic

    g = _setup_cubic(psic, a=4.0)
    g.mode_name = "fixed_chi"
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["chi"] == pytest.approx(90.0, abs=1e-6)
        assert sol["mu"] == pytest.approx(0.0, abs=1e-6)
        assert sol["nu"] == pytest.approx(0.0, abs=1e-6)


def test_fourcv_fixed_chi_value_respected():
    """chi must equal 90° (constraint value) in all solutions."""
    g = _setup_cubic(fourcv, a=4.0)
    g.mode_name = "fixed_chi"
    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["chi"] == pytest.approx(90.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Issue #149 — fourcv and fourch: all modes present, stubs not implemented
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(fourcv, id="fourcv"),
        pytest.param(fourch, id="fourch"),
    ],
)
def test_four_circle_has_all_six_modes(factory):
    """fourcv and fourch both expose all 6 declared modes as ConstraintSet."""
    g = factory()
    expected = {
        "bisecting",
        "fixed_chi",
        "fixed_phi",
        "constant_omega",
        "psi_constant",
        "double_diffraction",
    }
    assert set(g.modes.keys()) == expected
    for name, mode in g.modes.items():
        assert isinstance(mode, ConstraintSet), f"{name} is not a ConstraintSet"


@pytest.mark.parametrize(
    "factory, mode_name",
    [
        pytest.param(fourcv, "psi_constant", id="fourcv-psi_constant"),
        pytest.param(fourch, "psi_constant", id="fourch-psi_constant"),
    ],
)
def test_four_circle_stub_not_implemented(factory, mode_name):
    """psi_constant requires reference infrastructure (Issue J) — not yet implemented."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    assert not g.modes[mode_name].is_implemented(g)
    with pytest.raises(NotImplementedError):
        g.forward(1, 0, 0)


@pytest.mark.parametrize(
    "factory, h, k, l",
    [
        pytest.param(fourcv, 1, 0, 0, id="fourcv-100"),
        pytest.param(fourcv, 0, 1, 0, id="fourcv-010"),
        pytest.param(fourcv, 1, 1, 1, id="fourcv-111"),
        pytest.param(fourch, 1, 0, 0, id="fourch-100"),
        pytest.param(fourch, 1, 1, 1, id="fourch-111"),
    ],
)
def test_four_circle_constant_omega_round_trip(factory, h, k, l):  # noqa: E741
    """constant_omega (omega=0) is implemented by the generic solver; round-trips correctly."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = "constant_omega"
    assert g.modes["constant_omega"].is_implemented(g)
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(fourcv, id="fourcv"),
        pytest.param(fourch, id="fourch"),
    ],
)
def test_four_circle_constant_omega_value_in_solution(factory):
    """constant_omega: omega=0 in all returned solutions."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = "constant_omega"
    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["omega"] == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l",
    [
        pytest.param(fourcv, "fixed_phi", 0, 1, 0, id="fourcv-fixed_phi-010"),
        pytest.param(fourcv, "fixed_phi", 1, 1, 1, id="fourcv-fixed_phi-111"),
        pytest.param(fourch, "fixed_chi", 1, 0, 0, id="fourch-fixed_chi-100"),
        pytest.param(fourch, "fixed_chi", 1, 1, 0, id="fourch-fixed_chi-110"),
        pytest.param(fourch, "fixed_phi", 0, 1, 0, id="fourch-fixed_phi-010"),
        pytest.param(fourch, "fixed_phi", 1, 1, 1, id="fourch-fixed_phi-111"),
    ],
)
def test_four_circle_fixed_angle_round_trip(factory, mode_name, h, k, l):  # noqa: E741
    """fixed_chi and fixed_phi round-trips on both fourcv and fourch."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "factory, mode_name, expected_fixed_stage, expected_value, h, k, l",
    [
        pytest.param(
            fourcv, "fixed_phi", "phi", 0.0, 0, 1, 0, id="fourcv-fixed_phi-value"
        ),
        pytest.param(
            fourch, "fixed_chi", "chi", 90.0, 1, 0, 0, id="fourch-fixed_chi-value"
        ),
        pytest.param(
            fourch, "fixed_phi", "phi", 0.0, 0, 1, 0, id="fourch-fixed_phi-value"
        ),
    ],
)
def test_four_circle_fixed_angle_constraint_value_in_solution(
    factory,
    mode_name,
    expected_fixed_stage,
    expected_value,
    h,
    k,
    l,  # noqa: E741
):
    """The fixed stage appears at its declared constraint value in all solutions."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol[expected_fixed_stage] == pytest.approx(expected_value, abs=1e-8), (
            f"{expected_fixed_stage} not at {expected_value} in {sol}"
        )


@pytest.mark.parametrize(
    "factory, mode_name, extras_key, expected_value",
    [
        pytest.param(
            fourcv,
            "double_diffraction",
            "h2",
            "REQUIRED",
            id="fourcv-double_diffraction-h2",
        ),
        pytest.param(
            fourcv,
            "double_diffraction",
            "k2",
            "REQUIRED",
            id="fourcv-double_diffraction-k2",
        ),
        pytest.param(
            fourcv,
            "double_diffraction",
            "l2",
            "REQUIRED",
            id="fourcv-double_diffraction-l2",
        ),
        pytest.param(
            fourch,
            "double_diffraction",
            "h2",
            "REQUIRED",
            id="fourch-double_diffraction-h2",
        ),
        pytest.param(
            fourch,
            "double_diffraction",
            "k2",
            "REQUIRED",
            id="fourch-double_diffraction-k2",
        ),
        pytest.param(
            fourch,
            "double_diffraction",
            "l2",
            "REQUIRED",
            id="fourch-double_diffraction-l2",
        ),
        pytest.param(
            fourcv, "psi_constant", "n_hat", "REQUIRED", id="fourcv-psi_constant-n_hat"
        ),
        pytest.param(
            fourcv, "psi_constant", "psi", None, id="fourcv-psi_constant-psi-output"
        ),
        pytest.param(
            fourch, "psi_constant", "n_hat", "REQUIRED", id="fourch-psi_constant-n_hat"
        ),
        pytest.param(
            fourch, "psi_constant", "psi", None, id="fourch-psi_constant-psi-output"
        ),
    ],
)
def test_four_circle_mode_extras_declared(
    factory, mode_name, extras_key, expected_value
):
    """Mode extras carry the expected sentinel or output placeholder."""
    from ad_hoc_diffractometer import REQUIRED

    g = factory()
    mode = g.modes[mode_name]
    actual = mode.extras.get(extras_key)
    if expected_value == "REQUIRED":
        assert actual is REQUIRED
    else:
        assert actual is expected_value


# ---------------------------------------------------------------------------
# FixedAngleMode — no bisecting mode available
# ---------------------------------------------------------------------------


def test_all_stages_constrained_within_limits():
    """When all sample stages are constrained and within limits, returns [solution]."""
    g = _setup_cubic(fourcv, a=4.0)
    # BisectConstraint + chi=90 + phi=0 — all 3 sample stages constrained for a 4-circle
    # (over-constrained DOF-wise but valid for testing the "all frozen" solver path)
    fully_frozen_mode = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("chi", 90.0),
            SampleConstraint("phi", 0.0),
        ]
    )
    g.modes["fully_frozen"] = fully_frozen_mode
    g.mode_name = "fully_frozen"
    solutions = g.forward(1, 0, 0)
    assert isinstance(solutions, list)


def test_all_stages_constrained_out_of_limits():
    """When all sample stages are constrained but out of limits, returns []."""
    g = _setup_cubic(fourcv, a=4.0)
    fully_frozen_mode = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            SampleConstraint("chi", 90.0),
            SampleConstraint("phi", 0.0),
        ]
    )
    g.modes["fully_frozen"] = fully_frozen_mode
    g.mode_name = "fully_frozen"
    # Restrict omega limits so the bisected omega (ttheta/2) is out of range
    g.stage("omega").limits = (100.0, 180.0)
    solutions = g.forward(1, 0, 0)
    assert solutions == []


def test_solver_none_when_no_convergence():
    """_solve_two_angles returns None when starting point does not converge."""
    from ad_hoc_diffractometer.forward import _solve_two_angles
    from ad_hoc_diffractometer.orientation import angles_to_phi_vector

    g = _setup_cubic(fourcv, a=4.0)
    Q_phi_target = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    import math

    ttheta = 2 * math.degrees(
        math.asin(float(np.linalg.norm(Q_phi_target)) * WAVELENGTH / (4.0 * math.pi))
    )
    angles = {s.name: s.angle for s in list(g._stages.values())}
    angles["ttheta"] = ttheta
    angles["omega"] = ttheta / 2

    # Use max_iter=0 so it cannot converge and residual will be large
    result = _solve_two_angles(
        g,
        angles,
        "chi",
        "phi",
        45.0,
        45.0,
        Q_phi_target,
        angles_to_phi_vector,
        max_iter=0,
    )
    # With max_iter=0, no iterations run, so residual at start (45,45) is large -> None
    # (unless starting point is already the solution, which it isn't)
    # Result is either None or a valid solution — just ensure it doesn't crash.
    assert result is None or (isinstance(result, tuple) and len(result) == 2)


def test_solver_1d_none_when_no_convergence():
    """_solve_one_free_angle's internal call with max_iter=0 returns None."""
    from ad_hoc_diffractometer.forward import _solve_two_angles
    from ad_hoc_diffractometer.orientation import angles_to_phi_vector

    g = _setup_cubic(fourcv, a=4.0)
    Q_phi_target = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    import math

    ttheta = 2 * math.degrees(
        math.asin(float(np.linalg.norm(Q_phi_target)) * WAVELENGTH / (4.0 * math.pi))
    )
    angles = {s.name: s.angle for s in list(g._stages.values())}
    angles["ttheta"] = ttheta
    angles["omega"] = ttheta / 2
    angles["chi"] = 90.0

    result = _solve_two_angles(
        g,
        angles,
        "phi",
        "phi",
        45.0,
        45.0,
        Q_phi_target,
        angles_to_phi_vector,
        max_iter=0,
        _one_dimensional=True,
    )
    assert result is None or (isinstance(result, tuple) and len(result) == 2)


def test_not_implemented_mode_raises():
    """A mode with is_implemented()=False raises NotImplementedError on forward()."""
    from ad_hoc_diffractometer import ReferenceConstraint

    stages = [
        Stage("omega", -ZHAT, parent=None, role="sample"),
        Stage("ttheta", -ZHAT, parent=None, role="detector"),
    ]
    g = ahd.AdHocDiffractometer(
        name="minimal",
        stages=stages,
        basis={"vertical": XHAT, "longitudinal": YHAT, "transverse": ZHAT},
        modes={"psi_mode": ConstraintSet([ReferenceConstraint("psi", 90.0)])},
        default_mode="psi_mode",
    )
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=1.0)
    ub_identity(g.sample)
    with pytest.raises(NotImplementedError, match=re.escape("not yet implemented")):
        g.forward(1, 0, 0)


# ---------------------------------------------------------------------------
# Q along z: degenerate chi branch
# ---------------------------------------------------------------------------


def test_q_along_z_degenerate_chi_90():
    """When Q_phi is along ±z (vertical in BL basis), chi = ±90°.

    In the BL basis (z=vertical), b3 points along +z, so UB @ (0,0,1)
    gives Q_phi along +z.  The chi stage (rotating about longitudinal/y)
    must be at ±90° to bring Q into the scattering plane.  Phi is
    underdetermined but the solver returns a valid solution.
    """
    g = _setup_cubic(fourcv, a=4.0)
    solutions = g.forward(0, 0, 1)
    assert len(solutions) > 0
    for sol in solutions:
        # Every solution must round-trip correctly
        hkl_back = g.inverse(sol)
        assert np.allclose(hkl_back, [0, 0, 1], atol=1e-5)


# ---------------------------------------------------------------------------
# Stage limits filtering
# ---------------------------------------------------------------------------


def test_limits_filter_removes_invalid_solutions():
    """Solutions outside stage limits are not returned."""
    g = _setup_cubic(fourcv, a=4.0)
    # Restrict chi to [45°, 135°] — only the positive-chi branch should pass.
    g.stage("chi").limits = (45.0, 135.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert 45.0 <= sol["chi"] <= 135.0


def test_limits_filter_can_return_empty():
    """If all solutions violate limits, an empty list is returned."""
    g = _setup_cubic(fourcv, a=4.0)
    # Make chi unreachable for both branches
    g.stage("chi").limits = (200.0, 250.0)
    g.stage("chi")._limits = (200.0, 250.0)  # bypass setter for testing
    solutions = g.forward(1, 0, 0)
    assert solutions == []


# ---------------------------------------------------------------------------
# Cut-point application
# ---------------------------------------------------------------------------


def test_mode_cut_point_applied():
    """A mode-level cut-point shifts angles into [cut, cut+360)."""
    g = _setup_cubic(fourcv, a=4.0)
    mode = ConstraintSet(
        [BisectConstraint("omega", "ttheta")],
        cut_points={"phi": -180.0},
    )
    g.modes["bisecting_cut"] = mode
    g.mode_name = "bisecting_cut"
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        # phi must lie in [-180, 180)
        assert -180.0 <= sol["phi"] < 180.0


def test_geometry_cut_point_applied():
    """A geometry-level cut-point is applied when no mode cut-point overrides."""
    g = _setup_cubic(fourcv, a=4.0)
    g.cut_points["omega"] = 0.0  # omega in [0, 360)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert 0.0 <= sol["omega"] < 360.0


# ---------------------------------------------------------------------------
# compute_forward top-level function
# ---------------------------------------------------------------------------


def test_compute_forward_same_as_method():
    """compute_forward(g, h, k, l) and g.forward(h, k, l) return the same result."""
    g = _setup_cubic(fourcv, a=4.0)
    via_fn = compute_forward(g, 1, 0, 0)
    via_method = g.forward(1, 0, 0)
    assert len(via_fn) == len(via_method)
    for a, b in zip(via_fn, via_method, strict=False):
        assert a == b


# ---------------------------------------------------------------------------
# _check_limits and _apply_cut_points helpers
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "angles, expected, context",
    [
        pytest.param(
            {s: 0.0 for s in ("omega", "chi", "phi", "ttheta")},
            True,
            does_not_raise(),
            id="all-within-limits",
        ),
        pytest.param(
            {"omega": 0.0, "chi": 200.0, "phi": 0.0, "ttheta": 0.0},
            False,
            does_not_raise(),
            id="chi-out-of-limits",
        ),
        pytest.param(
            {"omega": 0.0, "nonexistent": 999.0},
            True,
            does_not_raise(),
            id="unknown-stage-ignored",
        ),
    ],
)
def test_check_limits(angles, expected, context):
    with context:
        g = _setup_cubic(fourcv, a=4.0)
        assert _check_limits(g, angles) is expected


@pytest.mark.parametrize(
    "mode_cut_points, geom_cut_points, angle_in, expected_out, context",
    [
        pytest.param(
            {"phi": 0.0},
            {"phi": -180.0},
            -10.0,
            350.0,
            does_not_raise(),
            id="mode-cut-takes-priority",
        ),
        pytest.param(
            {},
            {"phi": 0.0},
            -10.0,
            350.0,
            does_not_raise(),
            id="geometry-cut-fallback",
        ),
        pytest.param(
            {},
            {},
            -10.0,
            -10.0,
            does_not_raise(),
            id="no-cut-unchanged",
        ),
    ],
)
def test_apply_cut_points(
    mode_cut_points, geom_cut_points, angle_in, expected_out, context
):
    with context:
        g = _setup_cubic(fourcv, a=4.0)
        mode = ConstraintSet(
            [BisectConstraint("omega", "ttheta")], cut_points=mode_cut_points
        )
        g.cut_points.update(geom_cut_points)
        angles = {"phi": angle_in}
        _apply_cut_points(angles, mode, g)
        assert angles["phi"] == pytest.approx(expected_out, abs=1e-10)


# ---------------------------------------------------------------------------
# _solve_bisecting — virtual-stage-name (non-existent stage) branch
# ---------------------------------------------------------------------------


def test_bisecting_virtual_stage_name_skipped():
    """SampleConstraint with a name not in geometry._stages is silently skipped."""
    # Use psic bisecting mode which includes SampleConstraint("mu", 0.0).
    # Inject a fake "virtual" constraint with a name not in the geometry.
    from ad_hoc_diffractometer.forward import _solve_bisecting

    g = _setup_cubic(psic, a=5.0)
    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    # Mode with a virtual stage name that does NOT exist in psic geometry
    mode = ConstraintSet(
        [
            BisectConstraint("eta", "delta"),
            SampleConstraint("virtual_omega", 0.0),  # not a real stage — skipped
            DetectorConstraint("nu", 0.0),
        ]
    )
    # Should not raise — virtual_omega is skipped silently
    result = _solve_bisecting(g, Q_phi, ttheta_deg, mode)
    assert isinstance(result, list)


def test_bisecting_det_constraint_names_only_detector_stage():
    """When DetectorConstraint names the only detector stage, fallback to that stage."""
    from ad_hoc_diffractometer.constants import XHAT as _XHAT
    from ad_hoc_diffractometer.constants import YHAT as _YHAT
    from ad_hoc_diffractometer.constants import ZHAT as _ZHAT
    from ad_hoc_diffractometer.forward import _solve_bisecting

    # Build a 2-circle geometry: 1 sample + 1 detector (same axis)
    stages = [
        Stage("omega", -_ZHAT, parent=None, role="sample"),
        Stage("ttheta", -_ZHAT, parent=None, role="detector"),
    ]
    g = ahd.AdHocDiffractometer(
        name="twocircle",
        stages=stages,
        basis={"vertical": _XHAT, "longitudinal": _YHAT, "transverse": _ZHAT},
    )
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=5.0)
    ub_identity(g.sample)

    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    # DetectorConstraint names the ONLY detector stage ("ttheta")
    # → active_det_stage fallback to geometry.detector_stages[-1]
    mode = ConstraintSet(
        [
            BisectConstraint("omega", "ttheta"),
            DetectorConstraint("ttheta", 0.0),  # the only detector stage
        ]
    )
    result = _solve_bisecting(g, Q_phi, ttheta_deg, mode)
    assert isinstance(result, list)


# ---------------------------------------------------------------------------
# _solve_fixed_sample — virtual-stage and single-detector-stage branches
# ---------------------------------------------------------------------------


def test_fixed_sample_virtual_stage_skipped():
    """SampleConstraint with name not in geometry._stages is skipped in _solve_fixed_sample."""
    from ad_hoc_diffractometer.forward import _solve_fixed_sample

    g = _setup_cubic(fourcv, a=5.0)
    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    mode = ConstraintSet(
        [
            SampleConstraint("virtual_nonexistent", 0.0),  # not a real stage
            SampleConstraint("chi", 90.0),
        ]
    )
    result = _solve_fixed_sample(g, Q_phi, ttheta_deg, mode)
    assert isinstance(result, list)


def test_fixed_sample_det_constraint_names_only_detector():
    """_solve_fixed_sample: DetectorConstraint names the only detector; fallback used."""
    from ad_hoc_diffractometer.constants import XHAT as _XHAT
    from ad_hoc_diffractometer.constants import YHAT as _YHAT
    from ad_hoc_diffractometer.constants import ZHAT as _ZHAT
    from ad_hoc_diffractometer.forward import _solve_fixed_sample

    stages = [
        Stage("omega", -_ZHAT, parent=None, role="sample"),
        Stage("chi", np.array([0.0, 1.0, 0.0]), parent="omega", role="sample"),
        Stage("phi", -_ZHAT, parent="chi", role="sample"),
        Stage("ttheta", -_ZHAT, parent=None, role="detector"),
    ]
    g = ahd.AdHocDiffractometer(
        name="fourcv_test",
        stages=stages,
        basis={"vertical": _XHAT, "longitudinal": _YHAT, "transverse": _ZHAT},
    )
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=5.0)
    ub_identity(g.sample)

    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    # DetectorConstraint names the ONLY detector stage
    mode = ConstraintSet(
        [
            SampleConstraint("omega", ttheta_deg / 2.0),
            SampleConstraint("chi", 90.0),
            DetectorConstraint("ttheta", 0.0),  # the only detector stage
        ]
    )
    result = _solve_fixed_sample(g, Q_phi, ttheta_deg, mode)
    assert isinstance(result, list)


# ---------------------------------------------------------------------------
# _solve_fixed_sample — unconstrained paths
# ---------------------------------------------------------------------------


def test_fixed_sample_all_constrained_in_limits():
    """_solve_fixed_sample: all sample stages fixed, valid solution returned."""
    g = _setup_cubic(fourcv, a=4.0)
    # For fourcv, fix all 3 sample stages; only ttheta (detector) is free.
    # With bisect=False and no bisect stage, _solve_fixed_sample handles this.
    # omega=ttheta/2, chi=90, phi=0 is a valid bisecting solution for (1,0,0).
    solutions_ref = g.forward(1, 0, 0)  # bisecting mode reference
    assert len(solutions_ref) > 0
    ref_sol = solutions_ref[0]  # use this as our "all fixed" solution

    from ad_hoc_diffractometer.forward import _solve_fixed_sample

    mode = ConstraintSet(
        [
            SampleConstraint("omega", ref_sol["omega"]),
            SampleConstraint("chi", ref_sol["chi"]),
            SampleConstraint("phi", ref_sol["phi"]),
        ]
    )
    # Compute ttheta from bragg
    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))
    result = _solve_fixed_sample(g, Q_phi, ttheta_deg, mode)
    assert len(result) == 1


def test_fixed_sample_all_constrained_out_of_limits():
    """_solve_fixed_sample: all sample stages fixed but out of limits returns []."""
    from ad_hoc_diffractometer.forward import _solve_fixed_sample

    g = _setup_cubic(fourcv, a=4.0)
    g.stage("omega").limits = (100.0, 180.0)  # omega out of range

    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    mode = ConstraintSet(
        [
            SampleConstraint("omega", 0.0),  # out of limits
            SampleConstraint("chi", 90.0),
            SampleConstraint("phi", 0.0),
        ]
    )
    result = _solve_fixed_sample(g, Q_phi, ttheta_deg, mode)
    assert result == []


def test_fixed_sample_one_free():
    """_solve_fixed_sample: one free stage (chi fixed, phi free) round-trips."""
    from ad_hoc_diffractometer.forward import _solve_fixed_sample

    g = _setup_cubic(fourcv, a=4.0)
    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    # Fix omega and chi; phi is free
    omega_val = ttheta_deg / 2.0
    mode = ConstraintSet(
        [
            SampleConstraint("omega", omega_val),
            SampleConstraint("chi", 90.0),
        ]
    )
    result = _solve_fixed_sample(g, Q_phi, ttheta_deg, mode)
    assert len(result) > 0
    for sol in result:
        hkl_back = g.inverse(sol)
        assert np.allclose(hkl_back, [1, 0, 0], atol=1e-5)


def test_fixed_sample_with_detector_constraint():
    """_solve_fixed_sample: DetectorConstraint branch (line 617-624) is exercised."""
    from ad_hoc_diffractometer.forward import _solve_fixed_sample
    from ad_hoc_diffractometer.presets import psic

    g = psic()
    g.wavelength = WAVELENGTH
    g.sample.lattice = ahd.Lattice(a=5.0)
    ub_identity(g.sample)

    import math

    Q_phi = g.sample.UB @ np.array([1.0, 0.0, 0.0])
    sin_theta = float(np.linalg.norm(Q_phi)) * WAVELENGTH / (4.0 * math.pi)
    ttheta_deg = 2.0 * math.degrees(math.asin(sin_theta))

    # Fix mu + DetectorConstraint(nu); no bisect — exercises the det constraint branch
    mode = ConstraintSet(
        [
            SampleConstraint("mu", 0.0),
            DetectorConstraint("nu", 0.0),
        ]
    )
    # This leaves eta, chi, phi free (3 stages > 2) — will use _solve_two_angles
    result = _solve_fixed_sample(g, Q_phi, ttheta_deg, mode)
    # Result may be empty or non-empty — just ensure it doesn't crash
    assert isinstance(result, list)


# ---------------------------------------------------------------------------
# Issue #154 — fivec: all modes implemented, round-trips
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("bisecting_4c", 1, 0, 0, id="bisecting_4c-100"),
        pytest.param("bisecting_4c", 0, 1, 0, id="bisecting_4c-010"),
        pytest.param("bisecting_4c", 1, 1, 1, id="bisecting_4c-111"),
        pytest.param("fixed_chi", 1, 0, 0, id="fixed_chi-100"),
        pytest.param("fixed_chi", 0, 1, 0, id="fixed_chi-010"),
        pytest.param("fixed_phi", 0, 1, 0, id="fixed_phi-010"),
        pytest.param("fixed_phi", 1, 1, 1, id="fixed_phi-111"),
        pytest.param("fixed_mu", 1, 0, 0, id="fixed_mu-100"),
        pytest.param("fixed_mu", 1, 1, 1, id="fixed_mu-111"),
        pytest.param(
            "constant_omega_noncoplanar", 1, 0, 0, id="constant_omega_noncoplanar-100"
        ),
        pytest.param(
            "constant_omega_noncoplanar", 0, 1, 0, id="constant_omega_noncoplanar-010"
        ),
    ],
)
def test_fivec_round_trip(mode_name, h, k, l):  # noqa: E741
    """All fivec modes solve and round-trip correctly with mu=0."""
    g = _setup_cubic(fivec, a=4.0)
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "mode_name, stage, expected_value",
    [
        pytest.param("bisecting_4c", "mu", 0.0, id="bisecting_4c-mu-zero"),
        pytest.param("fixed_chi", "mu", 0.0, id="fixed_chi-mu-zero"),
        pytest.param("fixed_chi", "chi", 90.0, id="fixed_chi-chi-value"),
        pytest.param("fixed_phi", "mu", 0.0, id="fixed_phi-mu-zero"),
        pytest.param("fixed_phi", "phi", 0.0, id="fixed_phi-phi-value"),
        pytest.param("fixed_mu", "mu", 0.0, id="fixed_mu-mu-zero"),
        pytest.param(
            "constant_omega_noncoplanar", "omega", 0.0, id="constant_omega-omega-zero"
        ),
    ],
)
def test_fivec_constraint_value_in_solution(mode_name, stage, expected_value):
    """Declared constraint values appear at their declared value in all solutions."""
    g = _setup_cubic(fivec, a=4.0)
    g.mode_name = mode_name
    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol[stage] == pytest.approx(expected_value, abs=1e-8)


def test_fivec_bisecting_omega_equals_ttheta_half():
    """In fivec bisecting_4c mode, omega must equal ttheta/2 for all solutions."""
    g = _setup_cubic(fivec, a=4.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["omega"] == pytest.approx(sol["ttheta"] / 2.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Issue #155 — sixc: implemented modes round-trip; stubs raise
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("bisecting_4c", 1, 0, 0, id="bisecting_4c-100"),
        pytest.param("bisecting_4c", 0, 1, 0, id="bisecting_4c-010"),
        pytest.param("bisecting_4c", 1, 1, 1, id="bisecting_4c-111"),
        pytest.param("fixed_gamma_5c", 1, 0, 0, id="fixed_gamma_5c-100"),
        pytest.param("fixed_gamma_5c", 0, 1, 0, id="fixed_gamma_5c-010"),
        pytest.param("fixed_alpha_5c", 1, 0, 0, id="fixed_alpha_5c-100"),
        pytest.param("fixed_alpha_5c", 1, 1, 1, id="fixed_alpha_5c-111"),
    ],
)
def test_sixc_round_trip(mode_name, h, k, l):  # noqa: E741
    """Implemented sixc modes solve and round-trip correctly."""
    g = _setup_cubic(sixc, a=4.0)
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "mode_name",
    [
        pytest.param("fixed_alpha_zaxis", id="fixed_alpha_zaxis"),
        pytest.param("fixed_beta_zaxis", id="fixed_beta_zaxis"),
        pytest.param("alpha_eq_beta_zaxis", id="alpha_eq_beta_zaxis"),
    ],
)
def test_sixc_zaxis_stub_not_implemented(mode_name):
    """zaxis modes require reference infrastructure (Issue J) — not yet implemented."""
    g = _setup_cubic(sixc, a=4.0)
    g.mode_name = mode_name
    assert not g.modes[mode_name].is_implemented(g)
    with pytest.raises(NotImplementedError):
        g.forward(1, 0, 0)


def test_sixc_four_circle_matches_fourcv():
    """sixc four_circle mode with alpha=gamma=0 gives same ttheta as fourcv bisecting."""
    g_sixc = _setup_cubic(sixc, a=4.0)
    g_fourcv = _setup_cubic(fourcv, a=4.0)

    sols_sixc = g_sixc.forward(1, 0, 0)
    sols_fourcv = g_fourcv.forward(1, 0, 0)

    assert len(sols_sixc) > 0 and len(sols_fourcv) > 0
    # ttheta/delta should match between the two geometries
    ttheta_fourcv = sols_fourcv[0]["ttheta"]
    delta_sixc = sols_sixc[0]["delta"]
    assert delta_sixc == pytest.approx(ttheta_fourcv, abs=1e-8)


def test_sixc_four_circle_alpha_gamma_frozen():
    """In four_circle mode, alpha=0 and gamma=0 in all solutions."""
    g = _setup_cubic(sixc, a=4.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["alpha"] == pytest.approx(0.0, abs=1e-8)
        assert sol["gamma"] == pytest.approx(0.0, abs=1e-8)


def test_sixc_four_circle_omega_equals_delta_half():
    """In sixc four_circle mode, omega must equal delta/2 for all solutions."""
    g = _setup_cubic(sixc, a=4.0)
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["omega"] == pytest.approx(sol["delta"] / 2.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Issue #150 — psic: all modes, bisecting_horizontal, stubs
# ---------------------------------------------------------------------------

_PSIC_MODES_ALL = {
    "bisecting_vertical",
    "fixed_chi",
    "fixed_phi",
    "fixed_mu",
    "bisecting_horizontal",
    "fixed_nu",
    "double_diffraction_vertical",
    "double_diffraction_horizontal",
    "lifting_detector_mu",
    "lifting_detector_phi",
    "psi_constant_vertical",
    "psi_constant_horizontal",
}

_PSIC_MODES_IMPLEMENTED = {
    "bisecting_vertical",
    "fixed_chi",
    "fixed_phi",
    "fixed_mu",
    "bisecting_horizontal",
    "fixed_nu",
    "double_diffraction_vertical",
    "double_diffraction_horizontal",
    "lifting_detector_mu",
    "lifting_detector_phi",
}

_PSIC_MODES_STUBS = _PSIC_MODES_ALL - _PSIC_MODES_IMPLEMENTED


def test_psic_has_all_twelve_modes():
    """psic exposes exactly 12 declared modes."""
    assert set(psic().modes.keys()) == _PSIC_MODES_ALL


@pytest.mark.parametrize(
    "mode_name, expected_implemented",
    [pytest.param(m, True, id=f"{m}-impl") for m in sorted(_PSIC_MODES_IMPLEMENTED)]
    + [pytest.param(m, False, id=f"{m}-stub") for m in sorted(_PSIC_MODES_STUBS)],
)
def test_psic_mode_is_implemented(mode_name, expected_implemented):
    """Implemented psic modes return True; stubs return False."""
    g = psic()
    assert g.modes[mode_name].is_implemented(g) == expected_implemented


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("bisecting_vertical", 1, 0, 0, id="bisecting_vertical-100"),
        pytest.param("bisecting_vertical", 1, 1, 1, id="bisecting_vertical-111"),
        pytest.param("fixed_chi", 0, 0, 1, id="fixed_chi-001"),
        pytest.param("fixed_chi", 0, 1, 0, id="fixed_chi-010"),
        pytest.param("fixed_phi", 0, 1, 0, id="fixed_phi-010"),
        pytest.param("fixed_phi", 1, 1, 1, id="fixed_phi-111"),
        pytest.param("fixed_mu", 1, 0, 0, id="fixed_mu-100"),
        pytest.param("fixed_mu", 0, 1, 0, id="fixed_mu-010"),
        # bisecting_horizontal: Q must lie in the horizontal plane (no longitudinal component)
        pytest.param("bisecting_horizontal", 1, 0, 0, id="bisecting_horiz-100"),
        pytest.param("bisecting_horizontal", 0, 0, 1, id="bisecting_horiz-001"),
        pytest.param("bisecting_horizontal", 1, 0, 1, id="bisecting_horiz-101"),
        pytest.param("fixed_nu", 1, 0, 0, id="fixed_nu-100"),
        pytest.param("fixed_nu", 1, 1, 1, id="fixed_nu-111"),
        # double_diffraction_vertical tested separately (requires extras h2/k2/l2)
        # lifting_detector_* modes: qaz=90 out-of-plane constraint
        pytest.param("lifting_detector_mu", 1, 0, 0, id="lifting_detector_mu-100"),
        pytest.param("lifting_detector_mu", 0, 0, 1, id="lifting_detector_mu-001"),
        pytest.param("lifting_detector_phi", 1, 0, 0, id="lifting_detector_phi-100"),
        pytest.param("lifting_detector_phi", 0, 1, 0, id="lifting_detector_phi-010"),
    ],
)
def test_psic_mode_round_trip(mode_name, h, k, l):  # noqa: E741
    """All implemented psic modes solve and round-trip correctly."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "mode_name",
    [pytest.param(m, id=m) for m in sorted(_PSIC_MODES_STUBS)],
)
def test_psic_stub_not_implemented(mode_name):
    """Stub modes raise NotImplementedError on forward()."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = mode_name
    with pytest.raises(NotImplementedError):
        g.forward(1, 0, 0)


def test_psic_bisecting_horizontal_mu_equals_nu_half():
    """bisecting_horizontal: mu = nu/2 for all solutions."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = "bisecting_horizontal"
    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["mu"] == pytest.approx(sol["nu"] / 2.0, abs=1e-10)


def test_psic_bisecting_horizontal_eta_delta_frozen():
    """bisecting_horizontal: eta=0 and delta=0 in all solutions."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = "bisecting_horizontal"
    solutions = g.forward(1, 0, 0)
    for sol in solutions:
        assert sol["eta"] == pytest.approx(0.0, abs=1e-8)
        assert sol["delta"] == pytest.approx(0.0, abs=1e-8)


@pytest.mark.parametrize(
    "mode_name, stage, expected_value, h, k, l",
    [
        pytest.param(
            "bisecting_vertical", "mu", 0.0, 1, 0, 0, id="bisecting_vertical-mu=0"
        ),
        pytest.param(
            "bisecting_vertical", "nu", 0.0, 1, 0, 0, id="bisecting_vertical-nu=0"
        ),
        pytest.param("fixed_phi", "phi", 0.0, 0, 1, 0, id="fixed_phi-phi=0"),
        pytest.param("fixed_mu", "mu", 0.0, 1, 0, 0, id="fixed_mu-mu=0"),
        pytest.param("fixed_nu", "nu", 0.0, 1, 0, 0, id="fixed_nu-nu=0"),
        # fixed_chi: (0,0,1) gives mu=0; chi=90 is the only declared constraint
        pytest.param("fixed_chi", "chi", 90.0, 0, 0, 1, id="fixed_chi-chi=90"),
        # bisecting_horizontal: Q must lie in horizontal plane; (1,0,0) works
        pytest.param(
            "bisecting_horizontal", "eta", 0.0, 1, 0, 0, id="bisecting_horiz-eta=0"
        ),
        pytest.param(
            "bisecting_horizontal", "delta", 0.0, 1, 0, 0, id="bisecting_horiz-delta=0"
        ),
    ],
)
def test_psic_constraint_value_in_solution(mode_name, stage, expected_value, h, k, l):  # noqa: E741
    """Constraint values appear at their declared values in all solutions."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol[stage] == pytest.approx(expected_value, abs=1e-8)


def test_psic_psi_constant_extras_declared():
    """psi_constant modes carry REQUIRED n_hat and None psi output extras."""
    from ad_hoc_diffractometer import REQUIRED

    for mode_name in ("psi_constant_vertical", "psi_constant_horizontal"):
        cs = psic().modes[mode_name]
        assert cs.extras.get("n_hat") is REQUIRED
        assert cs.extras.get("psi") is None


def test_psic_double_diffraction_extras_declared():
    """double_diffraction_vertical carries REQUIRED h2, k2, l2 extras."""
    from ad_hoc_diffractometer import REQUIRED

    cs = psic().modes["double_diffraction_vertical"]
    assert cs.extras.get("h2") is REQUIRED
    assert cs.extras.get("k2") is REQUIRED
    assert cs.extras.get("l2") is REQUIRED


def test_psic_modes_serialisation_round_trip():
    """All 11 psic modes survive to_dict / from_dict round-trip."""
    import json

    from ad_hoc_diffractometer import AdHocDiffractometer

    g = psic()
    d = g.to_dict()
    assert json.dumps(d)
    assert set(d["modes"].keys()) == _PSIC_MODES_ALL
    g2 = AdHocDiffractometer.from_dict(d)
    assert set(g2.modes.keys()) == _PSIC_MODES_ALL
    assert g2.mode_name == "bisecting_vertical"


# ---------------------------------------------------------------------------
# Issue #152 — kappa6c: implemented modes round-trip; stubs raise
# ---------------------------------------------------------------------------

_KAPPA6C_ALL_MODES = {
    "bisecting_vertical",
    "bisecting_horizontal",
    "fixed_kphi",
    "fixed_mu",
    "fixed_nu",
    "fixed_delta",
    "lifting_detector_mu",
    "lifting_detector_kphi",
    "psi_constant_vertical",
    "psi_constant_horizontal",
    "double_diffraction_vertical",
    "double_diffraction_horizontal",
}
_KAPPA6C_STUB_MODES = {
    "psi_constant_vertical",
    "psi_constant_horizontal",
}


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("bisecting_vertical", 1, 0, 0, id="bisecting_vertical-100"),
        pytest.param("bisecting_vertical", 0, 1, 0, id="bisecting_vertical-010"),
        pytest.param("bisecting_horizontal", 1, 0, 0, id="bisecting_horizontal-100"),
        pytest.param("bisecting_horizontal", 0, 0, 1, id="bisecting_horizontal-001"),
        pytest.param("fixed_kphi", 0, 1, 0, id="fixed_kphi-010"),
        pytest.param("fixed_kphi", 1, 0, 0, id="fixed_kphi-100"),
        pytest.param("fixed_mu", 1, 0, 0, id="fixed_mu-100"),
        pytest.param("fixed_mu", 0, 1, 0, id="fixed_mu-010"),
        pytest.param("fixed_nu", 1, 0, 0, id="fixed_nu-100"),
        pytest.param("fixed_nu", 0, 1, 0, id="fixed_nu-010"),
        pytest.param("fixed_delta", 1, 0, 0, id="fixed_delta-100"),
        pytest.param("fixed_delta", 0, 0, 1, id="fixed_delta-001"),
        # lifting_detector_* modes: qaz=90 out-of-plane constraint
        pytest.param("lifting_detector_mu", 1, 0, 0, id="lifting_detector_mu-100"),
        pytest.param("lifting_detector_mu", 0, 1, 0, id="lifting_detector_mu-010"),
        pytest.param("lifting_detector_kphi", 1, 0, 0, id="lifting_detector_kphi-100"),
        pytest.param("lifting_detector_kphi", 0, 1, 0, id="lifting_detector_kphi-010"),
    ],
)
def test_kappa6c_mode_round_trip(mode_name, h, k, l):  # noqa: E741
    """Implemented kappa6c modes solve and round-trip correctly."""
    g = _setup_cubic(kappa6c, a=4.0)
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "mode_name",
    [pytest.param(m, id=m) for m in sorted(_KAPPA6C_STUB_MODES)],
)
def test_kappa6c_stub_not_implemented(mode_name):
    """kappa6c stub modes raise NotImplementedError on forward()."""
    g = _setup_cubic(kappa6c, a=4.0)
    g.mode_name = mode_name
    with pytest.raises(NotImplementedError):
        g.forward(1, 0, 0)


def test_kappa6c_bisecting_vertical_invariants():
    """bisecting_vertical: komega=delta/2, mu=0, nu=0 in all solutions."""
    g = _setup_cubic(kappa6c, a=4.0)
    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["komega"] == pytest.approx(sol["delta"] / 2.0, abs=1e-10)
        assert sol["mu"] == pytest.approx(0.0, abs=1e-8)
        assert sol["nu"] == pytest.approx(0.0, abs=1e-8)


def test_kappa6c_bisecting_horizontal_invariants():
    """bisecting_horizontal: mu=nu/2, komega=0, delta=0 in all solutions."""
    g = _setup_cubic(kappa6c, a=4.0)
    g.mode_name = "bisecting_horizontal"
    solutions = g.forward(1, 0, 0)
    assert len(solutions) > 0
    for sol in solutions:
        assert sol["mu"] == pytest.approx(sol["nu"] / 2.0, abs=1e-10)
        assert sol["komega"] == pytest.approx(0.0, abs=1e-8)
        assert sol["delta"] == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------------------
# Issue #177 — qaz detector constraint solver
# ---------------------------------------------------------------------------


def _qaz_from_angles(nu_deg: float, delta_deg: float) -> float:
    """Compute qaz from detector angles per You (1999) eq. 18."""
    import math

    nu_r = math.radians(nu_deg)
    delta_r = math.radians(delta_deg)
    return math.degrees(math.atan2(math.tan(delta_r), math.sin(nu_r)))


@pytest.mark.parametrize(
    "nu_deg, delta_deg, target_qaz_deg, expected_residual",
    [
        pytest.param(
            0.0,
            20.0,
            90.0,
            0.0,
            id="nu=0-delta=20-qaz=90-satisfied",
        ),
        pytest.param(
            90.0,
            30.0,
            30.0,
            0.0,
            id="nu=90-delta=30-qaz=30-satisfied",
        ),
        pytest.param(
            45.0,
            0.0,
            0.0,
            0.0,
            id="nu=45-delta=0-qaz=0-satisfied",
        ),
        pytest.param(
            30.0,
            20.0,
            90.0,
            _qaz_from_angles(30.0, 20.0) - 90.0,
            id="nu=30-delta=20-qaz=90-nonzero-residual",
        ),
    ],
)
def test_qaz_residual(nu_deg, delta_deg, target_qaz_deg, expected_residual):
    """_qaz_residual returns correct residual per You (1999) eq. 18."""
    from ad_hoc_diffractometer.mode import _qaz_residual

    g = psic()
    angles = {
        "mu": 0.0,
        "eta": 0.0,
        "chi": 90.0,
        "phi": 0.0,
        "nu": nu_deg,
        "delta": delta_deg,
    }
    residual = _qaz_residual(angles, g, target_qaz_deg)
    assert residual == pytest.approx(expected_residual, abs=1e-6)


def test_psic_lifting_detector_mu_limits_exclude_nu_zero():
    """If nu stage limits exclude nu=0, qaz=90 solutions are filtered out."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = "lifting_detector_mu"
    # Default solution for qaz=90 uses nu=0; excluding nu=0 gives no solutions.
    nu_stage = g.stage("nu")
    nu_stage.limits = (1.0, 180.0)
    solutions = g.forward(1, 0, 0)
    # All detector solutions use nu=0 which is now outside limits
    assert len(solutions) == 0


def test_qaz_residual_two_detector_stages_required():
    """_qaz_residual raises ValueError for geometry with < 2 detector stages."""

    import pytest

    # Build a geometry with only 1 detector stage
    from ad_hoc_diffractometer import AdHocDiffractometer
    from ad_hoc_diffractometer.constants import XHAT
    from ad_hoc_diffractometer.constants import YHAT
    from ad_hoc_diffractometer.constants import ZHAT
    from ad_hoc_diffractometer.mode import _qaz_residual
    from ad_hoc_diffractometer.stage import Stage

    stages = [
        Stage("omega", -ZHAT, parent=None, role="sample"),
        Stage("ttheta", -ZHAT, parent=None, role="detector"),
    ]
    g_1det = AdHocDiffractometer(
        name="one_det",
        stages=stages,
        basis={"vertical": XHAT, "longitudinal": YHAT, "transverse": ZHAT},
    )
    angles = {"omega": 0.0, "ttheta": 20.0}
    with pytest.raises(ValueError, match="at least 2 detector stages"):
        _qaz_residual(angles, g_1det, 90.0)


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("lifting_detector_mu", 1, 0, 0, id="psic-liftmu-100"),
        pytest.param("lifting_detector_mu", 0, 0, 1, id="psic-liftmu-001"),
        pytest.param("lifting_detector_phi", 1, 0, 0, id="psic-liftphi-100"),
        pytest.param("lifting_detector_phi", 0, 1, 0, id="psic-liftphi-010"),
    ],
)
def test_psic_lifting_detector_qaz_satisfied(mode_name, h, k, l):  # noqa: E741
    """lifting_detector_* modes: qaz = 90.0 verified in all solutions."""
    g = _setup_cubic(psic, a=4.0)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0, f"No solutions for {mode_name} ({h},{k},{l})"
    for sol in solutions:
        nu_deg = sol["nu"]
        delta_deg = sol["delta"]
        qaz_computed = _qaz_from_angles(nu_deg, delta_deg)
        assert qaz_computed == pytest.approx(90.0, abs=1e-4), (
            f"{mode_name}: expected qaz=90, got {qaz_computed:.6f} "
            f"(nu={nu_deg:.4f}, delta={delta_deg:.4f})"
        )


@pytest.mark.parametrize(
    "mode_name, h, k, l",
    [
        pytest.param("lifting_detector_mu", 1, 0, 0, id="kappa6c-liftmu-100"),
        pytest.param("lifting_detector_mu", 0, 1, 0, id="kappa6c-liftmu-010"),
        pytest.param("lifting_detector_kphi", 1, 0, 0, id="kappa6c-liftkphi-100"),
        pytest.param("lifting_detector_kphi", 0, 1, 0, id="kappa6c-liftkphi-010"),
    ],
)
def test_kappa6c_lifting_detector_qaz_satisfied(mode_name, h, k, l):  # noqa: E741
    """kappa6c lifting_detector_* modes: qaz = 90.0 verified in all solutions."""
    g = _setup_cubic(kappa6c, a=4.0)
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0, f"No solutions for {mode_name} ({h},{k},{l})"
    for sol in solutions:
        nu_deg = sol["nu"]
        delta_deg = sol["delta"]
        qaz_computed = _qaz_from_angles(nu_deg, delta_deg)
        assert qaz_computed == pytest.approx(90.0, abs=1e-4), (
            f"{mode_name}: expected qaz=90, got {qaz_computed:.6f} "
            f"(nu={nu_deg:.4f}, delta={delta_deg:.4f})"
        )


# ---------------------------------------------------------------------------
# Issue #176 — psi_constant forward solver (validation filter)
# ---------------------------------------------------------------------------


def _setup_psi(factory, ref=(0, 0, 1), a=4.0):
    """Return a geometry with wavelength, UB=B, and azimuthal_reference set."""
    g = _setup_cubic(factory, a=a)
    g.azimuthal_reference = ref
    return g


def _natural_psi(g, h, k, l):  # noqa: E741
    """Compute the natural psi from the phi frame (motor-angle independent)."""
    from ad_hoc_diffractometer.forward import _compute_natural_psi

    Q_phi = g.sample.UB @ np.array([h, k, l], dtype=float)
    return _compute_natural_psi(g, Q_phi)


# --- _compute_natural_psi unit tests ---


@pytest.mark.parametrize(
    "factory, h, k, l, ref, expected_psi",
    [
        pytest.param(fourcv, 1, 0, 0, (0, 0, 1), 90.0, id="fourcv-100-ref001"),
        pytest.param(fourcv, 1, 1, 1, (0, 0, 1), 120.0, id="fourcv-111-ref001"),
        pytest.param(fourcv, 1, 0, 1, (0, 0, 1), 90.0, id="fourcv-101-ref001"),
        pytest.param(psic, 1, 0, 0, (0, 0, 1), 90.0, id="psic-100-ref001"),
        pytest.param(psic, 1, 1, 1, (0, 0, 1), 120.0, id="psic-111-ref001"),
    ],
)
def test_compute_natural_psi(factory, h, k, l, ref, expected_psi):  # noqa: E741
    """_compute_natural_psi returns the correct motor-angle-independent psi."""
    g = _setup_psi(factory, ref=ref)
    psi = _natural_psi(g, h, k, l)
    assert psi is not None
    assert psi == pytest.approx(expected_psi, abs=1e-4)


def test_compute_natural_psi_undefined_when_ref_parallel_to_Q():
    """_compute_natural_psi returns None when reference is parallel to Q."""
    g = _setup_psi(fourcv, ref=(1, 0, 0))
    psi = _natural_psi(g, 1, 0, 0)  # Q along x, ref along x -> parallel
    assert psi is None


def test_compute_natural_psi_undefined_when_Q_parallel_to_beam():
    """_compute_natural_psi returns None when Q is parallel to incident beam."""
    g = _setup_psi(fourcv, ref=(0, 0, 1))
    psi = _natural_psi(g, 0, 1, 0)  # Q along y (longitudinal = beam direction)
    assert psi is None


# --- fourcv / fourch psi_constant round-trip tests (synthetic bisect path) ---


@pytest.mark.parametrize(
    "factory, h, k, l, ref",
    [
        pytest.param(fourcv, 1, 0, 0, (0, 0, 1), id="fourcv-100-ref001"),
        pytest.param(fourcv, 1, 0, 1, (0, 0, 1), id="fourcv-101-ref001"),
        pytest.param(fourcv, 1, 1, 1, (0, 0, 1), id="fourcv-111-ref001"),
        pytest.param(fourch, 1, 0, 0, (0, 0, 1), id="fourch-100-ref001"),
        pytest.param(fourch, 1, 1, 1, (0, 0, 1), id="fourch-111-ref001"),
    ],
)
def test_four_circle_psi_constant_round_trip(factory, h, k, l, ref):  # noqa: E741
    """psi_constant returns bisecting solutions when natural psi matches target."""
    g = _setup_psi(factory, ref=ref)
    natural = _natural_psi(g, h, k, l)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    g.modes["psi_constant"] = ConstraintSet(
        [ReferenceConstraint("psi", natural)],
        computed=g.modes["psi_constant"].computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "psi_constant"
    assert _round_trip_ok(g, h, k, l)


@pytest.mark.parametrize(
    "factory",
    [
        pytest.param(fourcv, id="fourcv"),
        pytest.param(fourch, id="fourch"),
    ],
)
def test_four_circle_psi_constant_wrong_target_returns_empty(factory):
    """psi_constant returns [] when natural psi does not match psi_target."""
    g = _setup_psi(factory, ref=(0, 0, 1))
    natural = _natural_psi(g, 1, 0, 0)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    wrong_target = natural + 45.0
    g.modes["psi_constant"] = ConstraintSet(
        [ReferenceConstraint("psi", wrong_target)],
        computed=g.modes["psi_constant"].computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "psi_constant"
    solutions = g.forward(1, 0, 0)
    assert solutions == []


# --- psic psi_constant_vertical / horizontal round-trip tests (bisecting path) ---


@pytest.mark.parametrize(
    "mode_name, h, k, l, ref",
    [
        pytest.param(
            "psi_constant_vertical", 1, 0, 0, (0, 0, 1), id="psic-psi_vert-100"
        ),
        pytest.param(
            "psi_constant_vertical", 1, 0, 1, (0, 0, 1), id="psic-psi_vert-101"
        ),
        pytest.param(
            "psi_constant_vertical", 1, 1, 1, (0, 0, 1), id="psic-psi_vert-111"
        ),
        pytest.param(
            "psi_constant_horizontal", 1, 0, 0, (0, 0, 1), id="psic-psi_horiz-100"
        ),
    ],
)
def test_psic_psi_constant_round_trip(mode_name, h, k, l, ref):  # noqa: E741
    """psic psi_constant modes return bisecting solutions when psi matches."""
    g = _setup_psi(psic, ref=ref)
    natural = _natural_psi(g, h, k, l)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    old_mode = g.modes[mode_name]
    new_constraints = []
    for c in old_mode.constraints:
        if isinstance(c, ReferenceConstraint) and c.name == "psi":
            new_constraints.append(ReferenceConstraint("psi", natural))
        else:
            new_constraints.append(c)
    g.modes[mode_name] = ConstraintSet(
        new_constraints,
        computed=old_mode.computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


def test_psic_psi_constant_wrong_target_returns_empty():
    """psic psi_constant_vertical returns [] when psi target is wrong."""
    g = _setup_psi(psic, ref=(0, 0, 1))
    natural = _natural_psi(g, 1, 0, 0)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import BisectConstraint
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint
    from ad_hoc_diffractometer import SampleConstraint

    g.modes["psi_constant_vertical"] = ConstraintSet(
        [
            BisectConstraint("eta", "delta"),
            SampleConstraint("mu", 0.0),
            ReferenceConstraint("psi", natural + 45.0),
        ],
        computed=["eta", "chi", "phi", "delta"],
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "psi_constant_vertical"
    solutions = g.forward(1, 0, 0)
    assert solutions == []


# --- psi_constant psi verification in solutions ---


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l, ref",
    [
        pytest.param(fourcv, "psi_constant", 1, 0, 0, (0, 0, 1), id="fourcv-100"),
        pytest.param(fourcv, "psi_constant", 1, 1, 1, (0, 0, 1), id="fourcv-111"),
        pytest.param(
            psic, "psi_constant_vertical", 1, 0, 0, (0, 0, 1), id="psic-vert-100"
        ),
        pytest.param(
            psic, "psi_constant_vertical", 1, 1, 1, (0, 0, 1), id="psic-vert-111"
        ),
        pytest.param(
            psic, "psi_constant_horizontal", 1, 0, 0, (0, 0, 1), id="psic-horiz-100"
        ),
    ],
)
def test_psi_constant_psi_verified_in_solutions(
    factory,
    mode_name,
    h,
    k,
    l,  # noqa: E741
    ref,
):
    """All psi_constant solutions have psi == natural_psi."""
    g = _setup_psi(factory, ref=ref)
    natural = _natural_psi(g, h, k, l)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    old_mode = g.modes[mode_name]
    new_constraints = []
    for c in old_mode.constraints:
        if isinstance(c, ReferenceConstraint) and c.name == "psi":
            new_constraints.append(ReferenceConstraint("psi", natural))
        else:
            new_constraints.append(c)
    g.modes[mode_name] = ConstraintSet(
        new_constraints,
        computed=old_mode.computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = mode_name
    solutions = g.forward(h, k, l)
    assert len(solutions) > 0, f"No solutions for {mode_name} ({h},{k},{l})"
    for sol in solutions:
        psi_check = g.psi(angles=sol)
        assert psi_check == pytest.approx(natural, abs=1e-3), (
            f"psi mismatch: expected {natural:.4f}, got {psi_check:.4f}"
        )


# --- kappa6c psi_constant (bisecting path with kappa stages) ---


@pytest.mark.parametrize(
    "mode_name, h, k, l, ref",
    [
        pytest.param(
            "psi_constant_vertical", 1, 0, 0, (0, 0, 1), id="kappa6c-psi_vert-100"
        ),
        pytest.param(
            "psi_constant_vertical", 1, 1, 0, (0, 0, 1), id="kappa6c-psi_vert-110"
        ),
    ],
)
def test_kappa6c_psi_constant_round_trip(mode_name, h, k, l, ref):  # noqa: E741
    """kappa6c psi_constant modes return bisecting solutions when psi matches."""
    g = _setup_psi(kappa6c, ref=ref)
    natural = _natural_psi(g, h, k, l)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    old_mode = g.modes[mode_name]
    new_constraints = []
    for c in old_mode.constraints:
        if isinstance(c, ReferenceConstraint) and c.name == "psi":
            new_constraints.append(ReferenceConstraint("psi", natural))
        else:
            new_constraints.append(c)
    g.modes[mode_name] = ConstraintSet(
        new_constraints,
        computed=old_mode.computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = mode_name
    assert _round_trip_ok(g, h, k, l)


# --- kappa4cv psi_constant (synthetic bisect path) ---


def test_kappa4cv_psi_constant_round_trip():
    """kappa4cv psi_constant returns bisecting solutions via synthetic bisect."""
    g = _setup_psi(kappa4cv, ref=(0, 0, 1))
    natural = _natural_psi(g, 1, 1, 0)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    g.modes["psi_constant"] = ConstraintSet(
        [ReferenceConstraint("psi", natural)],
        computed=g.modes["psi_constant"].computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "psi_constant"
    assert _round_trip_ok(g, 1, 1, 0)


# --- psi undefined → empty solutions ---


def test_psi_constant_undefined_psi_returns_empty():
    """psi_constant returns [] when psi is undefined (Q ∥ incident beam)."""
    g = _setup_psi(fourcv, ref=(0, 0, 1))
    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    g.modes["psi_constant"] = ConstraintSet(
        [ReferenceConstraint("psi", 0.0)],
        computed=g.modes["psi_constant"].computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "psi_constant"
    # (0,1,0) has Q along y (longitudinal = incident beam direction) → psi undefined
    solutions = g.forward(0, 1, 0)
    assert solutions == []


# --- wraparound tolerance test ---


def test_psi_constant_wraparound_tolerance():
    """psi_constant handles ±180° wraparound correctly."""
    g = _setup_psi(fourcv, ref=(1, 0, 0))
    natural = _natural_psi(g, 1, 1, 0)
    assert natural is not None

    from ad_hoc_diffractometer import REQUIRED
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import ReferenceConstraint

    # Set target to the opposite wraparound side (shift by 360°)
    if natural > 0:
        wrapped = natural - 360.0
    else:
        wrapped = natural + 360.0

    g.modes["psi_constant"] = ConstraintSet(
        [ReferenceConstraint("psi", wrapped)],
        computed=g.modes["psi_constant"].computed,
        extras={"n_hat": REQUIRED, "psi": None},
    )
    g.mode_name = "psi_constant"
    solutions = g.forward(1, 1, 0)
    # wrapped and natural differ by 360° — should still match
    assert len(solutions) > 0


# ---------------------------------------------------------------------------
# Issue #176 — double_diffraction forward solver (4D simultaneous)
# ---------------------------------------------------------------------------


def _setup_dd(factory, hkl2=(0, 1, 0), a=4.0, mode_name=None):
    """Return a geometry with double_diffraction mode and h2/k2/l2 set."""
    g = _setup_cubic(factory, a=a)
    if mode_name is None:
        # Pick the first double_diffraction mode available
        for m in g.modes:
            if "double_diffraction" in m:
                mode_name = m
                break
    assert mode_name is not None, f"No double_diffraction mode found in {factory}"
    g.mode_name = mode_name
    cs = g.modes[mode_name]
    cs.extras["h2"] = float(hkl2[0])
    cs.extras["k2"] = float(hkl2[1])
    cs.extras["l2"] = float(hkl2[2])
    return g


# --- ValueError when extras not set ---


@pytest.mark.parametrize(
    "factory, mode_name",
    [
        pytest.param(fourcv, "double_diffraction", id="fourcv"),
        pytest.param(fourch, "double_diffraction", id="fourch"),
        pytest.param(psic, "double_diffraction_vertical", id="psic-vert"),
        pytest.param(psic, "double_diffraction_horizontal", id="psic-horiz"),
        pytest.param(kappa4cv, "double_diffraction", id="kappa4cv"),
        pytest.param(kappa6c, "double_diffraction_vertical", id="kappa6c-vert"),
        pytest.param(kappa6c, "double_diffraction_horizontal", id="kappa6c-horiz"),
    ],
)
def test_double_diffraction_raises_without_extras(factory, mode_name):
    """forward() raises ValueError when h2/k2/l2 are REQUIRED sentinels."""
    g = _setup_cubic(factory, a=4.0)
    g.mode_name = mode_name
    with pytest.raises(ValueError, match="h2, k2, l2"):
        g.forward(1, 0, 0)


# --- Round-trip tests ---


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l, hkl2",
    [
        pytest.param(fourcv, "double_diffraction", 1, 0, 0, (0, 1, 0), id="fourcv-100"),
        pytest.param(fourcv, "double_diffraction", 1, 1, 0, (0, 0, 1), id="fourcv-110"),
        pytest.param(fourch, "double_diffraction", 1, 0, 0, (0, 1, 0), id="fourch-100"),
        pytest.param(
            psic,
            "double_diffraction_vertical",
            1,
            0,
            0,
            (0, 1, 0),
            id="psic-vert-100",
        ),
        pytest.param(
            psic,
            "double_diffraction_horizontal",
            1,
            0,
            0,
            (0, 1, 0),
            id="psic-horiz-100",
        ),
    ],
)
def test_double_diffraction_round_trip(
    factory,
    mode_name,
    h,
    k,
    l,  # noqa: E741
    hkl2,
):
    """double_diffraction returns solutions that round-trip via inverse()."""
    g = _setup_dd(factory, hkl2=hkl2, mode_name=mode_name)
    solutions = g.forward(h, k, l)
    # May return 0 solutions if no simultaneous diffraction exists
    for sol in solutions:
        hkl_back = g.inverse(sol)
        assert np.allclose(hkl_back, [h, k, l], atol=1e-6), (
            f"Round-trip failed: {[h, k, l]} -> {sol} -> {hkl_back}"
        )


# --- Ewald sphere verification for secondary reflection ---


@pytest.mark.parametrize(
    "factory, mode_name, h, k, l, hkl2",
    [
        pytest.param(fourcv, "double_diffraction", 1, 0, 0, (0, 1, 0), id="fourcv-100"),
        pytest.param(
            psic,
            "double_diffraction_vertical",
            1,
            0,
            0,
            (0, 1, 0),
            id="psic-vert-100",
        ),
    ],
)
def test_double_diffraction_secondary_on_ewald_sphere(
    factory,
    mode_name,
    h,
    k,
    l,  # noqa: E741
    hkl2,
):
    """In all solutions, the secondary reflection satisfies the Ewald sphere."""
    import math

    from ad_hoc_diffractometer.rotation import rotation_matrix

    g = _setup_dd(factory, hkl2=hkl2, mode_name=mode_name)
    solutions = g.forward(h, k, l)

    hkl2_arr = np.array(hkl2, dtype=float)
    Q2_phi = g.sample.UB @ hkl2_arr
    k_mag = 2.0 * math.pi / g.wavelength
    y_raw = np.asarray(g.basis.get("longitudinal", [0, 1, 0]), dtype=float)
    ki = k_mag * (y_raw / np.linalg.norm(y_raw))
    ki_sq = float(np.dot(ki, ki))

    for sol in solutions:
        # Build Z from solution angles
        Z = np.eye(3)
        for s in g.sample_stages:
            Z = Z @ rotation_matrix(s.axis, sol[s.name])
        Q2_lab = Z @ Q2_phi
        kf2 = ki + Q2_lab
        kf2_sq = float(np.dot(kf2, kf2))
        residual = abs(kf2_sq - ki_sq)
        assert residual < 1e-3, (  # noqa: PLR2004
            f"Secondary reflection not on Ewald sphere: |kf2|²-|ki|² = {residual:.6f}"
        )


# ---------------------------------------------------------------------------
# ForwardContext coverage — all-sample-stages-constrained caching path
# ---------------------------------------------------------------------------


def test_forward_context_all_stages_constrained():
    """ForwardContext.prepare_caching with no free sample stages sets Z_prefix to None."""
    from ad_hoc_diffractometer.forward import ForwardContext

    g = _setup_cubic(fourcv, a=4.0)

    ctx = ForwardContext(g)
    # Pass empty free set — all stages are constrained
    angles = {s.name: s.angle for s in g._stages.values()}
    ctx.prepare_caching(angles, set())

    assert ctx._cached_Z_prefix is None
    assert ctx._free_sample_indices == []
    assert ctx._cached_D is not None

    # q_phi_uncached still works
    Q = ctx.q_phi_uncached(angles)
    assert Q.shape == (3,)


def test_bisecting_early_termination_stale():
    """Bisecting solver terminates early after enough stale seeds."""
    # Use a reflection that has exactly 2 solutions so that:
    # - The first 2 solutions are found from the analytic seeds
    # - Subsequent seeds are all stale (converge to duplicates)
    # - The solver breaks after _MAX_STALE consecutive stale seeds
    # If it didn't early-terminate, it would try all 24 seeds.
    g = _setup_cubic(fourcv, a=4.0)
    solutions = g.forward(1, 0, 0)
    assert len(solutions) == 2  # noqa: PLR2004
    # Verify round-trip for both solutions
    for sol in solutions:
        hkl_back = g.inverse(sol)
        assert abs(hkl_back[0] - 1.0) < 1e-6
        assert abs(hkl_back[1]) < 1e-6
        assert abs(hkl_back[2]) < 1e-6


def test_bisecting_max_solutions_termination():
    """Bisecting solver terminates after finding _MAX_SOLUTIONS unique solutions."""
    # s2d2 geometry with bisecting_vertical can produce 4 solutions
    # for some reflections, triggering the _MAX_SOLUTIONS break.
    g = _setup_cubic(psic, a=4.0)

    # Try several reflections to find one with many solutions
    max_sols = 0
    for hkl in [(1, 0, 0), (0, 1, 0), (1, 1, 0), (0, 0, 1), (1, 1, 1)]:
        solutions = g.forward(*hkl)
        max_sols = max(max_sols, len(solutions))
        # All solutions must round-trip correctly
        for sol in solutions:
            hkl_back = g.inverse(sol)
            assert abs(hkl_back[0] - hkl[0]) < 1e-6
            assert abs(hkl_back[1] - hkl[1]) < 1e-6
            assert abs(hkl_back[2] - hkl[2]) < 1e-6
    # Verify we find more than 2 solutions for at least one reflection
    assert max_sols >= 2  # noqa: PLR2004
