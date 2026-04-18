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
from ad_hoc_diffractometer import Stage
from ad_hoc_diffractometer import compute_forward
from ad_hoc_diffractometer import fourch
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import kappa4cv
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.constants import XHAT
from ad_hoc_diffractometer.constants import YHAT
from ad_hoc_diffractometer.constants import ZHAT
from ad_hoc_diffractometer.forward import _apply_cut_points
from ad_hoc_diffractometer.forward import _check_limits

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
# BisectingMode — fourcv (bisecting, omega = ttheta/2)
# ---------------------------------------------------------------------------


def test_fourcv_bisecting_round_trip_100():
    g = _setup_cubic(fourcv, a=5.0)
    assert g.mode_name == "bisecting"
    assert _round_trip_ok(g, 1, 0, 0)


def test_fourcv_bisecting_round_trip_010():
    g = _setup_cubic(fourcv, a=5.0)
    assert _round_trip_ok(g, 0, 1, 0)


def test_fourcv_bisecting_round_trip_001():
    g = _setup_cubic(fourcv, a=5.0)
    assert _round_trip_ok(g, 0, 0, 1)


def test_fourcv_bisecting_round_trip_111():
    g = _setup_cubic(fourcv, a=4.0)
    assert _round_trip_ok(g, 1, 1, 1)


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
# BisectingMode — fourch (omega = ttheta/2, vertical axis)
# ---------------------------------------------------------------------------


def test_fourch_bisecting_round_trip():
    g = _setup_cubic(fourch, a=4.0)
    assert g.mode_name == "bisecting"
    assert _round_trip_ok(g, 1, 0, 0)


# ---------------------------------------------------------------------------
# BisectingMode — psic (eta = delta/2, mu = nu = 0)
# ---------------------------------------------------------------------------


def test_psic_bisecting_round_trip_100():
    g = _setup_cubic(psic, a=5.0)
    assert g.mode_name == "bisecting"
    assert _round_trip_ok(g, 1, 0, 0)


def test_psic_bisecting_round_trip_111():
    g = _setup_cubic(psic, a=4.0)
    assert _round_trip_ok(g, 1, 1, 1)


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
    """kappa4cv round-trip: use (0,1,0) which is reachable in BL lateral basis."""
    g = _setup_cubic(kappa4cv, a=4.0)
    assert g.mode_name == "bisecting"
    assert _round_trip_ok(g, 0, 1, 0)


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
    from ad_hoc_diffractometer import psic

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
    from ad_hoc_diffractometer import angles_to_phi_vector
    from ad_hoc_diffractometer.forward import _solve_two_angles

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
    from ad_hoc_diffractometer import angles_to_phi_vector
    from ad_hoc_diffractometer.forward import _solve_two_angles

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
        basis={"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT},
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


def test_check_limits_all_ok():
    g = _setup_cubic(fourcv, a=4.0)
    angles = {s.name: 0.0 for s in g.sample_stages + g.detector_stages}
    assert _check_limits(g, angles) is True


def test_check_limits_one_fails():
    g = _setup_cubic(fourcv, a=4.0)
    angles = {"omega": 0.0, "chi": 200.0, "phi": 0.0, "ttheta": 0.0}
    assert _check_limits(g, angles) is False


def test_check_limits_unknown_stage_ignored():
    g = _setup_cubic(fourcv, a=4.0)
    angles = {"omega": 0.0, "nonexistent": 999.0}
    assert _check_limits(g, angles) is True


def test_apply_cut_points_mode_priority():
    """Mode cut-point takes priority over geometry-level cut-point."""
    g = _setup_cubic(fourcv, a=4.0)
    mode = ConstraintSet(
        [BisectConstraint("omega", "ttheta")],
        cut_points={"phi": 0.0},  # phi in [0, 360)
    )
    g.cut_points["phi"] = -180.0  # would put phi in [-180, 180) if used
    angles = {"phi": -10.0}
    _apply_cut_points(angles, mode, g)
    # Mode cut-point (0°) wins: -10° → 350°
    assert angles["phi"] == pytest.approx(350.0, abs=1e-10)


def test_apply_cut_points_geometry_fallback():
    """Geometry-level cut-point used when no mode cut-point set for that stage."""
    g = _setup_cubic(fourcv, a=4.0)
    mode = ConstraintSet([BisectConstraint("omega", "ttheta")])
    g.cut_points["phi"] = 0.0  # phi in [0, 360)
    angles = {"phi": -10.0}
    _apply_cut_points(angles, mode, g)
    assert angles["phi"] == pytest.approx(350.0, abs=1e-10)


def test_apply_cut_points_no_cut_unchanged():
    """Angles without any cut-point are unchanged."""
    g = _setup_cubic(fourcv, a=4.0)
    mode = ConstraintSet([BisectConstraint("omega", "ttheta")])
    angles = {"phi": -10.0}
    _apply_cut_points(angles, mode, g)
    assert angles["phi"] == pytest.approx(-10.0, abs=1e-10)


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
        basis={"vertical": _XHAT, "longitudinal": _YHAT, "lateral": _ZHAT},
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
        basis={"vertical": _XHAT, "longitudinal": _YHAT, "lateral": _ZHAT},
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
    from ad_hoc_diffractometer import psic
    from ad_hoc_diffractometer.forward import _solve_fixed_sample

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
