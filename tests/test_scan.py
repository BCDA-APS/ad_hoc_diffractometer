# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.scan.

Covers:
  - NEAREST_ANGLES: scoring, previous=None, minimises distance
  - _hkl_points: line, radial, transverse; n_points<2; unknown type; zero direction
  - _pick_solution: None key, custom key, empty list
  - _euler_from_Z_standard: round-trip decomposition for fourcv and psic
  - _kappa_from_Z: round-trip decomposition for kappa4cv
  - hkl_trajectory: line/radial/transverse, round-trip invariant,
      inaccessible points, NEAREST_ANGLES continuity, custom solution_key,
      error cases
  - psi_trajectory: BL1967 psi round-trip on fourcv, kappa4cv, psic;
      NEAREST_ANGLES continuity; no-Bragg path; error cases
  - trajectory_plan: hkl/Q space, accessible/inaccessible flags,
      NEAREST_ANGLES continuity, space validation, n_points<2 raises
"""

import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import kappa4cv
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.scan import NEAREST_ANGLES
from ad_hoc_diffractometer.scan import _euler_from_Z_standard
from ad_hoc_diffractometer.scan import _hkl_points
from ad_hoc_diffractometer.scan import _kappa_from_Z
from ad_hoc_diffractometer.scan import _pick_solution
from ad_hoc_diffractometer.scan import hkl_trajectory
from ad_hoc_diffractometer.scan import psi_trajectory
from ad_hoc_diffractometer.scan import trajectory_plan

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

WAVELENGTH = 1.5406  # Cu Kα in Å
HKL_TEST = (1, 1, 0)  # reflection used throughout psi tests


def _setup(factory, a=4.0, *, mode="bisecting"):
    """Return a geometry with wavelength set, cubic lattice, UB=B."""
    g = factory()
    g.wavelength = WAVELENGTH
    g.sample.lattice = Lattice(a=a)
    ub_identity(g.sample)
    g.mode_name = mode
    return g


def _round_trip_ok(g, angles, hkl, atol=1e-4):
    """Return True if g.inverse(angles) ≈ hkl."""
    return np.allclose(g.inverse(angles), hkl, atol=atol)


# ---------------------------------------------------------------------------
# NEAREST_ANGLES
# ---------------------------------------------------------------------------


def test_nearest_angles_previous_none():
    """Score is 0.0 when previous is None."""
    candidate = {"omega": 10.0, "chi": 45.0, "phi": 90.0, "ttheta": 20.0}
    assert NEAREST_ANGLES(candidate, None) == 0.0


def test_nearest_angles_identical():
    """Score is 0.0 when candidate equals previous."""
    angles = {"omega": 10.0, "chi": 45.0, "phi": 90.0}
    assert NEAREST_ANGLES(angles, angles) == pytest.approx(0.0)


def test_nearest_angles_nonzero():
    """Score equals the expected sum of squares."""
    candidate = {"omega": 11.0, "chi": 46.0}
    previous = {"omega": 10.0, "chi": 45.0}
    assert NEAREST_ANGLES(candidate, previous) == pytest.approx(2.0)


def test_nearest_angles_picks_closer():
    """_pick_solution with NEAREST_ANGLES picks the nearer of two candidates."""
    previous = {"omega": 10.0, "chi": 45.0, "phi": 0.0}
    close = {"omega": 10.5, "chi": 45.5, "phi": 0.5}
    far = {"omega": 20.0, "chi": -45.0, "phi": 180.0}
    assert _pick_solution([far, close], previous, NEAREST_ANGLES) is close


# ---------------------------------------------------------------------------
# _pick_solution
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "solutions, previous, key, expected_idx, context",
    [
        pytest.param(
            [],
            None,
            NEAREST_ANGLES,
            None,
            does_not_raise(),
            id="empty-list-returns-none",
        ),
        pytest.param(
            [{"a": 1.0}, {"a": 2.0}],
            None,
            None,
            0,
            does_not_raise(),
            id="none-key-returns-first",
        ),
        pytest.param(
            [{"a": 5.0}, {"a": 1.1}],
            {"a": 1.0},
            NEAREST_ANGLES,
            1,  # {"a": 1.1} is much closer to {"a": 1.0} than {"a": 5.0}
            does_not_raise(),
            id="nearest-picks-closest",
        ),
    ],
)
def test_pick_solution(solutions, previous, key, expected_idx, context):
    with context:
        result = _pick_solution(solutions, previous, key)
        if expected_idx is None:
            assert result is None
        else:
            assert result is solutions[expected_idx]


# ---------------------------------------------------------------------------
# _hkl_points
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "trajectory, n_points, expected_first, expected_last, context",
    [
        pytest.param(
            {"type": "line", "start": (1, 0, 0), "end": (3, 0, 0)},
            5,
            (1.0, 0.0, 0.0),
            (3.0, 0.0, 0.0),
            does_not_raise(),
            id="line-endpoints",
        ),
        pytest.param(
            {
                "type": "radial",
                "center": (1, 0, 0),
                "direction": (1, 0, 0),
                "extent": 1.0,
            },
            5,
            (0.0, 0.0, 0.0),
            (2.0, 0.0, 0.0),
            does_not_raise(),
            id="radial-endpoints",
        ),
        pytest.param(
            {
                "type": "transverse",
                "center": (0, 0, 1),
                "Q_ref": (0, 0, 1),
                "extent": 0.5,
            },
            3,
            None,
            None,
            does_not_raise(),
            id="transverse-runs",
        ),
        pytest.param(
            {"type": "line", "start": (0, 0, 1), "end": (0, 0, 2)},
            1,
            None,
            None,
            pytest.raises(ValueError, match=re.escape("n_points must be at least 2")),
            id="n-points-too-small",
        ),
        pytest.param(
            {"type": "unknown"},
            3,
            None,
            None,
            pytest.raises(ValueError, match=re.escape("unknown trajectory type")),
            id="unknown-type",
        ),
        pytest.param(
            {
                "type": "radial",
                "center": (0, 0, 1),
                "direction": (0, 0, 0),
                "extent": 1.0,
            },
            3,
            None,
            None,
            pytest.raises(ValueError, match=re.escape("'direction' must be non-zero")),
            id="radial-zero-direction",
        ),
        pytest.param(
            {
                "type": "transverse",
                "center": (0, 0, 1),
                "Q_ref": (0, 0, 0),
                "extent": 0.5,
            },
            3,
            None,
            None,
            pytest.raises(ValueError, match=re.escape("'Q_ref' must be non-zero")),
            id="transverse-zero-qref",
        ),
    ],
)
def test_hkl_points(trajectory, n_points, expected_first, expected_last, context):
    with context:
        points = _hkl_points(trajectory, n_points)
        assert len(points) == n_points
        if expected_first is not None:
            assert np.allclose(points[0], expected_first, atol=1e-12)
        if expected_last is not None:
            assert np.allclose(points[-1], expected_last, atol=1e-12)


def test_hkl_points_transverse_perpendicular():
    """Transverse points lie perpendicular to Q_ref from center."""
    q_ref = np.array([0.0, 0.0, 1.0])
    center = np.array([1.0, 0.0, 0.0])
    points = _hkl_points(
        {
            "type": "transverse",
            "center": tuple(center),
            "Q_ref": tuple(q_ref),
            "extent": 0.5,
        },
        n_points=5,
    )
    for pt in points:
        disp = np.array(pt) - center
        assert abs(np.dot(disp, q_ref)) < 1e-12


# ---------------------------------------------------------------------------
# Decomposition round-trips
# ---------------------------------------------------------------------------


def test_euler_from_Z_round_trip_fourcv():
    """_euler_from_Z_standard reproduces fourcv base angles to 1e-10."""
    from ad_hoc_diffractometer.rotation import rotation_matrix as Rmat

    g = _setup(fourcv)
    base = g.forward(1, 1, 0)[0]
    for name, angle in base.items():
        try:
            g.set_angle(name, angle)
        except KeyError:
            pass
    Z0 = g.sample_rotation_matrix()
    omega_s, chi_s, phi_s = (
        g.sample_stages[-3],
        g.sample_stages[-2],
        g.sample_stages[-1],
    )
    for om, chi, phi in _euler_from_Z_standard(Z0, omega_s, chi_s, phi_s):
        Z_r = Rmat(phi_s.axis, phi) @ Rmat(chi_s.axis, chi) @ Rmat(omega_s.axis, om)
        assert np.linalg.norm(Z0 - Z_r) < 1e-10, (
            f"Recomposition error: {np.linalg.norm(Z0 - Z_r)}"
        )


def test_euler_from_Z_round_trip_psic():
    """_euler_from_Z_standard recomposes psic inner-3 Z to 1e-10."""
    from ad_hoc_diffractometer.rotation import rotation_matrix as Rmat

    g = _setup(psic)
    # Use a non-degenerate reflection (chi ≠ 0) to avoid the chi=0 degeneracy
    base = g.forward(1, 1, 0)[0]
    # Force a non-trivial chi by using the negative-chi branch
    solutions = g.forward(1, 1, 0)
    base = max(solutions, key=lambda s: abs(s.get("chi", 0)))
    for name, angle in base.items():
        try:
            g.set_angle(name, angle)
        except KeyError:
            pass
    Z0 = g.sample_rotation_matrix()
    omega_s, chi_s, phi_s = (
        g.sample_stages[-3],
        g.sample_stages[-2],
        g.sample_stages[-1],
    )
    outer = g.sample_stages[:-3]
    R_outer = np.eye(3)
    for s in outer:
        R_outer = Rmat(s.axis, base.get(s.name, s.angle)) @ R_outer
    Z_inner = R_outer.T @ Z0
    recompose_errors = []
    for om, chi, phi in _euler_from_Z_standard(Z_inner, omega_s, chi_s, phi_s):
        Z_r = Rmat(phi_s.axis, phi) @ Rmat(chi_s.axis, chi) @ Rmat(omega_s.axis, om)
        recompose_errors.append(np.linalg.norm(Z_inner - Z_r))
    assert min(recompose_errors) < 1e-10, f"Recomposition errors: {recompose_errors}"


def test_kappa_from_Z_round_trip_kappa4cv():
    """_kappa_from_Z reproduces kappa4cv base angles to 1e-10."""
    from ad_hoc_diffractometer.rotation import rotation_matrix as Rmat

    g = _setup(kappa4cv)
    base = g.forward(1, 1, 0)[0]
    for name, angle in base.items():
        try:
            g.set_angle(name, angle)
        except KeyError:
            pass
    Z0 = g.sample_rotation_matrix()
    kom_s, kap_s, kph_s = g.sample_stages[-3], g.sample_stages[-2], g.sample_stages[-1]
    sols = _kappa_from_Z(Z0, kom_s, kap_s, kph_s, g.kappa_alpha_deg)
    # At least one solution must reproduce Z0 and match base angles
    assert any(
        np.linalg.norm(
            Rmat(kph_s.axis, kph) @ Rmat(kap_s.axis, kap) @ Rmat(kom_s.axis, kom) - Z0
        )
        < 1e-10
        for kom, kap, kph in sols
    )
    # One solution must match the base motor angles
    base_match = any(
        abs(kom - base["komega"]) < 1e-4
        and abs(kap - base["kappa"]) < 1e-4
        and abs(kph - base["kphi"]) < 1e-4
        for kom, kap, kph in sols
    )
    assert base_match


# ---------------------------------------------------------------------------
# hkl_trajectory
# ---------------------------------------------------------------------------


class TestHklTrajectory:
    def test_line_round_trip(self):
        """Every accessible point satisfies inverse(angles) ≈ hkl."""
        g = _setup(fourcv)
        result = hkl_trajectory(
            g,
            {"type": "line", "start": (1, 0, 0), "end": (2, 0, 0)},
            n_points=5,
        )
        assert len(result) == 5
        for pt in result:
            if pt["angles"] is not None:
                assert _round_trip_ok(g, pt["angles"], pt["hkl"])
                assert pt["warning"] is None

    def test_radial_trajectory_endpoints(self):
        """Radial trajectory: endpoints match expected hkl."""
        g = _setup(fourcv)
        result = hkl_trajectory(
            g,
            {
                "type": "radial",
                "center": (1, 0, 0),
                "direction": (1, 0, 0),
                "extent": 0.5,
            },
            n_points=3,
        )
        assert len(result) == 3
        assert np.allclose(result[0]["hkl"], (0.5, 0.0, 0.0), atol=1e-12)
        assert np.allclose(result[2]["hkl"], (1.5, 0.0, 0.0), atol=1e-12)

    def test_transverse_perpendicular(self):
        """Transverse points lie perpendicular to Q_ref from center."""
        g = _setup(fourcv)
        result = hkl_trajectory(
            g,
            {
                "type": "transverse",
                "center": (1, 0, 0),
                "Q_ref": (1, 0, 0),
                "extent": 0.1,
            },
            n_points=5,
        )
        center = np.array([1.0, 0.0, 0.0])
        q_ref = np.array([1.0, 0.0, 0.0])
        for pt in result:
            disp = np.array(pt["hkl"]) - center
            assert abs(np.dot(disp, q_ref)) < 1e-12

    def test_inaccessible_gives_none_and_warning(self):
        """When Q exceeds Ewald sphere, angles=None and warning is non-empty."""
        g = _setup(fourcv, a=0.1)  # tiny lattice → very large |Q|
        result = hkl_trajectory(
            g,
            {"type": "line", "start": (1, 0, 0), "end": (2, 0, 0)},
            n_points=3,
        )
        inaccessible = [pt for pt in result if pt["angles"] is None]
        assert len(inaccessible) > 0
        for pt in inaccessible:
            assert pt["warning"] is not None

    def test_nearest_angles_prevents_branch_flip(self):
        """chi values stay on the same branch (no sign flip) with NEAREST_ANGLES."""
        g = _setup(fourcv)
        result = hkl_trajectory(
            g,
            {"type": "line", "start": (0.5, 0, 0), "end": (1.5, 0, 0)},
            n_points=7,
            solution_key=NEAREST_ANGLES,
        )
        chi_values = [pt["angles"]["chi"] for pt in result if pt["angles"] is not None]
        assert len(chi_values) >= 2
        signs = {1 if v >= 0 else -1 for v in chi_values}
        assert len(signs) == 1, f"Branch flip detected: chi = {chi_values}"

    def test_custom_solution_key(self):
        """A caller-supplied solution_key controls solution selection."""
        g = _setup(fourcv)
        # Always pick the solution with the largest (most positive) chi
        max_chi_key = lambda c, p: -c.get("chi", 0.0)  # noqa: E731
        result = hkl_trajectory(
            g,
            {"type": "line", "start": (1, 0, 0), "end": (1, 0, 0)},
            n_points=2,
            solution_key=max_chi_key,
        )
        for pt in result:
            if pt["angles"] is not None:
                assert pt["angles"]["chi"] >= 0.0

    def test_none_solution_key_returns_result(self):
        """solution_key=None runs and returns angles dicts."""
        g = _setup(fourcv)
        result = hkl_trajectory(
            g,
            {"type": "line", "start": (1, 0, 0), "end": (1, 0, 0)},
            n_points=2,
            solution_key=None,
        )
        assert all("angles" in pt for pt in result)

    @pytest.mark.parametrize(
        "trajectory, n_points, context",
        [
            pytest.param(
                {"type": "line", "start": (1, 0, 0), "end": (2, 0, 0)},
                1,
                pytest.raises(
                    ValueError, match=re.escape("n_points must be at least 2")
                ),
                id="n-points-too-small",
            ),
            pytest.param(
                {"type": "bad_type"},
                3,
                pytest.raises(ValueError, match=re.escape("unknown trajectory type")),
                id="unknown-trajectory-type",
            ),
        ],
    )
    def test_hkl_trajectory_raises(self, trajectory, n_points, context):
        g = _setup(fourcv)
        with context:
            hkl_trajectory(g, trajectory, n_points)

    def test_no_wavelength_gives_warnings(self):
        """When wavelength is not set all points have warnings and angles=None."""
        g = fourcv()
        g.sample.lattice = Lattice(a=4.0)
        ub_identity(g.sample)
        g.mode_name = "bisecting"
        result = hkl_trajectory(
            g, {"type": "line", "start": (1, 0, 0), "end": (2, 0, 0)}, n_points=3
        )
        assert all(pt["angles"] is None for pt in result)
        assert all(pt["warning"] is not None for pt in result)


# ---------------------------------------------------------------------------
# psi_trajectory
# ---------------------------------------------------------------------------


class TestPsiTrajectory:
    """
    psi_trajectory() uses the BL1967 operational ψ definition.
    ψ = 0 at the base forward() solution; each psi_target specifies the
    rotation of the sample about Q relative to that reference.
    """

    PSI_TARGETS = list(range(-90, 91, 30))

    @staticmethod
    def _psi_err(actual, target):
        """Angular error mod 360, in (-180, 180]."""
        diff = (actual - target + 180.0) % 360.0 - 180.0
        return abs(diff)

    def _psi_round_trip(self, geometry, hkl=HKL_TEST, targets=None, atol=0.05):
        """Assert psi_actual ≈ psi_target (mod 360) for all accessible points."""
        if targets is None:
            targets = self.PSI_TARGETS
        result = psi_trajectory(geometry, *hkl, targets)
        assert len(result) == len(targets)
        accessible = [pt for pt in result if pt["angles"] is not None]
        assert len(accessible) >= len(targets) // 2, (
            f"Too few accessible psi points: {[(pt['psi_target'], pt['warning']) for pt in result]}"
        )
        for pt in accessible:
            err = self._psi_err(pt["psi_actual"], pt["psi_target"])
            assert err <= atol, (
                f"psi mismatch: target={pt['psi_target']}, actual={pt['psi_actual']}, err={err}"
            )
            assert pt["warning"] is None

    def test_psi_round_trip_fourcv(self):
        """BL1967 psi round-trip on fourcv."""
        g = _setup(fourcv)
        self._psi_round_trip(g)

    def test_psi_round_trip_psic(self):
        """BL1967 psi round-trip on psic using (0,1,1) which has non-zero chi."""
        g = _setup(psic)
        # (1,1,0) gives chi≈0 which is a degenerate case for Euler decomposition;
        # use (0,1,1) instead which has chi ≈ ±75° (well-conditioned).
        self._psi_round_trip(g, hkl=(0, 1, 1), targets=list(range(-60, 61, 30)))

    def test_psi_round_trip_kappa4cv(self):
        """BL1967 psi round-trip on kappa4cv (mid-range, avoids arm limits)."""
        g = _setup(kappa4cv)
        self._psi_round_trip(g, targets=list(range(-60, 61, 30)), atol=0.1)

    def test_psi_zero_at_base(self):
        """psi_actual = 0 when psi_target = 0 (base forward() solution)."""
        g = _setup(fourcv)
        result = psi_trajectory(g, *HKL_TEST, [0.0])
        assert result[0]["psi_actual"] == pytest.approx(0.0, abs=1e-6)

    def test_psi_smooth_motion_nearest_angles(self):
        """phi values vary smoothly across a dense ψ sweep (no large jumps)."""
        g = _setup(fourcv)
        targets = list(range(-60, 61, 10))
        result = psi_trajectory(g, *HKL_TEST, targets, solution_key=NEAREST_ANGLES)
        phi_values = [pt["angles"]["phi"] for pt in result if pt["angles"] is not None]
        assert len(phi_values) >= 2
        for i in range(1, len(phi_values)):
            diff = abs(phi_values[i] - phi_values[i - 1])
            assert diff < 45.0, (
                f"Large phi jump at index {i}: "
                f"{phi_values[i - 1]:.1f}° → {phi_values[i]:.1f}°"
            )

    def test_psi_no_bragg_solution(self):
        """When reflection is inaccessible all psi points have warning."""
        g = _setup(fourcv)
        g.stage("ttheta").limits = (0.0, 1.0)  # block the detector
        result = psi_trajectory(g, *HKL_TEST, [0.0, 30.0])
        for pt in result:
            assert pt["angles"] is None
            assert pt["warning"] is not None

    def test_psi_no_wavelength_raises(self):
        """Raises ValueError when wavelength is not set."""
        g = fourcv()
        g.sample.lattice = Lattice(a=4.0)
        ub_identity(g.sample)
        g.mode_name = "bisecting"
        with pytest.raises(ValueError, match=re.escape("wavelength")):
            psi_trajectory(g, *HKL_TEST, [0.0])

    def test_psi_no_ub_raises(self):
        """Raises ValueError when UB is not set."""
        g = fourcv()
        g.wavelength = WAVELENGTH
        g.mode_name = "bisecting"
        with pytest.raises(ValueError, match=re.escape("UB")):
            psi_trajectory(g, *HKL_TEST, [0.0])

    def test_psi_zero_hkl_raises(self):
        """Raises ValueError for hkl = (0,0,0)."""
        g = _setup(fourcv)
        with pytest.raises(ValueError, match=re.escape("(0, 0, 0)")):
            psi_trajectory(g, 0, 0, 0, [0.0])

    def test_psi_result_keys(self):
        """Every result entry has exactly the expected keys."""
        g = _setup(fourcv)
        result = psi_trajectory(g, *HKL_TEST, [0.0])
        assert set(result[0].keys()) == {
            "psi_target",
            "psi_actual",
            "angles",
            "warning",
        }

    def test_psi_none_solution_key(self):
        """solution_key=None runs without error."""
        g = _setup(fourcv)
        result = psi_trajectory(g, *HKL_TEST, [0.0], solution_key=None)
        assert len(result) == 1


# ---------------------------------------------------------------------------
# trajectory_plan
# ---------------------------------------------------------------------------


class TestTrajectoryPlan:
    def test_hkl_space_endpoints(self):
        """Start and end hkl are exactly the requested values."""
        g = _setup(fourcv)
        plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=5)
        assert len(plan) == 5
        assert np.allclose(plan[0]["hkl"], (1, 0, 0), atol=1e-12)
        assert np.allclose(plan[-1]["hkl"], (2, 0, 0), atol=1e-12)

    def test_hkl_space_round_trip(self):
        """All accessible points satisfy inverse(angles) ≈ hkl."""
        g = _setup(fourcv)
        plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=5)
        for pt in plan:
            if pt["accessible"]:
                assert _round_trip_ok(g, pt["angles"], pt["hkl"])
                assert pt["warnings"] == []

    def test_q_space_endpoints_exact(self):
        """space='Q': endpoints are exact despite Q interpolation."""
        g = _setup(fourcv)
        plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=5, space="Q")
        assert np.allclose(plan[0]["hkl"], (1, 0, 0), atol=1e-12)
        assert np.allclose(plan[-1]["hkl"], (2, 0, 0), atol=1e-12)

    def test_q_space_equal_dq_steps(self):
        """space='Q': consecutive Q-vector differences are equal in magnitude."""
        g = _setup(fourcv)
        # Tetragonal lattice: a≠c, so hkl-space is not uniform in Q-space.
        g.sample.lattice = ahd.Lattice(a=3.0, c=6.0)
        ub_identity(g.sample)

        n = 5
        plan_q = trajectory_plan(g, (1, 0, 0), (0, 0, 2), n_points=n, space="Q")

        # Compute Q vectors at each trajectory point
        Q_vecs = [g.sample.UB @ np.array(pt["hkl"]) for pt in plan_q]
        # Consecutive differences in Q should all be equal vectors
        dQ = [Q_vecs[i + 1] - Q_vecs[i] for i in range(n - 1)]
        dQ_mags = [np.linalg.norm(d) for d in dQ]
        # All step sizes equal (Q-space interpolation is uniform)
        assert np.std(dQ_mags) < 1e-10, f"Q step sizes not equal: {dQ_mags}"

        # Verify both plans share the same exact endpoints regardless of space.
        plan_hkl = trajectory_plan(g, (1, 0, 0), (0, 0, 2), n_points=n, space="hkl")
        assert np.allclose(plan_q[0]["hkl"], (1, 0, 0), atol=1e-12)
        assert np.allclose(plan_q[-1]["hkl"], (0, 0, 2), atol=1e-12)
        assert np.allclose(plan_hkl[0]["hkl"], (1, 0, 0), atol=1e-12)
        assert np.allclose(plan_hkl[-1]["hkl"], (0, 0, 2), atol=1e-12)

    def test_inaccessible_points_flagged(self):
        """Points with no valid solution are accessible=False."""
        g = _setup(fourcv, a=0.1)
        plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=3)
        inaccessible = [pt for pt in plan if not pt["accessible"]]
        assert len(inaccessible) > 0
        for pt in inaccessible:
            assert pt["angles"] is None
            assert len(pt["warnings"]) > 0

    def test_nearest_angles_no_branch_flip(self):
        """chi values stay on the same branch across the plan."""
        g = _setup(fourcv)
        plan = trajectory_plan(g, (0.5, 0, 0), (1.5, 0, 0), n_points=9)
        chi_values = [pt["angles"]["chi"] for pt in plan if pt["accessible"]]
        assert len(chi_values) >= 2
        signs = {1 if v >= 0 else -1 for v in chi_values}
        assert len(signs) == 1, f"Branch flip: chi = {chi_values}"

    def test_n_points_too_small_raises(self):
        g = _setup(fourcv)
        with pytest.raises(ValueError, match=re.escape("n_points must be at least 2")):
            trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=1)

    def test_invalid_space_raises(self):
        g = _setup(fourcv)
        with pytest.raises(ValueError, match=re.escape("space must be 'hkl' or 'Q'")):
            trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=3, space="crystal")

    def test_q_space_no_ub_raises(self):
        g = fourcv()
        g.wavelength = WAVELENGTH
        g.mode_name = "bisecting"
        with pytest.raises(ValueError, match=re.escape("UB")):
            trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=3, space="Q")

    def test_result_keys(self):
        g = _setup(fourcv)
        plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=3)
        for pt in plan:
            assert set(pt.keys()) == {"hkl", "angles", "accessible", "warnings"}

    def test_n_points_respected(self):
        g = _setup(fourcv)
        for n in [2, 5, 11]:
            plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=n)
            assert len(plan) == n

    def test_psic_accessible(self):
        """psic geometry: most points on a simple trajectory should be accessible."""
        g = _setup(psic)
        plan = trajectory_plan(g, (0, 0, 1), (0, 0, 2), n_points=5)
        assert sum(pt["accessible"] for pt in plan) >= 3

    def test_no_wavelength_gives_inaccessible(self):
        """When wavelength is not set all points are inaccessible."""
        g = fourcv()
        g.sample.lattice = Lattice(a=4.0)
        ub_identity(g.sample)
        g.mode_name = "bisecting"
        plan = trajectory_plan(g, (1, 0, 0), (2, 0, 0), n_points=3)
        assert all(not pt["accessible"] for pt in plan)
        assert all(len(pt["warnings"]) > 0 for pt in plan)


# ---------------------------------------------------------------------------
# Public API surface
# ---------------------------------------------------------------------------


def test_public_api_exports():
    """All scan symbols are importable from the top-level package."""
    assert ahd.NEAREST_ANGLES is NEAREST_ANGLES
    assert ahd.hkl_trajectory is hkl_trajectory
    assert ahd.psi_trajectory is psi_trajectory
    assert ahd.trajectory_plan is trajectory_plan
