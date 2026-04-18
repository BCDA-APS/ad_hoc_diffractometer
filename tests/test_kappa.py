# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.kappa — kappa-to-Eulerian conversion.

Covers:
  - kappa_to_eulerian: forward direction, both branches, singularity
  - eulerian_to_kappa: inverse direction, both branches, singularity
  - Round-trip: eulerian -> kappa -> eulerian
  - Round-trip: kappa -> eulerian -> kappa
  - Branch selection: +1 gives kappa >= 0, -1 gives kappa <= 0
  - Invalid branch raises ValueError
  - Chi outside reachable range raises ValueError
  - kappa_to_eulerian unreachable chi raises ValueError
  - Default alpha_deg = 50.0
  - Non-default alpha_deg
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import pytest

from ad_hoc_diffractometer import eulerian_to_kappa
from ad_hoc_diffractometer import kappa_to_eulerian

# ---------------------------------------------------------------------------
# Known values from Walko (2016) eq. [16] at alpha_0 = 50°
# At chi=90, omega=0, phi=0:
#   kappa = 2*arcsin(sin(45)/sin(50)) = 134.7559...
#   offset = arccos(cos(kappa/2)/cos(45)) = 57.0452...
#   komega = kphi = 0 + offset = 57.0452...
# ---------------------------------------------------------------------------

ALPHA = 50.0
_k90 = 2.0 * math.degrees(
    math.asin(math.sin(math.radians(45.0)) / math.sin(math.radians(ALPHA)))
)
_off90 = math.degrees(
    math.acos(math.cos(math.radians(_k90 / 2.0)) / math.cos(math.radians(45.0)))
)
KNOWN_CHI90 = (
    _off90,  # komega
    _k90,  # kappa
    _off90,  # kphi
)


# ---------------------------------------------------------------------------
# kappa_to_eulerian — forward direction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "komega, kappa, kphi, alpha_deg, expected_omega, expected_chi, expected_phi",
    [
        # chi=0: kappa=0, offset=0
        pytest.param(0.0, 0.0, 0.0, 50.0, 0.0, 0.0, 0.0, id="chi0-all-zero"),
        # chi=90: known values
        pytest.param(
            KNOWN_CHI90[0],
            KNOWN_CHI90[1],
            KNOWN_CHI90[2],
            50.0,
            0.0,
            90.0,
            0.0,
            id="chi90-omega0-phi0",
        ),
        # Negative kappa branch: negate komega, kappa, kphi
        pytest.param(
            -KNOWN_CHI90[0],
            -KNOWN_CHI90[1],
            -KNOWN_CHI90[2],
            50.0,
            0.0,
            -90.0,
            0.0,
            id="chi90-negative-branch",
        ),
        # Non-zero omega and phi
        pytest.param(
            KNOWN_CHI90[0] + 30.0,
            KNOWN_CHI90[1],
            KNOWN_CHI90[2] + 15.0,
            50.0,
            30.0,
            90.0,
            15.0,
            id="chi90-nonzero-omega-phi",
        ),
        # Different alpha
        pytest.param(0.0, 0.0, 0.0, 60.0, 0.0, 0.0, 0.0, id="alpha60-zero"),
    ],
)
def test_kappa_to_eulerian(
    komega,
    kappa,
    kphi,
    alpha_deg,
    expected_omega,
    expected_chi,
    expected_phi,
):
    omega, chi, phi = kappa_to_eulerian(komega, kappa, kphi, alpha_deg=alpha_deg)
    assert omega == pytest.approx(expected_omega, abs=1e-8)
    assert chi == pytest.approx(expected_chi, abs=1e-8)
    assert phi == pytest.approx(expected_phi, abs=1e-8)


# ---------------------------------------------------------------------------
# eulerian_to_kappa — inverse direction
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "omega, chi, phi, alpha_deg, branch, expected_komega, expected_kappa, expected_kphi",
    [
        # chi=0: kappa=0, offset=0
        pytest.param(0.0, 0.0, 0.0, 50.0, +1, 0.0, 0.0, 0.0, id="chi0-branch+1"),
        pytest.param(0.0, 0.0, 0.0, 50.0, -1, 0.0, -0.0, 0.0, id="chi0-branch-1"),
        # chi=90, omega=0, phi=0 — positive branch
        pytest.param(
            0.0,
            90.0,
            0.0,
            50.0,
            +1,
            KNOWN_CHI90[0],
            KNOWN_CHI90[1],
            KNOWN_CHI90[2],
            id="chi90-branch+1",
        ),
        # chi=90, omega=0, phi=0 — negative branch
        pytest.param(
            0.0,
            90.0,
            0.0,
            50.0,
            -1,
            -KNOWN_CHI90[0],
            -KNOWN_CHI90[1],
            -KNOWN_CHI90[2],
            id="chi90-branch-1",
        ),
        # Non-zero omega and phi
        pytest.param(
            30.0,
            90.0,
            15.0,
            50.0,
            +1,
            KNOWN_CHI90[0] + 30.0,
            KNOWN_CHI90[1],
            KNOWN_CHI90[2] + 15.0,
            id="chi90-nonzero-omega-phi",
        ),
    ],
)
def test_eulerian_to_kappa(
    omega,
    chi,
    phi,
    alpha_deg,
    branch,
    expected_komega,
    expected_kappa,
    expected_kphi,
):
    komega, kappa, kphi = eulerian_to_kappa(
        omega, chi, phi, alpha_deg=alpha_deg, branch=branch
    )
    assert komega == pytest.approx(expected_komega, abs=1e-8)
    assert abs(kappa) == pytest.approx(abs(expected_kappa), abs=1e-8)
    assert kphi == pytest.approx(expected_kphi, abs=1e-8)


# ---------------------------------------------------------------------------
# Round-trips
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "omega, chi, phi, alpha_deg, branch",
    [
        pytest.param(0.0, 0.0, 0.0, 50.0, +1, id="zero-branch+1"),
        # branch=-1 gives negative kappa which maps to negative chi
        # round-trip: eulerian(+chi) -> kappa(neg) -> eulerian(-chi) != +chi
        # So round-trips with branch=-1 use negative chi as input
        pytest.param(0.0, 0.0, 0.0, 50.0, -1, id="zero-branch-1"),
        pytest.param(0.0, 90.0, 0.0, 50.0, +1, id="chi90-branch+1"),
        pytest.param(0.0, -90.0, 0.0, 50.0, -1, id="neg-chi90-branch-1"),
        pytest.param(30.0, 45.0, -20.0, 50.0, +1, id="general-branch+1"),
        pytest.param(30.0, -45.0, -20.0, 50.0, -1, id="general-neg-chi-branch-1"),
        pytest.param(0.0, 60.0, 0.0, 50.0, +1, id="pos-chi-branch+1"),
        pytest.param(0.0, -60.0, 0.0, 50.0, -1, id="neg-chi-branch-1"),
        pytest.param(10.0, 30.0, 15.0, 60.0, +1, id="alpha60-branch+1"),
        pytest.param(10.0, -30.0, -15.0, 45.0, -1, id="alpha45-neg-chi-branch-1"),
    ],
)
def test_round_trip_eulerian_kappa_eulerian(omega, chi, phi, alpha_deg, branch):
    """eulerian -> kappa -> eulerian recovers original (omega, chi, phi)."""
    komega, kappa, kphi = eulerian_to_kappa(
        omega, chi, phi, alpha_deg=alpha_deg, branch=branch
    )
    o2, c2, p2 = kappa_to_eulerian(komega, kappa, kphi, alpha_deg=alpha_deg)
    assert o2 == pytest.approx(omega, abs=1e-10)
    assert c2 == pytest.approx(chi, abs=1e-10)
    assert p2 == pytest.approx(phi, abs=1e-10)


@pytest.mark.parametrize(
    "komega, kappa, kphi, alpha_deg",
    [
        pytest.param(0.0, 0.0, 0.0, 50.0, id="zero"),
        pytest.param(57.045, 134.756, 57.045, 50.0, id="chi90-pos"),
        pytest.param(-57.045, -134.756, -57.045, 50.0, id="chi90-neg"),
        pytest.param(20.0, 80.0, 10.0, 50.0, id="general-pos"),
        pytest.param(-20.0, -80.0, -10.0, 50.0, id="general-neg-kappa"),
        pytest.param(0.0, 60.0, 0.0, 60.0, id="alpha60"),
    ],
)
def test_round_trip_kappa_eulerian_kappa(komega, kappa, kphi, alpha_deg):
    """kappa -> eulerian -> kappa recovers original (komega, kappa, kphi)."""
    branch = +1 if kappa >= 0 else -1
    omega, chi, phi = kappa_to_eulerian(komega, kappa, kphi, alpha_deg=alpha_deg)
    ko2, k2, kp2 = eulerian_to_kappa(
        omega, chi, phi, alpha_deg=alpha_deg, branch=branch
    )
    assert ko2 == pytest.approx(komega, abs=1e-8)
    assert k2 == pytest.approx(kappa, abs=1e-8)
    assert kp2 == pytest.approx(kphi, abs=1e-8)


# ---------------------------------------------------------------------------
# Branch selection
# ---------------------------------------------------------------------------


def test_branch_plus1_gives_positive_kappa():
    """branch=+1 gives kappa >= 0."""
    _, kappa, _ = eulerian_to_kappa(0.0, 60.0, 0.0, branch=+1)
    assert kappa >= 0.0


def test_branch_minus1_gives_negative_kappa():
    """branch=-1 gives kappa <= 0."""
    _, kappa, _ = eulerian_to_kappa(0.0, 60.0, 0.0, branch=-1)
    assert kappa <= 0.0


def test_both_branches_give_opposite_chi():
    """Positive branch gives +chi; negative branch gives -chi."""
    _, k1, _ = eulerian_to_kappa(0.0, 45.0, 0.0, branch=+1)
    _, k2, _ = eulerian_to_kappa(0.0, 45.0, 0.0, branch=-1)
    _, c1, _ = kappa_to_eulerian(0.0, k1, 0.0)
    _, c2, _ = kappa_to_eulerian(0.0, k2, 0.0)
    assert c1 == pytest.approx(45.0, abs=1e-10)
    assert c2 == pytest.approx(-45.0, abs=1e-10)


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "branch, context",
    [
        pytest.param(+1, does_not_raise(), id="branch+1-valid"),
        pytest.param(-1, does_not_raise(), id="branch-1-valid"),
        pytest.param(
            0,
            pytest.raises(ValueError, match=re.escape("branch must be +1 or -1")),
            id="branch0-raises",
        ),
        pytest.param(
            2,
            pytest.raises(ValueError, match=re.escape("branch must be +1 or -1")),
            id="branch2-raises",
        ),
    ],
)
def test_eulerian_to_kappa_invalid_branch(branch, context):
    with context:
        eulerian_to_kappa(0.0, 45.0, 0.0, branch=branch)


@pytest.mark.parametrize(
    "chi, alpha_deg, context",
    [
        pytest.param(
            0.0,
            50.0,
            does_not_raise(),
            id="chi0-ok",
        ),
        pytest.param(
            89.0,
            50.0,
            does_not_raise(),
            id="chi89-ok",
        ),
        pytest.param(
            100.0,
            50.0,
            pytest.raises(ValueError, match=re.escape("kappa geometry limit")),
            id="chi100-at-limit-raises",
        ),
        pytest.param(
            120.0,
            50.0,
            pytest.raises(ValueError, match=re.escape("kappa geometry limit")),
            id="chi120-beyond-limit-raises",
        ),
        pytest.param(
            -100.0,
            50.0,
            pytest.raises(ValueError, match=re.escape("kappa geometry limit")),
            id="neg-chi100-at-limit-raises",
        ),
    ],
)
def test_eulerian_to_kappa_chi_limit(chi, alpha_deg, context):
    """chi ≥ 2·alpha_deg raises ValueError."""
    with context:
        eulerian_to_kappa(0.0, chi, 0.0, alpha_deg=alpha_deg)


def test_kappa_to_eulerian_chi_half_near_90():
    """kappa=180 is the physical maximum; chi approaches 2*alpha_0."""
    # kappa=180, alpha=50 -> chi = 2*arcsin(sin(90)*sin(50)) = 2*arcsin(0.766) = ~100
    # This is ABOVE the chi limit (2*50=100), but kappa_to_eulerian accepts it
    # because it's the forward direction (real -> virtual).
    omega, chi, phi = kappa_to_eulerian(0.0, 180.0, 0.0, alpha_deg=50.0)
    assert abs(chi) == pytest.approx(2 * 50.0, abs=0.1)  # chi ≈ 100°


# ---------------------------------------------------------------------------
# Default alpha_deg
# ---------------------------------------------------------------------------


def test_default_alpha_is_50():
    """Default alpha_deg=50 matches explicit alpha_deg=50."""
    o1, c1, p1 = kappa_to_eulerian(10.0, 20.0, 10.0)
    o2, c2, p2 = kappa_to_eulerian(10.0, 20.0, 10.0, alpha_deg=50.0)
    assert o1 == pytest.approx(o2)
    assert c1 == pytest.approx(c2)
    assert p1 == pytest.approx(p2)


def test_default_alpha_inverse():
    """eulerian_to_kappa default alpha_deg=50 matches explicit."""
    k1 = eulerian_to_kappa(10.0, 45.0, 5.0)
    k2 = eulerian_to_kappa(10.0, 45.0, 5.0, alpha_deg=50.0)
    for v1, v2 in zip(k1, k2, strict=False):
        assert v1 == pytest.approx(v2)


# ---------------------------------------------------------------------------
# Issue #174 — kappa virtual-angle mode helpers and coverage
# ---------------------------------------------------------------------------


def test_is_kappa_virtual_mode_no_kappa_alpha():
    """is_kappa_virtual_mode returns False when geometry has no kappa_alpha_deg."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import SampleConstraint
    from ad_hoc_diffractometer.kappa import is_kappa_virtual_mode

    g = ahd.fourcv()  # no kappa_alpha_deg
    mode = ConstraintSet([SampleConstraint("omega", 0.0)])
    assert is_kappa_virtual_mode(g, mode) is False


def test_is_kappa_virtual_mode_kappa_alpha_but_no_kappa_stage():
    """is_kappa_virtual_mode returns False when kappa_alpha_deg is set but no 'kappa' stage."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import SampleConstraint
    from ad_hoc_diffractometer.kappa import is_kappa_virtual_mode

    g = ahd.fourcv()
    g._kappa_alpha_deg = 50.0  # noqa: SLF001  — set kappa_alpha_deg manually
    mode = ConstraintSet([SampleConstraint("omega", 0.0)])
    # Has kappa_alpha_deg but no stage named "kappa" → False (line 247)
    assert is_kappa_virtual_mode(g, mode) is False


def test_is_kappa_virtual_mode_no_virtual_constraint():
    """is_kappa_virtual_mode returns False when mode has no virtual angle."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import SampleConstraint
    from ad_hoc_diffractometer.kappa import is_kappa_virtual_mode

    g = ahd.kappa4cv()
    # fixed_kphi uses real stage name — not virtual
    mode = ConstraintSet([SampleConstraint("kphi", 0.0)])
    assert is_kappa_virtual_mode(g, mode) is False


def test_solve_kappa_virtual_no_kappa_stage_returns_empty():
    """solve_kappa_virtual returns [] when geometry has no kappa stage."""
    import math

    import numpy as np

    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import ConstraintSet
    from ad_hoc_diffractometer import SampleConstraint
    from ad_hoc_diffractometer.kappa import solve_kappa_virtual

    # Build a fake geometry with no "kappa" stage
    g = ahd.fourcv()
    g.wavelength = 1.5406
    g.sample.lattice = ahd.Lattice(a=4.0)
    ahd.ub_identity(g.sample)
    # Override kappa_alpha_deg to non-None so the function proceeds
    g._kappa_alpha_deg = 50.0  # noqa: SLF001

    mode = ConstraintSet([SampleConstraint("omega", 0.0)])
    Q_phi = g.sample.UB @ np.array([0.0, 1.0, 0.0])
    ttheta = 2 * math.degrees(
        math.asin(float(np.linalg.norm(Q_phi)) * g.wavelength / (4 * math.pi))
    )
    result = solve_kappa_virtual(g, Q_phi, ttheta, mode)
    # No "kappa" stage found → returns []
    assert result == []


def test_sample_constraint_is_implemented_kappa_geometry():
    """SampleConstraint with virtual name returns True on kappa geometry."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import SampleConstraint

    g = ahd.kappa4cv()
    for vname in ("omega", "chi", "phi"):
        assert SampleConstraint(vname, 0.0).is_implemented(g) is True


def test_sample_constraint_is_implemented_virtual_on_non_kappa():
    """Virtual kappa names not in a non-kappa geometry return False if not real stages."""
    import ad_hoc_diffractometer as ahd
    from ad_hoc_diffractometer import SampleConstraint

    # psic has mu, eta, chi, phi — 'omega' is NOT a real stage, and psic
    # has no kappa_alpha_deg, so SampleConstraint('omega') returns False.
    g = ahd.psic()
    assert SampleConstraint("omega", 0.0).is_implemented(g) is False
