# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.conversions.

Covers:
  - hkl_to_Q: basic, no-UB raises, (0,0,0) gives zero vector
  - Q_to_hkl: round-trip with hkl_to_Q, no-UB raises, singular-UB raises
  - Q_to_d: basic, zero-Q raises
  - d_to_Q_mag: basic, zero/negative raises, round-trip with Q_to_d
  - hkl_to_d: cubic round-trip, (0,0,0) raises, no-UB raises
  - d_to_two_theta: Bragg round-trips, bad-d raises, bad-wl raises, unreachable raises
  - two_theta_to_d: round-trip with d_to_two_theta, bad-angle raises, bad-wl raises
  - hkl_to_two_theta: cubic round-trip via known value, propagates raises
  - two_theta_to_Q_mag: basic, bad-angle raises, bad-wl raises
  - Q_mag_to_two_theta: round-trip with two_theta_to_Q_mag, unreachable raises
  - Public API exports
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import numpy as np
import pytest
from helpers import fourcv

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import Lattice
from ad_hoc_diffractometer import ub_identity
from ad_hoc_diffractometer.conversions import Q_mag_to_two_theta
from ad_hoc_diffractometer.conversions import Q_to_d
from ad_hoc_diffractometer.conversions import Q_to_hkl
from ad_hoc_diffractometer.conversions import d_to_Q_mag
from ad_hoc_diffractometer.conversions import d_to_two_theta
from ad_hoc_diffractometer.conversions import hkl_to_d
from ad_hoc_diffractometer.conversions import hkl_to_Q
from ad_hoc_diffractometer.conversions import hkl_to_two_theta
from ad_hoc_diffractometer.conversions import two_theta_to_d
from ad_hoc_diffractometer.conversions import two_theta_to_Q_mag

# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

WAVELENGTH = 1.5406  # Cu Kα in Å
A_CUBIC = 3.905  # SrTiO₃ lattice parameter (Å)
TWO_PI = 2.0 * math.pi


def _sample(a=A_CUBIC):
    """Return a Sample with a cubic UB=B matrix."""
    g = fourcv()
    g.sample.lattice = Lattice(a=a)
    ub_identity(g.sample)
    return g.sample


# ---------------------------------------------------------------------------
# hkl_to_Q
# ---------------------------------------------------------------------------


class TestHklToQ:
    def test_cubic_100(self):
        """(1,0,0) gives Q along x with magnitude 2π/a."""
        s = _sample()
        Q = hkl_to_Q(s, 1, 0, 0)
        expected = TWO_PI / A_CUBIC
        assert Q == pytest.approx([expected, 0.0, 0.0], abs=1e-10)

    def test_cubic_110(self):
        """(1,1,0) gives |Q| = 2π/a * sqrt(2)."""
        s = _sample()
        Q = hkl_to_Q(s, 1, 1, 0)
        assert float(np.linalg.norm(Q)) == pytest.approx(
            TWO_PI / A_CUBIC * math.sqrt(2), rel=1e-10
        )

    def test_zero_hkl_gives_zero_vector(self):
        """(0,0,0) maps to the zero vector."""
        s = _sample()
        Q = hkl_to_Q(s, 0, 0, 0)
        assert np.allclose(Q, [0.0, 0.0, 0.0])

    def test_no_ub_raises(self):
        """Raises ValueError when UB is not set."""
        g = fourcv()
        g.sample.lattice = Lattice(a=A_CUBIC)
        with pytest.raises(ValueError, match=re.escape("sample.UB")):
            hkl_to_Q(g.sample, 1, 0, 0)

    def test_returns_ndarray(self):
        """Return type is numpy.ndarray."""
        s = _sample()
        Q = hkl_to_Q(s, 1, 0, 0)
        assert isinstance(Q, np.ndarray)
        assert Q.shape == (3,)


# ---------------------------------------------------------------------------
# Q_to_hkl
# ---------------------------------------------------------------------------


class TestQToHkl:
    def test_round_trip(self):
        """hkl → Q → hkl reproduces the original indices."""
        s = _sample()
        for hkl in [(1, 0, 0), (0, 1, 0), (0, 0, 1), (1, 1, 0), (2, 1, 3)]:
            Q = hkl_to_Q(s, *hkl)
            result = Q_to_hkl(s, *Q)
            assert result == pytest.approx(hkl, abs=1e-10)

    def test_no_ub_raises(self):
        """Raises ValueError when UB is not set."""
        g = fourcv()
        g.sample.lattice = Lattice(a=A_CUBIC)
        with pytest.raises(ValueError, match=re.escape("sample.UB")):
            Q_to_hkl(g.sample, 1.0, 0.0, 0.0)

    def test_singular_ub_raises(self):
        """Raises ValueError when UB is singular."""
        s = _sample()
        s.UB = np.zeros((3, 3))
        with pytest.raises(ValueError, match=re.escape("singular")):
            Q_to_hkl(s, 1.0, 0.0, 0.0)


# ---------------------------------------------------------------------------
# Q_to_d
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "Qx, Qy, Qz, expected_d, context",
    [
        pytest.param(
            TWO_PI,
            0.0,
            0.0,
            1.0,
            does_not_raise(),
            id="Q=2pi-along-x-gives-d=1",
        ),
        pytest.param(
            0.0,
            TWO_PI,
            0.0,
            1.0,
            does_not_raise(),
            id="Q=2pi-along-y-gives-d=1",
        ),
        pytest.param(
            0.0,
            0.0,
            0.0,
            None,
            pytest.raises(ValueError, match=re.escape("|Q| is zero")),
            id="zero-Q-raises",
        ),
    ],
)
def test_Q_to_d(Qx, Qy, Qz, expected_d, context):
    with context:
        result = Q_to_d(Qx, Qy, Qz)
        assert result == pytest.approx(expected_d, rel=1e-10)


# ---------------------------------------------------------------------------
# d_to_Q_mag
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "d, expected_Q_mag, context",
    [
        pytest.param(1.0, TWO_PI, does_not_raise(), id="d=1-gives-2pi"),
        pytest.param(A_CUBIC, TWO_PI / A_CUBIC, does_not_raise(), id="cubic-lattice"),
        pytest.param(
            0.0,
            None,
            pytest.raises(ValueError, match=re.escape("d must be > 0")),
            id="zero-d",
        ),
        pytest.param(
            -1.0,
            None,
            pytest.raises(ValueError, match=re.escape("d must be > 0")),
            id="negative-d",
        ),
    ],
)
def test_d_to_Q_mag(d, expected_Q_mag, context):
    with context:
        result = d_to_Q_mag(d)
        assert result == pytest.approx(expected_Q_mag, rel=1e-10)


def test_d_Q_mag_round_trip():
    """d → |Q| → d recovers the original d."""
    for d in [0.5, 1.0, 2.0, 3.905, 10.0]:
        Q_mag = d_to_Q_mag(d)
        # Q_to_d expects a Q-vector; supply as (Q_mag, 0, 0) so |Q| = Q_mag
        assert Q_to_d(Q_mag, 0.0, 0.0) == pytest.approx(d, rel=1e-10)

    def test_zero_hkl_raises(self):
        """Raises ValueError for (0,0,0)."""
        s = _sample()
        with pytest.raises(ValueError, match=re.escape("(0, 0, 0)")):
            hkl_to_d(s, 0, 0, 0)

    def test_no_ub_raises(self):
        """Raises ValueError when UB is not set."""
        g = fourcv()
        g.sample.lattice = Lattice(a=A_CUBIC)
        with pytest.raises(ValueError, match=re.escape("sample.UB")):
            hkl_to_d(g.sample, 1, 0, 0)


# ---------------------------------------------------------------------------
# d_to_two_theta / two_theta_to_d
# ---------------------------------------------------------------------------


class TestBraggLaw:
    def test_d_to_two_theta_known(self):
        """2θ for Si (111) at Cu Kα ≈ 28.44°."""
        d_si_111 = 3.1356  # Å
        two_theta = d_to_two_theta(d_si_111, WAVELENGTH)
        assert two_theta == pytest.approx(28.44, abs=0.01)

    def test_round_trip(self):
        """d → 2θ → d recovers the original d."""
        for d in [1.0, 2.0, 3.905]:
            two_theta = d_to_two_theta(d, WAVELENGTH)
            assert two_theta_to_d(two_theta, WAVELENGTH) == pytest.approx(d, rel=1e-10)

    def test_two_theta_to_d_known(self):
        """d from round-trip of d_to_two_theta recovers the original d."""
        d_original = 2.0
        two_theta = d_to_two_theta(d_original, WAVELENGTH)
        assert two_theta_to_d(two_theta, WAVELENGTH) == pytest.approx(
            d_original, rel=1e-10
        )

    @pytest.mark.parametrize(
        "d, wl, context",
        [
            pytest.param(
                0.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("d must be > 0")),
                id="zero-d",
            ),
            pytest.param(
                -1.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("d must be > 0")),
                id="negative-d",
            ),
            pytest.param(
                1.0,
                0.0,
                pytest.raises(ValueError, match=re.escape("wavelength must be > 0")),
                id="zero-wavelength",
            ),
            pytest.param(
                0.1,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("unreachable")),
                id="reflection-unreachable",
            ),
        ],
    )
    def test_d_to_two_theta_raises(self, d, wl, context):
        with context:
            d_to_two_theta(d, wl)

    @pytest.mark.parametrize(
        "two_theta, wl, context",
        [
            pytest.param(
                0.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("(0°, 180°)")),
                id="zero-angle",
            ),
            pytest.param(
                180.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("(0°, 180°)")),
                id="180-angle",
            ),
            pytest.param(
                30.0,
                0.0,
                pytest.raises(ValueError, match=re.escape("wavelength must be > 0")),
                id="zero-wavelength",
            ),
        ],
    )
    def test_two_theta_to_d_raises(self, two_theta, wl, context):
        with context:
            two_theta_to_d(two_theta, wl)


# ---------------------------------------------------------------------------
# hkl_to_two_theta
# ---------------------------------------------------------------------------


class TestHklToTwoTheta:
    def test_cubic_100_cu_kalpha(self):
        """2θ for SrTiO₃ (100) at Cu Kα matches Bragg's law."""
        s = _sample(A_CUBIC)
        two_theta = hkl_to_two_theta(s, 1, 0, 0, WAVELENGTH)
        expected = d_to_two_theta(A_CUBIC, WAVELENGTH)
        assert two_theta == pytest.approx(expected, rel=1e-10)

    def test_zero_hkl_raises(self):
        """Propagates ValueError from hkl_to_d for (0,0,0)."""
        s = _sample()
        with pytest.raises(ValueError, match=re.escape("(0, 0, 0)")):
            hkl_to_two_theta(s, 0, 0, 0, WAVELENGTH)

    def test_no_ub_raises(self):
        """Propagates ValueError when UB is not set."""
        g = fourcv()
        g.sample.lattice = Lattice(a=A_CUBIC)
        with pytest.raises(ValueError, match=re.escape("sample.UB")):
            hkl_to_two_theta(g.sample, 1, 0, 0, WAVELENGTH)

    def test_unreachable_raises(self):
        """Propagates ValueError when reflection is unreachable."""
        s = _sample(a=0.05)  # tiny lattice → large |Q| → unreachable
        with pytest.raises(ValueError, match=re.escape("unreachable")):
            hkl_to_two_theta(s, 1, 0, 0, WAVELENGTH)


# ---------------------------------------------------------------------------
# two_theta_to_Q_mag / Q_mag_to_two_theta
# ---------------------------------------------------------------------------


class TestTwoThetaQMag:
    def test_two_theta_to_Q_mag_known(self):
        """|Q| = 4π sin(θ)/λ for 2θ = 30°, Cu Kα."""
        expected = 4.0 * math.pi * math.sin(math.radians(15.0)) / WAVELENGTH
        assert two_theta_to_Q_mag(30.0, WAVELENGTH) == pytest.approx(
            expected, rel=1e-10
        )

    def test_Q_mag_to_two_theta_round_trip(self):
        """2θ → |Q| → 2θ recovers the original angle."""
        for two_theta in [10.0, 30.0, 60.0, 90.0, 120.0]:
            Q_mag = two_theta_to_Q_mag(two_theta, WAVELENGTH)
            assert Q_mag_to_two_theta(Q_mag, WAVELENGTH) == pytest.approx(
                two_theta, rel=1e-10
            )

    @pytest.mark.parametrize(
        "two_theta, wl, context",
        [
            pytest.param(
                0.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("(0°, 180°)")),
                id="zero-angle",
            ),
            pytest.param(
                180.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("(0°, 180°)")),
                id="180-angle",
            ),
            pytest.param(
                30.0,
                -1.0,
                pytest.raises(ValueError, match=re.escape("wavelength must be > 0")),
                id="negative-wavelength",
            ),
        ],
    )
    def test_two_theta_to_Q_mag_raises(self, two_theta, wl, context):
        with context:
            two_theta_to_Q_mag(two_theta, wl)

    @pytest.mark.parametrize(
        "Q_mag, wl, context",
        [
            pytest.param(
                -1.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("|Q| must be ≥ 0")),
                id="negative-Q-mag",
            ),
            pytest.param(
                100.0,
                WAVELENGTH,
                pytest.raises(ValueError, match=re.escape("unreachable")),
                id="unreachable",
            ),
            pytest.param(
                1.0,
                -1.0,
                pytest.raises(ValueError, match=re.escape("wavelength must be > 0")),
                id="negative-wavelength",
            ),
        ],
    )
    def test_Q_mag_to_two_theta_raises(self, Q_mag, wl, context):
        with context:
            Q_mag_to_two_theta(Q_mag, wl)

    def test_consistency_with_d(self):
        """two_theta_to_Q_mag and d_to_Q_mag agree via Bragg's law."""
        d = 2.0
        two_theta = d_to_two_theta(d, WAVELENGTH)
        Q_from_tth = two_theta_to_Q_mag(two_theta, WAVELENGTH)
        Q_from_d = d_to_Q_mag(d)
        assert Q_from_tth == pytest.approx(Q_from_d, rel=1e-10)


# ---------------------------------------------------------------------------
# Public API exports
# ---------------------------------------------------------------------------


def test_public_api_exports():
    """All engine functions are importable from the conversions submodule."""
    assert ahd.conversions.hkl_to_Q is hkl_to_Q
    assert ahd.conversions.Q_to_hkl is Q_to_hkl
    assert ahd.conversions.Q_to_d is Q_to_d
    assert ahd.conversions.d_to_Q_mag is d_to_Q_mag
    assert ahd.conversions.hkl_to_d is hkl_to_d
    assert ahd.conversions.d_to_two_theta is d_to_two_theta
    assert ahd.conversions.two_theta_to_d is two_theta_to_d
    assert ahd.conversions.hkl_to_two_theta is hkl_to_two_theta
    assert ahd.conversions.two_theta_to_Q_mag is two_theta_to_Q_mag
    assert ahd.conversions.Q_mag_to_two_theta is Q_mag_to_two_theta
