"""
Unit tests for ad_hoc_diffractometer.spec — SPEC #G1 format for fourc geometry.

Covers:
  - parse_fourc_g1(): field count, field values, tag stripping, error cases
  - emit_fourc_g1(): starts with '#G1', correct field count, round-trip
  - g1_to_sample(): lattice, reflections, wavelength, or1/or2 designation
  - sample_to_g1(): reciprocal lattice parameters, angle mapping, error cases
  - Round-trip: parse → g1_to_sample → sample_to_g1 → emit reproduces key fields
  - Compatibility with real data from Align4Pete.spec (three historical #G1 lines)
  - UB round-trip: parse #G1 → g1_to_sample → ub_from_two_reflections_bl1967
    → verify UB @ or1_hkl is parallel to angles_to_phi_vector(or1_angles)
"""

import re

import numpy as np
import pytest

from ad_hoc_diffractometer import angles_to_phi_vector
from ad_hoc_diffractometer import emit_fourc_g1
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import g1_to_sample
from ad_hoc_diffractometer import parse_fourc_g1
from ad_hoc_diffractometer import sample_to_g1
from ad_hoc_diffractometer import ub_from_two_reflections_bl1967
from ad_hoc_diffractometer.spec import FourcG1

# ---------------------------------------------------------------------------
# Reference #G1 lines from Align4Pete.spec (three historical variants)
# ---------------------------------------------------------------------------

# Scans 1–7: initial fake orientation, or2 = (1,0,0) at 60/30/0/0
_G1_LINE_A = (
    "#G1 4.785 4.785 12.991 90 90 120 "
    "1.516237713 1.516237713 0.483656786 90 90 60 "
    "0 0 6  1 0 0  "
    "41.94188 20.97 90 0  0 0  "
    "60 30 0 0  0 0  "
    "1.549802558 1.549802558  0 0"
)

# Scan 8+: measured (006) position entered, or2 still (1,0,0)
_G1_LINE_B = (
    "#G1 4.785 4.785 12.991 90 90 120 "
    "1.516237713 1.516237713 0.483656786 90 90 60 "
    "0 0 6  1 0 0  "
    "41.939375 20.3653625 89.32 0  0 0  "
    "60 30 0 0  0 0  "
    "1.549802558 1.549802558  0 0"
)

# Later: or2 changed to (1,0,4)
_G1_LINE_C = (
    "#G1 4.785 4.785 12.991 90 90 120 "
    "1.516237713 1.516237713 0.483656786 90 90 60 "
    "0 0 6  1 0 4  "
    "41.939375 20.3653625 89.32 0  0 0  "
    "35.392375 17.6428 50.8925 29.95  0 0  "
    "1.549802558 1.549802558  0 0"
)

# ---------------------------------------------------------------------------
# parse_fourc_g1()
# ---------------------------------------------------------------------------


class TestParseFourcG1:
    def test_returns_fourcg1_instance(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        assert isinstance(g1, FourcG1)

    def test_direct_lattice_a(self):
        assert parse_fourc_g1(_G1_LINE_A).a == pytest.approx(4.785)

    def test_direct_lattice_b(self):
        assert parse_fourc_g1(_G1_LINE_A).b == pytest.approx(4.785)

    def test_direct_lattice_c(self):
        assert parse_fourc_g1(_G1_LINE_A).c == pytest.approx(12.991)

    def test_direct_lattice_alpha(self):
        assert parse_fourc_g1(_G1_LINE_A).alpha == pytest.approx(90.0)

    def test_direct_lattice_gamma(self):
        assert parse_fourc_g1(_G1_LINE_A).gamma == pytest.approx(120.0)

    def test_reciprocal_a_star(self):
        assert parse_fourc_g1(_G1_LINE_A).a_star == pytest.approx(1.516237713)

    def test_reciprocal_c_star(self):
        assert parse_fourc_g1(_G1_LINE_A).c_star == pytest.approx(0.483656786)

    def test_reciprocal_gamma_star(self):
        assert parse_fourc_g1(_G1_LINE_A).gamma_star == pytest.approx(60.0)

    def test_or1_hkl_line_a(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        assert (g1.or1_h, g1.or1_k, g1.or1_l) == (0.0, 0.0, 6.0)

    def test_or2_hkl_line_a(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        assert (g1.or2_h, g1.or2_k, g1.or2_l) == (1.0, 0.0, 0.0)

    def test_or2_hkl_line_c(self):
        g1 = parse_fourc_g1(_G1_LINE_C)
        assert (g1.or2_h, g1.or2_k, g1.or2_l) == (1.0, 0.0, 4.0)

    def test_or1_angles_line_a(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        assert g1.or1_two_theta == pytest.approx(41.94188)
        assert g1.or1_omega == pytest.approx(20.97)
        assert g1.or1_chi == pytest.approx(90.0)
        assert g1.or1_phi == pytest.approx(0.0)

    def test_or2_angles_line_c(self):
        g1 = parse_fourc_g1(_G1_LINE_C)
        assert g1.or2_two_theta == pytest.approx(35.392375)
        assert g1.or2_omega == pytest.approx(17.6428)
        assert g1.or2_chi == pytest.approx(50.8925)
        assert g1.or2_phi == pytest.approx(29.95)

    def test_lambda1(self):
        assert parse_fourc_g1(_G1_LINE_A).lambda1 == pytest.approx(1.549802558)

    def test_lambda2(self):
        assert parse_fourc_g1(_G1_LINE_A).lambda2 == pytest.approx(1.549802558)

    def test_unused_zeros(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        assert g1.unused1 == (0.0, 0.0)
        assert g1.unused2 == (0.0, 0.0)
        assert g1.unused3 == (0.0, 0.0)

    def test_accepts_line_without_tag(self):
        """parse_fourc_g1 works whether or not the '#G1' prefix is present."""
        raw = _G1_LINE_A.replace("#G1 ", "", 1)
        g1 = parse_fourc_g1(raw)
        assert g1.a == pytest.approx(4.785)

    def test_wrong_field_count_raises(self):
        bad = "#G1 1 2 3"  # only 3 fields
        with pytest.raises(ValueError, match=re.escape("34 numeric fields")):
            parse_fourc_g1(bad)

    def test_line_b_or1_omega(self):
        """Scan 8+ measured omega = 20.3653625."""
        g1 = parse_fourc_g1(_G1_LINE_B)
        assert g1.or1_omega == pytest.approx(20.3653625)

    def test_line_b_or1_chi(self):
        """Scan 8+ measured chi = 89.32."""
        g1 = parse_fourc_g1(_G1_LINE_B)
        assert g1.or1_chi == pytest.approx(89.32)


# ---------------------------------------------------------------------------
# emit_fourc_g1()
# ---------------------------------------------------------------------------


class TestEmitFourcG1:
    def test_starts_with_g1_tag(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        assert emit_fourc_g1(g1).startswith("#G1")

    def test_field_count(self):
        g1 = parse_fourc_g1(_G1_LINE_A)
        parts = emit_fourc_g1(g1).split()
        assert len(parts) == 35  # '#G1' + 34 fields

    def test_round_trip_a(self):
        """parse → emit → parse gives identical FourcG1."""
        g1 = parse_fourc_g1(_G1_LINE_A)
        g1_rt = parse_fourc_g1(emit_fourc_g1(g1))
        assert g1_rt.a == pytest.approx(g1.a)
        assert g1_rt.c == pytest.approx(g1.c)
        assert g1_rt.a_star == pytest.approx(g1.a_star, rel=1e-6)
        assert g1_rt.or1_h == g1.or1_h
        assert g1_rt.or1_chi == pytest.approx(g1.or1_chi)
        assert g1_rt.lambda1 == pytest.approx(g1.lambda1)

    def test_round_trip_c(self):
        """Round-trip preserves or2 hkl and angles for line C."""
        g1 = parse_fourc_g1(_G1_LINE_C)
        g1_rt = parse_fourc_g1(emit_fourc_g1(g1))
        assert g1_rt.or2_l == pytest.approx(g1.or2_l)
        assert g1_rt.or2_chi == pytest.approx(g1.or2_chi)
        assert g1_rt.or2_phi == pytest.approx(g1.or2_phi)


# ---------------------------------------------------------------------------
# g1_to_sample()
# ---------------------------------------------------------------------------


class TestG1ToSample:
    @pytest.fixture
    def geom_a(self):
        g = fourcv()
        g1_to_sample(parse_fourc_g1(_G1_LINE_A), g)
        return g

    def test_wavelength_set(self, geom_a):
        assert geom_a.wavelength == pytest.approx(1.549802558)

    def test_lattice_a(self, geom_a):
        assert geom_a.sample.lattice.a == pytest.approx(4.785)

    def test_lattice_c(self, geom_a):
        assert geom_a.sample.lattice.c == pytest.approx(12.991)

    def test_lattice_gamma(self, geom_a):
        assert geom_a.sample.lattice.gamma == pytest.approx(120.0)

    def test_or1_designated(self, geom_a):
        ors = geom_a.sample.reflections.orienting_reflections
        assert len(ors) >= 1
        assert ors[0].name == "or1"

    def test_or2_designated(self, geom_a):
        ors = geom_a.sample.reflections.orienting_reflections
        assert len(ors) == 2
        assert ors[1].name == "or2"

    def test_or1_hkl(self, geom_a):
        or1 = geom_a.sample.reflections["or1"]
        assert or1.hkl == pytest.approx((0.0, 0.0, 6.0))

    def test_or2_hkl(self, geom_a):
        or2 = geom_a.sample.reflections["or2"]
        assert or2.hkl == pytest.approx((1.0, 0.0, 0.0))

    def test_or1_angles(self, geom_a):
        ang = geom_a.sample.reflections["or1"].angles
        assert ang["two_theta"] == pytest.approx(41.94188)
        assert ang["omega"] == pytest.approx(20.97)
        assert ang["chi"] == pytest.approx(90.0)
        assert ang["phi"] == pytest.approx(0.0)

    def test_or1_wavelength_stored(self, geom_a):
        assert geom_a.sample.reflections["or1"].wavelength == pytest.approx(1.549802558)

    def test_line_c_or2_hkl(self):
        """Line C: or2 = (1,0,4)."""
        g = fourcv()
        g1_to_sample(parse_fourc_g1(_G1_LINE_C), g)
        or2 = g.sample.reflections["or2"]
        assert or2.hkl == pytest.approx((1.0, 0.0, 4.0))

    def test_line_c_or2_angles(self):
        """Line C: or2 angles match the file values."""
        g = fourcv()
        g1_to_sample(parse_fourc_g1(_G1_LINE_C), g)
        ang = g.sample.reflections["or2"].angles
        assert ang["two_theta"] == pytest.approx(35.392375)
        assert ang["chi"] == pytest.approx(50.8925)
        assert ang["phi"] == pytest.approx(29.95)

    def test_replaces_existing_or1(self):
        """Calling g1_to_sample twice overwrites or1/or2 without error."""
        g = fourcv()
        g1_to_sample(parse_fourc_g1(_G1_LINE_A), g)
        g1_to_sample(parse_fourc_g1(_G1_LINE_B), g)
        assert g.sample.reflections["or1"].angles["omega"] == pytest.approx(20.3653625)


# ---------------------------------------------------------------------------
# sample_to_g1()
# ---------------------------------------------------------------------------


class TestSampleToG1:
    @pytest.fixture
    def g1_rt(self):
        """Parse line A, apply to geometry, emit back."""
        g = fourcv()
        g1_to_sample(parse_fourc_g1(_G1_LINE_A), g)
        return sample_to_g1(g), parse_fourc_g1(_G1_LINE_A)

    def test_direct_lattice_a(self, g1_rt):
        rt, orig = g1_rt
        assert rt.a == pytest.approx(orig.a)

    def test_direct_lattice_c(self, g1_rt):
        rt, orig = g1_rt
        assert rt.c == pytest.approx(orig.c)

    def test_direct_lattice_gamma(self, g1_rt):
        rt, orig = g1_rt
        assert rt.gamma == pytest.approx(orig.gamma)

    def test_recip_a_star(self, g1_rt):
        rt, orig = g1_rt
        assert rt.a_star == pytest.approx(orig.a_star, rel=1e-5)

    def test_recip_c_star(self, g1_rt):
        rt, orig = g1_rt
        assert rt.c_star == pytest.approx(orig.c_star, rel=1e-5)

    def test_recip_gamma_star(self, g1_rt):
        rt, orig = g1_rt
        assert rt.gamma_star == pytest.approx(orig.gamma_star, abs=0.01)

    def test_or1_hkl(self, g1_rt):
        rt, orig = g1_rt
        assert (rt.or1_h, rt.or1_k, rt.or1_l) == pytest.approx(
            (orig.or1_h, orig.or1_k, orig.or1_l)
        )

    def test_or1_angles(self, g1_rt):
        rt, orig = g1_rt
        assert rt.or1_two_theta == pytest.approx(orig.or1_two_theta)
        assert rt.or1_omega == pytest.approx(orig.or1_omega)
        assert rt.or1_chi == pytest.approx(orig.or1_chi)
        assert rt.or1_phi == pytest.approx(orig.or1_phi)

    def test_lambda1(self, g1_rt):
        rt, orig = g1_rt
        assert rt.lambda1 == pytest.approx(orig.lambda1)

    def test_no_wavelength_raises(self):
        g = fourcv()
        g1_to_sample(parse_fourc_g1(_G1_LINE_A), g)
        g.wavelength = None
        with pytest.raises(ValueError, match=re.escape("wavelength")):
            sample_to_g1(g)

    def test_no_or1_raises(self):
        g = fourcv()
        g.wavelength = 1.5406
        with pytest.raises(ValueError, match=re.escape("orienting")):
            sample_to_g1(g)


# ---------------------------------------------------------------------------
# UB round-trip: parse → g1_to_sample → ub_from_two_reflections_bl1967
# ---------------------------------------------------------------------------


class TestUBRoundTrip:
    """
    Verify that after loading a #G1 line and computing UB from the two
    orienting reflections, UB @ or1_hkl is parallel to the Q_phi vector
    obtained from the or1 motor angles.  This is the fundamental BL1967
    guarantee for the primary reflection.
    """

    @pytest.mark.parametrize(
        "line,desc",
        [
            pytest.param(_G1_LINE_A, "line-A-fake-orientation", id="line-A"),
            pytest.param(_G1_LINE_B, "line-B-measured-006", id="line-B"),
            pytest.param(_G1_LINE_C, "line-C-second-refl-104", id="line-C"),
        ],
    )
    def test_ub_or1_direction_parallel_to_q_phi(self, line, desc):
        """UB @ or1_hkl must be parallel to Q_phi of or1 (BL1967 guarantee)."""
        g = fourcv()
        g1_to_sample(parse_fourc_g1(line), g)
        UB = ub_from_two_reflections_bl1967(g.sample)
        or1 = g.sample.reflections["or1"]
        h_phi = angles_to_phi_vector(g, **or1.angles)
        q_calc = UB @ np.array(or1.hkl, dtype=float)
        # directions must agree (ignore magnitude)
        h_hat = h_phi / np.linalg.norm(h_phi)
        q_hat = q_calc / np.linalg.norm(q_calc)
        np.testing.assert_allclose(
            np.abs(np.dot(h_hat, q_hat)),
            1.0,
            atol=1e-6,
            err_msg=f"or1 direction mismatch for {desc}",
        )

    @pytest.mark.parametrize(
        "line,desc",
        [
            pytest.param(_G1_LINE_A, "line-A-fake-orientation", id="line-A"),
            pytest.param(_G1_LINE_B, "line-B-measured-006", id="line-B"),
            pytest.param(_G1_LINE_C, "line-C-second-refl-104", id="line-C"),
        ],
    )
    def test_ub_is_3x3(self, line, desc):
        g = fourcv()
        g1_to_sample(parse_fourc_g1(line), g)
        UB = ub_from_two_reflections_bl1967(g.sample)
        assert UB.shape == (3, 3), f"UB shape wrong for {desc}"

    @pytest.mark.parametrize(
        "line,desc",
        [
            pytest.param(_G1_LINE_A, "line-A-fake-orientation", id="line-A"),
            pytest.param(_G1_LINE_B, "line-B-measured-006", id="line-B"),
            pytest.param(_G1_LINE_C, "line-C-second-refl-104", id="line-C"),
        ],
    )
    def test_u_is_orthonormal(self, line, desc):
        """U from BL1967 must be an orthonormal rotation matrix."""
        g = fourcv()
        g1_to_sample(parse_fourc_g1(line), g)
        ub_from_two_reflections_bl1967(g.sample)
        U = g.sample.U
        np.testing.assert_allclose(
            U.T @ U,
            np.eye(3),
            atol=1e-8,
            err_msg=f"U not orthonormal for {desc}",
        )
        assert abs(np.linalg.det(U) - 1.0) < 1e-8, f"det(U) != 1 for {desc}"
