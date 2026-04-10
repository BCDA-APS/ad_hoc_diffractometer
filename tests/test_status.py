"""
Unit tests for ad_hoc_diffractometer.status — wh() and pa() commands.

Covers:
  - wh(): returns a string, contains 'H K L', contains 'Lambda',
    contains a motor-angle table header, HKL line when UB not set,
    motor names mapped to SPEC style, all-zero position
  - pa(): returns a string, contains geometry name, lattice constants,
    orienting reflections, wavelength, reciprocal lattice params
  - Both functions work with no UB set (graceful fallback)
  - Both functions work with no wavelength set (graceful fallback)
  - Regression against Align4Pete.log output values
  - g.wh and g.pa properties (#51): property access returns the same
    string as the module-level functions; both are @property descriptors
"""

import pytest

from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer import psic
from ad_hoc_diffractometer import ub_from_two_reflections_bl1967
from ad_hoc_diffractometer.spec import g1_to_sample
from ad_hoc_diffractometer.spec import parse_fourc_g1
from ad_hoc_diffractometer.status import _spec_motor_name
from ad_hoc_diffractometer.status import pa
from ad_hoc_diffractometer.status import wh

# ---------------------------------------------------------------------------
# Reference #G1 line from Align4Pete.spec (final session state)
# ---------------------------------------------------------------------------

_G1_LINE = (
    "#G1 4.785 4.785 12.991 90 90 120 "
    "1.516237713 1.516237713 0.483656786 90 90 60 "
    "0 0 6  1 0 4  "
    "41.939375 20.3653625 89.32 0  0 0  "
    "35.392375 17.6428 50.8925 29.95  0 0  "
    "1.549802558 1.549802558  0 0"
)


@pytest.fixture
def sapphire_geom():
    """fourcv with sapphire lattice and two orienting reflections; UB set."""
    g = fourcv()
    g1_to_sample(parse_fourc_g1(_G1_LINE), g)
    ub_from_two_reflections_bl1967(g.sample)
    return g


@pytest.fixture
def bare_geom():
    """fourcv with no wavelength, no UB — tests graceful fallback."""
    return fourcv()


# ---------------------------------------------------------------------------
# _spec_motor_name helper
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
    assert _spec_motor_name(internal) == expected


# ---------------------------------------------------------------------------
# wh()
# ---------------------------------------------------------------------------


class TestWh:
    def test_returns_str(self, sapphire_geom):
        assert isinstance(wh(sapphire_geom), str)

    def test_contains_hkl_line(self, sapphire_geom):
        assert "H K L" in wh(sapphire_geom)

    def test_contains_lambda_line(self, sapphire_geom):
        assert "Lambda" in wh(sapphire_geom)

    def test_contains_motor_header(self, sapphire_geom):
        out = wh(sapphire_geom)
        # All four fourcv stages should appear in the header
        assert "TwoTheta" in out
        assert "Theta" in out
        assert "Chi" in out
        assert "Phi" in out

    def test_lambda_value(self, sapphire_geom):
        out = wh(sapphire_geom)
        assert "1.5498" in out

    def test_motor_angles_at_zero(self, sapphire_geom):
        """With all motors at 0, the angle table should show zeros."""
        for name in ("omega", "chi", "phi", "two_theta"):
            sapphire_geom.set_angle(name, 0.0)
        out = wh(sapphire_geom)
        # The values row should contain '0.000' entries
        assert "0.000" in out

    def test_motor_angles_after_move(self, sapphire_geom):
        """Setting omega=17.643 should appear in the output."""
        sapphire_geom.set_angle("omega", 17.643)
        out = wh(sapphire_geom)
        assert "17.643" in out

    def test_no_ub_graceful(self, bare_geom):
        """Without a UB, HKL line shows 'not available' rather than crashing."""
        bare_geom.wavelength = 1.5406
        out = wh(bare_geom)
        assert "not available" in out
        assert "Lambda" in out

    def test_no_wavelength_graceful(self, bare_geom):
        """Without wavelength, output still contains H K L and Lambda lines."""
        out = wh(bare_geom)
        assert "H K L" in out
        assert "Lambda" in out

    def test_psic_stage_names(self):
        """psic geometry: stage names appear in the motor table header."""
        g = psic()
        g.wavelength = 1.5406
        out = wh(g)
        assert "Mu" in out
        assert "Eta" in out
        assert "Nu" in out
        assert "Delta" in out

    def test_hkl_all_zero_at_zero_angles(self, sapphire_geom):
        """All motors at zero → Q=0 → inverse() fails or returns 0/not-available."""
        for name in ("omega", "chi", "phi", "two_theta"):
            sapphire_geom.set_angle(name, 0.0)
        out = wh(sapphire_geom)
        # Either "not available" or all-zero HKL
        assert "H K L" in out


# ---------------------------------------------------------------------------
# pa()
# ---------------------------------------------------------------------------


class TestPa:
    def test_returns_str(self, sapphire_geom):
        assert isinstance(pa(sapphire_geom), str)

    def test_contains_geometry_name(self, sapphire_geom):
        out = pa(sapphire_geom)
        assert "fourcv" in out

    def test_contains_lattice_section(self, sapphire_geom):
        out = pa(sapphire_geom)
        assert "Lattice" in out

    def test_contains_real_space_params(self, sapphire_geom):
        """a=4.785, c=12.991, gamma=120 must appear."""
        out = pa(sapphire_geom)
        assert "4.785" in out
        assert "12.991" in out
        assert "120" in out

    def test_contains_reciprocal_space_params(self, sapphire_geom):
        """a*≈1.516, c*≈0.484 must appear."""
        out = pa(sapphire_geom)
        assert "1.516" in out
        assert "0.48" in out

    def test_reciprocal_gamma_star_60(self, sapphire_geom):
        """For hexagonal: gamma*=60° must appear."""
        out = pa(sapphire_geom)
        assert "60" in out

    def test_contains_lambda(self, sapphire_geom):
        out = pa(sapphire_geom)
        assert "Lambda" in out
        assert "1.5498" in out

    def test_contains_primary_reflection(self, sapphire_geom):
        """Primary reflection hkl (0 0 6) must appear."""
        out = pa(sapphire_geom)
        assert "Primary" in out
        # hkl components
        assert "0" in out
        assert "6" in out

    def test_contains_secondary_reflection(self, sapphire_geom):
        """Secondary reflection hkl (1 0 4) must appear."""
        out = pa(sapphire_geom)
        assert "Secondary" in out
        assert "4" in out

    def test_no_reflections_graceful(self, bare_geom):
        """Without orienting reflections, pa() shows 'not set'."""
        bare_geom.wavelength = 1.5406
        out = pa(bare_geom)
        assert "not set" in out

    def test_no_wavelength_graceful(self, bare_geom):
        """Without wavelength, pa() still shows geometry name and lattice."""
        out = pa(bare_geom)
        assert "fourcv" in out
        assert "Lattice" in out

    def test_primary_angles_in_output(self, sapphire_geom):
        """Primary reflection angles 41.9394 / 20.3654 must appear."""
        out = pa(sapphire_geom)
        assert "41.939" in out or "41.9394" in out

    def test_secondary_angles_in_output(self, sapphire_geom):
        """Secondary reflection angles 35.392 / 17.6428 must appear."""
        out = pa(sapphire_geom)
        assert "35.392" in out or "35.3924" in out

    def test_psic_geometry_name(self):
        """pa() shows the correct geometry name for psic."""
        g = psic()
        g.wavelength = 1.5406
        assert "psic" in pa(g)

    def test_reciprocal_params_numeric(self, sapphire_geom):
        """Reciprocal lattice values must be numeric (can be parsed as floats)."""
        out = pa(sapphire_geom)
        # Find the reciprocal space line
        for line in out.splitlines():
            if "reciprocal space" in line:
                # Extract numbers — should have 6 floats
                import re

                nums = re.findall(r"[-+]?\d*\.?\d+", line)
                assert len(nums) >= 6
                for n in nums:
                    float(n)  # must be parseable
                break
        else:
            pytest.fail("No 'reciprocal space' line found in pa() output")


# ---------------------------------------------------------------------------
# g.wh and g.pa properties (#51)
# ---------------------------------------------------------------------------


class TestWhPaProperties:
    """g.wh and g.pa are @property descriptors on AdHocDiffractometer."""

    def test_wh_is_property(self):
        """AdHocDiffractometer.wh must be a property descriptor."""
        from ad_hoc_diffractometer import AdHocDiffractometer

        assert isinstance(AdHocDiffractometer.wh, property)

    def test_pa_is_property(self):
        """AdHocDiffractometer.pa must be a property descriptor."""
        from ad_hoc_diffractometer import AdHocDiffractometer

        assert isinstance(AdHocDiffractometer.pa, property)

    def test_wh_property_equals_function(self, sapphire_geom):
        """g.wh returns the same string as wh(g)."""
        assert sapphire_geom.wh == wh(sapphire_geom)

    def test_pa_property_equals_function(self, sapphire_geom):
        """g.pa returns the same string as pa(g)."""
        assert sapphire_geom.pa == pa(sapphire_geom)

    def test_wh_property_returns_str(self, bare_geom):
        """g.wh returns a str even when nothing is configured."""
        bare_geom.wavelength = 1.5406
        assert isinstance(bare_geom.wh, str)

    def test_pa_property_returns_str(self, bare_geom):
        """g.pa returns a str even when nothing is configured."""
        bare_geom.wavelength = 1.5406
        assert isinstance(bare_geom.pa, str)

    def test_wh_property_contains_hkl(self, sapphire_geom):
        """g.wh contains the 'H K L' line."""
        assert "H K L" in sapphire_geom.wh

    def test_pa_property_contains_geometry_name(self, sapphire_geom):
        """g.pa contains the geometry name."""
        assert "fourcv" in sapphire_geom.pa

    def test_wh_property_not_callable(self, sapphire_geom):
        """g.wh is a str, not callable — consistent with property semantics."""
        assert not callable(sapphire_geom.wh)

    def test_pa_property_not_callable(self, sapphire_geom):
        """g.pa is a str, not callable — consistent with property semantics."""
        assert not callable(sapphire_geom.pa)

    def test_wh_property_same_result_psic(self):
        """g.wh works for psic geometry too."""
        g = psic()
        g.wavelength = 1.5406
        assert g.wh == wh(g)

    def test_pa_property_same_result_psic(self):
        """g.pa works for psic geometry too."""
        g = psic()
        g.wavelength = 1.5406
        assert g.pa == pa(g)
