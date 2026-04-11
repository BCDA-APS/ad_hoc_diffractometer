# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for neutron radiation support (issue #8).

Covers:
  - NEUTRON_MEV_ANGSTROM2: value and uncertainty
  - SOURCE_TYPES: valid values
  - neutron_wavelength_to_energy: known value, round-trip, zero/negative raises
  - neutron_energy_to_wavelength: round-trip, zero/negative raises
  - AdHocDiffractometer.source_type: default, set valid, set invalid raises
  - AdHocDiffractometer.energy_units: "keV" for xray, "meV" for neutron
  - AdHocDiffractometer.energy: xray vs neutron values, set neutron, None
  - summary() reports meV for neutron source_type
  - X-ray and neutron energy conversions do not mix
  - Public API exports
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer.radiation import NEUTRON_MEV_ANGSTROM2
from ad_hoc_diffractometer.radiation import NEUTRON_MEV_ANGSTROM2_UNCERTAINTY
from ad_hoc_diffractometer.radiation import SOURCE_TYPES
from ad_hoc_diffractometer.radiation import neutron_energy_to_wavelength
from ad_hoc_diffractometer.radiation import neutron_wavelength_to_energy

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


def test_NEUTRON_MEV_ANGSTROM2_value():
    """h²/(2 m_n) = 81.804 210 235 2 meV·Å² (NIST CODATA 2022)."""
    assert NEUTRON_MEV_ANGSTROM2 == pytest.approx(81.8042102352, rel=1e-10)


def test_NEUTRON_MEV_ANGSTROM2_uncertainty():
    """Uncertainty propagated from CODATA 2022 m_n uncertainty."""
    assert NEUTRON_MEV_ANGSTROM2_UNCERTAINTY == pytest.approx(4.15e-8, rel=1e-2)
    assert NEUTRON_MEV_ANGSTROM2_UNCERTAINTY > 0.0


def test_SOURCE_TYPES():
    """SOURCE_TYPES contains exactly 'xray' and 'neutron'."""
    assert set(SOURCE_TYPES) == {"xray", "neutron"}


# ---------------------------------------------------------------------------
# neutron_wavelength_to_energy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "wavelength, expected_mev, context",
    [
        pytest.param(
            1.0,
            NEUTRON_MEV_ANGSTROM2,
            does_not_raise(),
            id="lambda=1-gives-const",
        ),
        pytest.param(
            1.8,
            NEUTRON_MEV_ANGSTROM2 / 1.8**2,
            does_not_raise(),
            id="thermal-neutron-1.8A",
        ),
        pytest.param(
            0.0,
            None,
            pytest.raises(ValueError, match=re.escape("wavelength must be > 0")),
            id="zero-wavelength",
        ),
        pytest.param(
            -1.0,
            None,
            pytest.raises(ValueError, match=re.escape("wavelength must be > 0")),
            id="negative-wavelength",
        ),
    ],
)
def test_neutron_wavelength_to_energy(wavelength, expected_mev, context):
    with context:
        result = neutron_wavelength_to_energy(wavelength)
        assert result == pytest.approx(expected_mev, rel=1e-10)


def test_neutron_round_trip():
    """λ → E (meV) → λ recovers the original wavelength."""
    for wl in [0.5, 1.0, 1.8, 2.5]:
        assert neutron_energy_to_wavelength(
            neutron_wavelength_to_energy(wl)
        ) == pytest.approx(wl, rel=1e-10)


def test_neutron_energy_to_wavelength_known():
    """E = 25 meV → λ ≈ 1.809 Å."""
    lam = neutron_energy_to_wavelength(25.0)
    assert lam == pytest.approx(math.sqrt(NEUTRON_MEV_ANGSTROM2 / 25.0), rel=1e-10)


@pytest.mark.parametrize(
    "energy_mev, context",
    [
        pytest.param(
            0.0,
            pytest.raises(ValueError, match=re.escape("energy must be > 0 meV")),
            id="zero-energy",
        ),
        pytest.param(
            -1.0,
            pytest.raises(ValueError, match=re.escape("energy must be > 0 meV")),
            id="negative-energy",
        ),
    ],
)
def test_neutron_energy_to_wavelength_raises(energy_mev, context):
    with context:
        neutron_energy_to_wavelength(energy_mev)


# ---------------------------------------------------------------------------
# AdHocDiffractometer.source_type
# ---------------------------------------------------------------------------


class TestSourceType:
    def test_default_is_xray(self):
        """Default source_type is 'xray'."""
        g = fourcv()
        assert g.source_type == "xray"

    def test_set_neutron(self):
        """source_type can be set to 'neutron'."""
        g = fourcv()
        g.source_type = "neutron"
        assert g.source_type == "neutron"

    def test_set_back_to_xray(self):
        """source_type can be switched back to 'xray'."""
        g = fourcv()
        g.source_type = "neutron"
        g.source_type = "xray"
        assert g.source_type == "xray"

    @pytest.mark.parametrize(
        "bad_type, context",
        [
            pytest.param(
                "gamma",
                pytest.raises(
                    ValueError, match=re.escape("source_type must be one of")
                ),
                id="invalid-gamma",
            ),
            pytest.param(
                "spallation",
                pytest.raises(
                    ValueError, match=re.escape("source_type must be one of")
                ),
                id="invalid-spallation",
            ),
            pytest.param(
                "",
                pytest.raises(
                    ValueError, match=re.escape("source_type must be one of")
                ),
                id="empty-string",
            ),
        ],
    )
    def test_invalid_source_type_raises(self, bad_type, context):
        g = fourcv()
        with context:
            g.source_type = bad_type


# ---------------------------------------------------------------------------
# AdHocDiffractometer.energy_units
# ---------------------------------------------------------------------------


class TestEnergyUnits:
    def test_xray_gives_keV(self):
        """energy_units is 'keV' for xray source."""
        g = fourcv()
        assert g.energy_units == "keV"

    def test_neutron_gives_meV(self):
        """energy_units is 'meV' for neutron source."""
        g = fourcv()
        g.source_type = "neutron"
        assert g.energy_units == "meV"


# ---------------------------------------------------------------------------
# AdHocDiffractometer.energy — neutron branch
# ---------------------------------------------------------------------------


class TestEnergyNeutron:
    def test_none_when_wavelength_not_set(self):
        """energy is None for neutron source when wavelength not set."""
        g = fourcv()
        g.source_type = "neutron"
        assert g.energy is None

    def test_neutron_energy_value(self):
        """energy = h²/(2 m_n λ²) in meV for neutron source."""
        g = fourcv()
        g.source_type = "neutron"
        g.wavelength = 1.8
        assert g.energy == pytest.approx(neutron_wavelength_to_energy(1.8), rel=1e-10)

    def test_neutron_energy_units_are_meV(self):
        """energy is in meV for neutron — much smaller than keV X-ray energy."""
        g_xray = fourcv()
        g_xray.wavelength = 1.5406

        g_n = fourcv()
        g_n.source_type = "neutron"
        g_n.wavelength = 1.5406

        # Same wavelength: X-ray gives ~8 keV, neutron gives ~35 meV
        # They are six orders of magnitude apart in different units.
        assert g_xray.energy == pytest.approx(8.0, abs=0.5)  # keV
        assert g_n.energy == pytest.approx(35.0, abs=5.0)  # meV

    def test_set_neutron_energy_updates_wavelength(self):
        """Setting energy (meV) for neutron updates wavelength via de Broglie."""
        g = fourcv()
        g.source_type = "neutron"
        g.energy = 25.0  # meV
        assert g.wavelength == pytest.approx(
            neutron_energy_to_wavelength(25.0), rel=1e-10
        )

    def test_set_neutron_energy_none_clears_wavelength(self):
        """Setting energy = None clears wavelength."""
        g = fourcv()
        g.source_type = "neutron"
        g.wavelength = 1.8
        g.energy = None
        assert g.wavelength is None
        assert g.energy is None

    def test_neutron_round_trip_via_energy(self):
        """wavelength → energy → wavelength round-trip for neutron."""
        g = fourcv()
        g.source_type = "neutron"
        g.wavelength = 1.8
        e = g.energy
        g2 = fourcv()
        g2.source_type = "neutron"
        g2.energy = e
        assert g2.wavelength == pytest.approx(1.8, rel=1e-10)

    def test_switching_source_type_changes_energy(self):
        """Switching source_type changes energy (and energy_units) immediately."""
        g = fourcv()
        g.wavelength = 1.5406
        e_xray = g.energy  # keV

        g.source_type = "neutron"
        e_neutron = g.energy  # meV

        # Same wavelength, different energy formulas — results differ by ~10^6
        assert g.energy_units == "meV"
        assert e_xray != pytest.approx(e_neutron, rel=0.01)


# ---------------------------------------------------------------------------
# summary() reports meV for neutron
# ---------------------------------------------------------------------------


def test_summary_neutron_reports_meV(capsys):
    """summary() shows meV when source_type is neutron."""
    g = fourcv()
    g.source_type = "neutron"
    g.wavelength = 1.8
    g.summary()
    out = capsys.readouterr().out
    assert "meV" in out
    assert "keV" not in out


def test_summary_xray_reports_keV(capsys):
    """summary() shows keV when source_type is xray (default)."""
    g = fourcv()
    g.wavelength = 1.5406
    g.summary()
    out = capsys.readouterr().out
    assert "keV" in out


# ---------------------------------------------------------------------------
# X-ray and neutron conversions must not be mixed
# ---------------------------------------------------------------------------


def test_xray_and_neutron_energies_differ_at_same_wavelength():
    """X-ray and neutron energy formulas give different results (different units)."""
    wl = 1.5406
    e_xray = ahd.wavelength_to_energy(wl)  # keV
    e_neutron = ahd.neutron_wavelength_to_energy(wl)  # meV
    # X-ray: ~8 keV; neutron: ~34 meV — not simply related by a factor of 1000
    assert abs(e_xray * 1000 - e_neutron) > 1.0  # they are NOT equal


# ---------------------------------------------------------------------------
# Public API exports
# ---------------------------------------------------------------------------


def test_public_api_exports():
    """All neutron symbols are importable from the top-level package."""
    assert ahd.NEUTRON_MEV_ANGSTROM2 is NEUTRON_MEV_ANGSTROM2
    assert ahd.NEUTRON_MEV_ANGSTROM2_UNCERTAINTY is NEUTRON_MEV_ANGSTROM2_UNCERTAINTY
    assert ahd.SOURCE_TYPES is SOURCE_TYPES
    assert ahd.neutron_wavelength_to_energy is neutron_wavelength_to_energy
    assert ahd.neutron_energy_to_wavelength is neutron_energy_to_wavelength
