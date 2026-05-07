# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
Unit tests for ad_hoc_diffractometer.radiation (issues #21 and #8).

Covers:
  - HC_KEV_ANGSTROM, HC_KEV_ANGSTROM_UNCERTAINTY: value and exactness
  - NEUTRON_MEV_ANGSTROM2, NEUTRON_MEV_ANGSTROM2_UNCERTAINTY: value and uncertainty
  - SOURCE_TYPES: valid values
  - XRAY_LINES: all keys present, values positive, Cu Kα known value
  - wavelength_to_energy: known value, round-trip, zero/negative raises
  - energy_to_wavelength: round-trip, zero/negative raises
  - wavelength_to_wavenumber: known value, round-trip, zero/negative raises
  - wavenumber_to_wavelength: round-trip, zero/negative raises
  - neutron_wavelength_to_energy: known value, round-trip, zero/negative raises
  - neutron_energy_to_wavelength: round-trip, zero/negative raises
  - AdHocDiffractometer.source_type: default, set valid, set invalid raises
  - AdHocDiffractometer.energy_units: "keV" for xray, "meV" for neutron
  - AdHocDiffractometer.energy: xray and neutron values, set, None
  - AdHocDiffractometer.wavenumber property: get, set, None
  - summary() includes energy and correct units for each source type
  - X-ray and neutron energy conversions do not mix
  - Public API exports
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import pytest
from helpers import fourcv

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.radiation import HC_KEV_ANGSTROM
from ad_hoc_diffractometer.radiation import HC_KEV_ANGSTROM_UNCERTAINTY
from ad_hoc_diffractometer.radiation import NEUTRON_MEV_ANGSTROM2
from ad_hoc_diffractometer.radiation import NEUTRON_MEV_ANGSTROM2_UNCERTAINTY
from ad_hoc_diffractometer.radiation import SOURCE_TYPES
from ad_hoc_diffractometer.radiation import XRAY_LINES
from ad_hoc_diffractometer.radiation import energy_to_wavelength
from ad_hoc_diffractometer.radiation import neutron_energy_to_wavelength
from ad_hoc_diffractometer.radiation import neutron_wavelength_to_energy
from ad_hoc_diffractometer.radiation import wavelength_to_energy
from ad_hoc_diffractometer.radiation import wavelength_to_wavenumber
from ad_hoc_diffractometer.radiation import wavenumber_to_wavelength

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


def test_HC_KEV_ANGSTROM_value():
    """hc = 12.398 419 843 320 026 keV·Å (exact, NIST CODATA 2022)."""
    assert HC_KEV_ANGSTROM == pytest.approx(12.398419843320026, rel=1e-15)


def test_HC_KEV_ANGSTROM_uncertainty():
    """HC_KEV_ANGSTROM is exact — uncertainty is zero."""
    assert HC_KEV_ANGSTROM_UNCERTAINTY == 0.0


def test_XRAY_LINES_keys():
    """XRAY_LINES contains all expected emission line keys."""
    expected_keys = {
        "Cu_Ka",
        "Cu_Ka1",
        "Cu_Ka2",
        "Mo_Ka",
        "Mo_Ka1",
        "Ag_Ka",
        "Ag_Ka1",
        "Co_Ka",
        "Co_Ka1",
    }
    assert expected_keys.issubset(set(XRAY_LINES.keys()))


def test_XRAY_LINES_positive():
    """All XRAY_LINES wavelengths are positive."""
    for name, wl in XRAY_LINES.items():
        assert wl > 0.0, f"{name} wavelength must be > 0"


def test_Cu_Ka_known_value():
    """Cu Kα weighted mean ≈ 1.5406 Å."""
    assert XRAY_LINES["Cu_Ka"] == pytest.approx(1.5406, abs=0.0001)


def test_Cu_Ka1_less_than_Ka2():
    """Cu Kα1 < Cu Kα2 (Kα1 is the shorter-wavelength component)."""
    assert XRAY_LINES["Cu_Ka1"] < XRAY_LINES["Cu_Ka2"]


def test_Ka_between_Ka1_Ka2():
    """Weighted mean Kα lies between Kα1 and Kα2."""
    for element in ("Cu", "Mo"):
        ka1 = XRAY_LINES[f"{element}_Ka1"]
        ka = XRAY_LINES[f"{element}_Ka"]
        assert ka1 <= ka, f"{element} Kα should be ≥ Kα1"


# ---------------------------------------------------------------------------
# wavelength_to_energy
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "wavelength, expected_kev, context",
    [
        pytest.param(
            1.0,
            HC_KEV_ANGSTROM,
            does_not_raise(),
            id="lambda=1-gives-hc",
        ),
        pytest.param(
            XRAY_LINES["Cu_Ka"],
            None,
            does_not_raise(),
            id="Cu-Ka-positive-energy",
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
def test_wavelength_to_energy(wavelength, expected_kev, context):
    with context:
        result = wavelength_to_energy(wavelength)
        if expected_kev is not None:
            assert result == pytest.approx(expected_kev, rel=1e-10)
        else:
            assert result > 0.0


def test_wavelength_to_energy_round_trip():
    """wavelength → energy → wavelength recovers the original."""
    for wl in [0.5, 1.0, 1.5406, 2.0]:
        assert energy_to_wavelength(wavelength_to_energy(wl)) == pytest.approx(
            wl, rel=1e-10
        )


# ---------------------------------------------------------------------------
# energy_to_wavelength
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "energy, context",
    [
        pytest.param(8.048, does_not_raise(), id="Cu-Ka-energy"),
        pytest.param(
            0.0,
            pytest.raises(ValueError, match=re.escape("energy must be > 0")),
            id="zero-energy",
        ),
        pytest.param(
            -1.0,
            pytest.raises(ValueError, match=re.escape("energy must be > 0")),
            id="negative-energy",
        ),
    ],
)
def test_energy_to_wavelength(energy, context):
    with context:
        result = energy_to_wavelength(energy)
        assert result > 0.0


# ---------------------------------------------------------------------------
# wavelength_to_wavenumber
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "wavelength, expected_k, context",
    [
        pytest.param(
            1.0,
            2.0 * math.pi,
            does_not_raise(),
            id="lambda=1-gives-2pi",
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
def test_wavelength_to_wavenumber(wavelength, expected_k, context):
    with context:
        result = wavelength_to_wavenumber(wavelength)
        if expected_k is not None:
            assert result == pytest.approx(expected_k, rel=1e-10)


def test_wavelength_to_wavenumber_round_trip():
    """wavelength → wavenumber → wavelength recovers the original."""
    for wl in [0.5, 1.0, 1.5406, 2.0]:
        assert wavenumber_to_wavelength(wavelength_to_wavenumber(wl)) == pytest.approx(
            wl, rel=1e-10
        )


# ---------------------------------------------------------------------------
# wavenumber_to_wavelength
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "k, context",
    [
        pytest.param(2.0 * math.pi, does_not_raise(), id="k=2pi"),
        pytest.param(
            0.0,
            pytest.raises(ValueError, match=re.escape("wavenumber must be > 0")),
            id="zero-k",
        ),
        pytest.param(
            -1.0,
            pytest.raises(ValueError, match=re.escape("wavenumber must be > 0")),
            id="negative-k",
        ),
    ],
)
def test_wavenumber_to_wavelength(k, context):
    with context:
        result = wavenumber_to_wavelength(k)
        assert result > 0.0


# ---------------------------------------------------------------------------
# AdHocDiffractometer.energy property
# ---------------------------------------------------------------------------


class TestEnergyKevProperty:
    def test_none_when_wavelength_not_set(self):
        """energy is None when wavelength is not set."""
        g = fourcv()
        assert g.energy is None

    def test_value_matches_conversion(self):
        """energy = wavelength_to_energy(wavelength)."""
        g = fourcv()
        g.wavelength = 1.5406
        assert g.energy == pytest.approx(wavelength_to_energy(1.5406), rel=1e-10)

    def test_set_updates_wavelength(self):
        """Setting energy updates wavelength via λ = hc/E."""
        g = fourcv()
        g.energy = 8.048
        assert g.wavelength == pytest.approx(energy_to_wavelength(8.048), rel=1e-10)

    def test_set_none_clears_wavelength(self):
        """Setting energy = None clears wavelength."""
        g = fourcv()
        g.wavelength = 1.5406
        g.energy = None
        assert g.wavelength is None
        assert g.energy is None

    def test_round_trip_via_wavelength(self):
        """wavelength → energy → wavelength recovers the original."""
        g = fourcv()
        g.wavelength = 1.5406
        e = g.energy
        g2 = fourcv()
        g2.energy = e
        assert g2.wavelength == pytest.approx(1.5406, rel=1e-6)


# ---------------------------------------------------------------------------
# AdHocDiffractometer.wavenumber property
# ---------------------------------------------------------------------------


class TestWavenumberProperty:
    def test_none_when_wavelength_not_set(self):
        """wavenumber is None when wavelength is not set."""
        g = fourcv()
        assert g.wavenumber is None

    def test_value_matches_conversion(self):
        """wavenumber = 2π / wavelength."""
        g = fourcv()
        g.wavelength = 1.5406
        assert g.wavenumber == pytest.approx(
            wavelength_to_wavenumber(1.5406), rel=1e-10
        )

    def test_set_updates_wavelength(self):
        """Setting wavenumber updates wavelength via λ = 2π/k."""
        g = fourcv()
        k = 4.0
        g.wavenumber = k
        assert g.wavelength == pytest.approx(wavenumber_to_wavelength(k), rel=1e-10)

    def test_set_none_clears_wavelength(self):
        """Setting wavenumber = None clears wavelength."""
        g = fourcv()
        g.wavelength = 1.5406
        g.wavenumber = None
        assert g.wavelength is None
        assert g.wavenumber is None

    def test_round_trip_via_wavelength(self):
        """wavelength → wavenumber → wavelength recovers the original."""
        g = fourcv()
        g.wavelength = 1.5406
        k = g.wavenumber
        g2 = fourcv()
        g2.wavenumber = k
        assert g2.wavelength == pytest.approx(1.5406, rel=1e-10)


# ---------------------------------------------------------------------------
# summary() includes energy
# ---------------------------------------------------------------------------


def test_summary_includes_energy(capsys):
    """summary() reports energy (keV) alongside wavelength when set."""
    g = fourcv()
    g.wavelength = 1.5406
    g.summary()
    captured = capsys.readouterr()
    assert "keV" in captured.out


def test_summary_no_energy_when_wavelength_not_set(capsys):
    """summary() reports 'not set' when wavelength is None."""
    g = fourcv()
    g.summary()
    captured = capsys.readouterr()
    assert "not set" in captured.out
    assert "keV" not in captured.out


# ---------------------------------------------------------------------------
# Public API exports
# ---------------------------------------------------------------------------


def test_public_api_exports():
    """All radiation symbols are importable from the top-level package."""
    assert ahd.radiation.HC_KEV_ANGSTROM is HC_KEV_ANGSTROM
    assert ahd.radiation.HC_KEV_ANGSTROM_UNCERTAINTY is HC_KEV_ANGSTROM_UNCERTAINTY
    assert ahd.radiation.NEUTRON_MEV_ANGSTROM2 is NEUTRON_MEV_ANGSTROM2
    assert (
        ahd.radiation.NEUTRON_MEV_ANGSTROM2_UNCERTAINTY
        is NEUTRON_MEV_ANGSTROM2_UNCERTAINTY
    )
    assert ahd.radiation.SOURCE_TYPES is SOURCE_TYPES
    assert ahd.radiation.XRAY_LINES is XRAY_LINES
    assert ahd.radiation.wavelength_to_energy is wavelength_to_energy
    assert ahd.radiation.energy_to_wavelength is energy_to_wavelength
    assert ahd.radiation.wavelength_to_wavenumber is wavelength_to_wavenumber
    assert ahd.radiation.wavenumber_to_wavelength is wavenumber_to_wavelength
    assert ahd.radiation.neutron_wavelength_to_energy is neutron_wavelength_to_energy
    assert ahd.radiation.neutron_energy_to_wavelength is neutron_energy_to_wavelength


# ---------------------------------------------------------------------------
# Neutron constants
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
    """E = 25 meV → λ = sqrt(81.8042102352 / 25)."""
    assert neutron_energy_to_wavelength(25.0) == pytest.approx(
        math.sqrt(NEUTRON_MEV_ANGSTROM2 / 25.0), rel=1e-10
    )


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
        assert g_xray.energy == pytest.approx(8.0, abs=0.5)  # keV
        assert g_n.energy == pytest.approx(35.0, abs=5.0)  # meV

    def test_set_neutron_energy_updates_wavelength(self):
        """Setting energy (meV) for neutron updates wavelength via de Broglie."""
        g = fourcv()
        g.source_type = "neutron"
        g.energy = 25.0
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
        e_xray = g.energy
        g.source_type = "neutron"
        e_neutron = g.energy
        assert g.energy_units == "meV"
        assert e_xray != pytest.approx(e_neutron, rel=0.01)


# ---------------------------------------------------------------------------
# summary() — neutron source
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


# ---------------------------------------------------------------------------
# X-ray and neutron conversions do not mix
# ---------------------------------------------------------------------------


def test_xray_and_neutron_energies_differ_at_same_wavelength():
    """X-ray (keV) and neutron (meV) formulas give numerically different results."""
    wl = 1.5406
    e_xray = ahd.radiation.wavelength_to_energy(wl)  # ~8 keV
    e_neutron = ahd.radiation.neutron_wavelength_to_energy(wl)  # ~35 meV
    # Not equal even after ×1000 (different physical formulas)
    assert abs(e_xray * 1000 - e_neutron) > 1.0
