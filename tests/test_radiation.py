# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Unit tests for ad_hoc_diffractometer.radiation (issue #21).

Covers:
  - HC_KEV_ANGSTROM constant value
  - XRAY_LINES: all keys present, values positive, Cu Kα known value
  - wavelength_to_energy: known value, round-trip, zero/negative raises
  - energy_to_wavelength: round-trip, zero/negative raises
  - wavelength_to_wavenumber: known value, round-trip, zero/negative raises
  - wavenumber_to_wavelength: round-trip, zero/negative raises
  - AdHocDiffractometer.energy_kev property: get, set, None
  - AdHocDiffractometer.wavenumber property: get, set, None
  - summary() includes energy when wavelength is set
  - Public API exports
"""

import math
import re
from contextlib import nullcontext as does_not_raise

import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer import fourcv
from ad_hoc_diffractometer.radiation import HC_KEV_ANGSTROM
from ad_hoc_diffractometer.radiation import XRAY_LINES
from ad_hoc_diffractometer.radiation import energy_to_wavelength
from ad_hoc_diffractometer.radiation import wavelength_to_energy
from ad_hoc_diffractometer.radiation import wavelength_to_wavenumber
from ad_hoc_diffractometer.radiation import wavenumber_to_wavelength

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------


def test_HC_KEV_ANGSTROM_value():
    """hc = 12.39842 keV·Å (NIST CODATA 2018)."""
    assert HC_KEV_ANGSTROM == pytest.approx(12.39842, rel=1e-5)


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
    "energy_kev, context",
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
def test_energy_to_wavelength(energy_kev, context):
    with context:
        result = energy_to_wavelength(energy_kev)
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
# AdHocDiffractometer.energy_kev property
# ---------------------------------------------------------------------------


class TestEnergyKevProperty:
    def test_none_when_wavelength_not_set(self):
        """energy_kev is None when wavelength is not set."""
        g = fourcv()
        assert g.energy_kev is None

    def test_value_matches_conversion(self):
        """energy_kev = wavelength_to_energy(wavelength)."""
        g = fourcv()
        g.wavelength = 1.5406
        assert g.energy_kev == pytest.approx(wavelength_to_energy(1.5406), rel=1e-10)

    def test_set_updates_wavelength(self):
        """Setting energy_kev updates wavelength via λ = hc/E."""
        g = fourcv()
        g.energy_kev = 8.048
        assert g.wavelength == pytest.approx(energy_to_wavelength(8.048), rel=1e-10)

    def test_set_none_clears_wavelength(self):
        """Setting energy_kev = None clears wavelength."""
        g = fourcv()
        g.wavelength = 1.5406
        g.energy_kev = None
        assert g.wavelength is None
        assert g.energy_kev is None

    def test_round_trip_via_wavelength(self):
        """wavelength → energy_kev → wavelength recovers the original."""
        g = fourcv()
        g.wavelength = 1.5406
        e = g.energy_kev
        g2 = fourcv()
        g2.energy_kev = e
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
    assert ahd.HC_KEV_ANGSTROM is HC_KEV_ANGSTROM
    assert ahd.XRAY_LINES is XRAY_LINES
    assert ahd.wavelength_to_energy is wavelength_to_energy
    assert ahd.energy_to_wavelength is energy_to_wavelength
    assert ahd.wavelength_to_wavenumber is wavelength_to_wavenumber
    assert ahd.wavenumber_to_wavelength is wavenumber_to_wavelength
