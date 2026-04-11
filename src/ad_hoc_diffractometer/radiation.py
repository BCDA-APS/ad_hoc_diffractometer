# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
radiation.py — Wavelength, energy, and wave-number conversions for X-rays.

Provides named constants for common laboratory X-ray emission lines and
standalone conversion functions between wavelength, photon energy (keV),
and wave number k (Å⁻¹).

Constants
---------
XRAY_LINES : dict[str, float]
    Named X-ray emission line wavelengths in Å.  Weighted means (Kα) or
    exact values (Kα1) as noted.

    ``"Cu_Ka"``   Cu Kα  ≈ 1.5406 Å  (weighted mean of Kα1=1.54056 Å, Kα2=1.54439 Å)
    ``"Cu_Ka1"``  Cu Kα1 = 1.54056 Å
    ``"Cu_Ka2"``  Cu Kα2 = 1.54439 Å
    ``"Mo_Ka"``   Mo Kα  ≈ 0.7107 Å  (weighted mean)
    ``"Mo_Ka1"``  Mo Kα1 = 0.70930 Å
    ``"Ag_Ka"``   Ag Kα  ≈ 0.5594 Å  (weighted mean)
    ``"Ag_Ka1"``  Ag Kα1 = 0.55941 Å
    ``"Co_Ka"``   Co Kα  ≈ 1.7902 Å  (weighted mean)
    ``"Co_Ka1"``  Co Kα1 = 1.78897 Å

Functions
---------
wavelength_to_energy(wavelength) -> float
    Wavelength (Å) → photon energy E (keV): E = hc/λ = 12.39842 / λ.

energy_to_wavelength(energy_kev) -> float
    Photon energy (keV) → wavelength (Å): λ = hc/E = 12.39842 / E.

wavelength_to_wavenumber(wavelength) -> float
    Wavelength (Å) → wave number k (Å⁻¹): k = 2π / λ.

wavenumber_to_wavelength(wavenumber) -> float
    Wave number k (Å⁻¹) → wavelength (Å): λ = 2π / k.

Notes
-----
The conversion constant hc = 12.39842 keV·Å is the standard value used
in X-ray crystallography (NIST CODATA 2018).

These conversions apply to **X-rays only**.  For neutrons, the de Broglie
relation gives E (meV) = 81.8042 / λ² (Å) — see the `radiation_neutron`
module (issue #8).

References
----------
NIST, "X-Ray Transition Energies Database",
https://physics.nist.gov/PhysRefData/XrayTrans/Html/search.html

NIST CODATA 2018 — hc = 12398.42 eV·Å = 12.39842 keV·Å
"""

from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Physical constant
# ---------------------------------------------------------------------------

#: hc in keV·Å (NIST CODATA 2018).
HC_KEV_ANGSTROM: float = 12.39842

# ---------------------------------------------------------------------------
# Named X-ray emission lines (wavelengths in Å)
# ---------------------------------------------------------------------------

#: Common laboratory X-ray emission line wavelengths in Å.
#:
#: Weighted means (Kα) use the standard 2:1 weighting of Kα1 and Kα2.
#: Values from NIST X-Ray Transition Energies Database.
XRAY_LINES: dict[str, float] = {
    "Cu_Ka": 1.54060,  # Cu Kα  weighted mean (2·Kα1 + Kα2) / 3
    "Cu_Ka1": 1.54056,  # Cu Kα1
    "Cu_Ka2": 1.54439,  # Cu Kα2
    "Mo_Ka": 0.71073,  # Mo Kα  weighted mean
    "Mo_Ka1": 0.70930,  # Mo Kα1
    "Ag_Ka": 0.56087,  # Ag Kα  weighted mean
    "Ag_Ka1": 0.55941,  # Ag Kα1
    "Co_Ka": 1.79021,  # Co Kα  weighted mean
    "Co_Ka1": 1.78897,  # Co Kα1
}

# ---------------------------------------------------------------------------
# Conversion functions
# ---------------------------------------------------------------------------


def wavelength_to_energy(wavelength: float) -> float:
    """
    Convert wavelength to photon energy for X-rays.

    E (keV) = hc / λ = 12.39842 / λ (Å)

    Parameters
    ----------
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        Photon energy in keV.

    Raises
    ------
    ValueError
        If ``wavelength`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.wavelength_to_energy(1.5406), 4)
    8.0478
    """
    if wavelength <= 0.0:
        raise ValueError(
            f"wavelength_to_energy(): wavelength must be > 0 Å; got {wavelength}."
        )
    return HC_KEV_ANGSTROM / wavelength


def energy_to_wavelength(energy_kev: float) -> float:
    """
    Convert photon energy to wavelength for X-rays.

    λ (Å) = hc / E = 12.39842 / E (keV)

    Parameters
    ----------
    energy_kev : float
        Photon energy in keV.

    Returns
    -------
    float
        Wavelength in Å.

    Raises
    ------
    ValueError
        If ``energy_kev`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.energy_to_wavelength(8.047), 4)
    1.5408
    """
    if energy_kev <= 0.0:
        raise ValueError(
            f"energy_to_wavelength(): energy must be > 0 keV; got {energy_kev}."
        )
    return HC_KEV_ANGSTROM / energy_kev


def wavelength_to_wavenumber(wavelength: float) -> float:
    """
    Convert wavelength to wave number k.

    k (Å⁻¹) = 2π / λ (Å)

    Parameters
    ----------
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        Wave number k in Å⁻¹.

    Raises
    ------
    ValueError
        If ``wavelength`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.wavelength_to_wavenumber(1.0), 6)
    6.283185
    """
    if wavelength <= 0.0:
        raise ValueError(
            f"wavelength_to_wavenumber(): wavelength must be > 0 Å; got {wavelength}."
        )
    return 2.0 * math.pi / wavelength


def wavenumber_to_wavelength(wavenumber: float) -> float:
    """
    Convert wave number k to wavelength.

    λ (Å) = 2π / k (Å⁻¹)

    Parameters
    ----------
    wavenumber : float
        Wave number k in Å⁻¹.

    Returns
    -------
    float
        Wavelength in Å.

    Raises
    ------
    ValueError
        If ``wavenumber`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.wavenumber_to_wavelength(6.283185), 6)
    1.0
    """
    if wavenumber <= 0.0:
        raise ValueError(
            f"wavenumber_to_wavelength(): wavenumber must be > 0 Å⁻¹; got {wavenumber}."
        )
    return 2.0 * math.pi / wavenumber
