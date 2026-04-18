# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
radiation.py — Wavelength, energy, and wave-number conversions.

Provides named constants for common laboratory X-ray emission lines,
standalone conversion functions for both X-rays and fixed-wavelength
neutrons (reactor sources), and the ``SOURCE_TYPES`` sentinel used by
:class:`AdHocDiffractometer` to select the correct energy formula.

**Source type scope**

This module supports two source types:

``"xray"`` (default)
    Photon energy: E (keV) = hc/λ = 12.39842 / λ (Å).

``"neutron"``
    Fixed-wavelength **reactor** neutrons only.  Energy via the de Broglie
    relation: E (meV) = h²/(2 m_n λ²) = 81.8042 / λ² (Å).
    Spallation/time-of-flight sources sweep the full spectrum with each
    pulse and require a fundamentally different instrument description;
    they are out of scope for this project.

Constants
---------
HC_KEV_ANGSTROM : float
    hc = 12.398 419 843 320 026 keV·Å.  Exact (NIST CODATA 2022).

HC_KEV_ANGSTROM_UNCERTAINTY : float
    0.0 — exactly zero.  The 2019 SI redefinition fixed both *h* and *c*
    exactly, so *hc* has no uncertainty by definition.

NEUTRON_MEV_ANGSTROM2 : float
    h²/(2 m_n) = 81.804 210 235 2 meV·Å² (NIST CODATA 2022).

NEUTRON_MEV_ANGSTROM2_UNCERTAINTY : float
    0.0000000415 meV·Å² — propagated from CODATA 2022 m_n uncertainty.

SOURCE_TYPES : tuple[str, ...]
    Valid source-type strings: ``("xray", "neutron")``.

XRAY_LINES : dict[str, float]
    Named X-ray emission line wavelengths in Å.

X-ray functions
---------------
wavelength_to_energy(wavelength) -> float
    λ (Å) → E (keV) = hc/λ.

energy_to_wavelength(energy_kev) -> float
    E (keV) → λ (Å) = hc/E.

wavelength_to_wavenumber(wavelength) -> float
    λ (Å) → k (Å⁻¹) = 2π/λ.

wavenumber_to_wavelength(wavenumber) -> float
    k (Å⁻¹) → λ (Å) = 2π/k.

Neutron functions
-----------------
neutron_wavelength_to_energy(wavelength) -> float
    λ (Å) → E (meV) = 81.8042 / λ².

neutron_energy_to_wavelength(energy_mev) -> float
    E (meV) → λ (Å) = sqrt(81.8042 / E).

References
----------
NIST CODATA 2018 — hc = 12398.42 eV·Å = 12.39842 keV·Å
NIST neutron scattering — h²/(2 m_n) = 81.8042 meV·Å²
de Broglie relation — λ = h/p, E = p²/(2 m_n)
"""

from __future__ import annotations

import math

# ---------------------------------------------------------------------------
# Physical constant
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# Physical constants — NIST CODATA 2022 adjustment
# https://physics.nist.gov/cuu/Constants/Table/allascii.txt
#
# HC_KEV_ANGSTROM: hc in keV·Å.
#   Derived from h (exact), c (exact), e (exact) per SI 2019 redefinition.
#   hc = 6.62607015e-34 J·Hz⁻¹ × 299792458 m·s⁻¹ / 1.602176634e-19 J·eV⁻¹
#      = 1.239841984332003e-6 eV·m  (exact, listed in CODATA table)
#      = 12.398419843320026 keV·Å  (exact — no uncertainty)
#
# NEUTRON_MEV_ANGSTROM2: h²/(2 m_n) in meV·Å².
#   h (exact); m_n = 1.674 927 500 56 e-27 kg ± 0.000 000 000 85 e-27 kg
#   CODATA 2022 relative uncertainty of m_n: 5.075e-10
#   NEUTRON_MEV_ANGSTROM2_UNCERTAINTY = 0.0000000415 meV·Å²
# ---------------------------------------------------------------------------

HC_KEV_ANGSTROM: float = 12.398419843320026
"""hc = 12.398 419 843 320 026 keV·Å.

Exact (h, c, and e are all defined SI constants since 2019).
Source: NIST CODATA 2022.
"""

HC_KEV_ANGSTROM_UNCERTAINTY: float = 0.0
"""Uncertainty of :data:`HC_KEV_ANGSTROM` in keV·Å.

Exactly **0.0** — not a placeholder.  The 2019 redefinition of the SI
fixed the numerical values of Planck's constant *h* and the speed of light
*c* exactly; their product *hc* therefore has zero uncertainty by
definition.  Reference: BIPM SI Brochure, 9th edition (2019); NIST CODATA 2022.
"""

NEUTRON_MEV_ANGSTROM2: float = 81.8042102352
"""h²/(2 m_n) = 81.804 210 235 2 meV·Å².

Used for reactor-neutron kinetic energy via the de Broglie relation.
Source: NIST CODATA 2022; m_n = 1.674 927 500 56(85) × 10⁻²⁷ kg.
"""

NEUTRON_MEV_ANGSTROM2_UNCERTAINTY: float = 0.0000000415
"""Uncertainty of :data:`NEUTRON_MEV_ANGSTROM2` in meV·Å².

Propagated from the CODATA 2022 uncertainty in the neutron mass m_n.
"""

SOURCE_TYPES: tuple[str, ...] = ("xray", "neutron")
"""Valid source-type strings for :attr:`AdHocDiffractometer.source_type`.

``"xray"`` uses hc/λ (keV); ``"neutron"`` uses de Broglie h²/(2 m_n λ²) (meV).
"""

# ---------------------------------------------------------------------------
# Named X-ray emission lines (wavelengths in Å)
# ---------------------------------------------------------------------------

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
"""Common laboratory X-ray emission line wavelengths in Å.

Weighted means (Kα) use the standard 2:1 weighting of Kα1 and Kα2.
Values from NIST X-Ray Transition Energies Database.
"""

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


# ---------------------------------------------------------------------------
# Neutron conversion functions  (de Broglie, reactor / fixed-wavelength only)
# ---------------------------------------------------------------------------


def neutron_wavelength_to_energy(wavelength: float) -> float:
    """
    Convert wavelength to neutron kinetic energy (de Broglie relation).

    E (meV) = h² / (2 m_n λ²) = 81.8042 / λ² (Å)

    Applies to **fixed-wavelength reactor neutrons** only.  Spallation
    (time-of-flight) sources are out of scope for this project.

    Parameters
    ----------
    wavelength : float
        Neutron wavelength in Å.

    Returns
    -------
    float
        Neutron kinetic energy in meV.

    Raises
    ------
    ValueError
        If ``wavelength`` ≤ 0.

    See Also
    --------
    wavelength_to_energy : X-ray photon energy in keV.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.neutron_wavelength_to_energy(1.8), 4)
    25.2482
    """
    if wavelength <= 0.0:
        raise ValueError(
            "neutron_wavelength_to_energy(): wavelength must be > 0 Å; "
            f"got {wavelength}."
        )
    return NEUTRON_MEV_ANGSTROM2 / (wavelength**2)


def neutron_energy_to_wavelength(energy_mev: float) -> float:
    """
    Convert neutron kinetic energy to wavelength (de Broglie relation).

    λ (Å) = sqrt(h² / (2 m_n E)) = sqrt(81.8042 / E (meV))

    Applies to **fixed-wavelength reactor neutrons** only.

    Parameters
    ----------
    energy_mev : float
        Neutron kinetic energy in meV.

    Returns
    -------
    float
        Neutron wavelength in Å.

    Raises
    ------
    ValueError
        If ``energy_mev`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.neutron_energy_to_wavelength(25.0), 4)
    1.809
    """
    if energy_mev <= 0.0:
        raise ValueError(
            f"neutron_energy_to_wavelength(): energy must be > 0 meV; got {energy_mev}."
        )
    return math.sqrt(NEUTRON_MEV_ANGSTROM2 / energy_mev)
