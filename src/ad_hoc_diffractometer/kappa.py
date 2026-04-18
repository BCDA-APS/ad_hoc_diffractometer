# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
kappa.py — Kappa-to-Eulerian angle conversion.

Kappa diffractometers replace the Eulerian chi circle with a kappa arm
tilted at angle alpha_0 from the omega axis.  The real kappa motor angles
(komega, kappa, kphi) map to virtual Eulerian pseudoangles (omega, chi, phi)
via the relations in Walko (2016), eq. [16]:

    chi   = 2 arcsin[sin(kappa/2) · sin(alpha_0)]
    offset = arccos[cos(kappa/2) / cos(chi/2)]
    omega = komega − offset
    phi   = kphi   − offset

Both directions are implemented here as pure-NumPy functions.

Functions
---------
:func:`kappa_to_eulerian`
    Forward direction: real (komega, kappa, kphi) → virtual (omega, chi, phi).

:func:`eulerian_to_kappa`
    Inverse direction: virtual (omega, chi, phi) → real (komega, kappa, kphi).

Notes
-----
The default kappa tilt angle is 50° (Walko 2016; Enraf-Nonius convention;
ITC Vol. C §2.2.6).

The inverse has two solutions distinguished by the sign of kappa (the *branch*
parameter).  The positive branch (kappa ≥ 0) is the default; the negative
branch gives kappa ≤ 0.  The caller should choose the branch that keeps all
angles within the hardware limits.

Singularity: when ``|chi|`` approaches ±2·alpha_0, ``|kappa|`` approaches 180° and
the offset formula becomes numerically ill-conditioned.  A :class:`ValueError`
is raised when ``|chi|`` ≥ 2·alpha_0 − epsilon.

References
----------
* Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), eq. [16].
* ITC Vol. C §2.2.6 (2006).
"""

from __future__ import annotations

import math


def kappa_to_eulerian(
    komega: float,
    kappa: float,
    kphi: float,
    alpha_deg: float = 50.0,
) -> tuple[float, float, float]:
    """
    Convert real kappa angles to virtual Eulerian pseudoangles.

    Implements Walko (2016), eq. [16] — forward direction.

    Parameters
    ----------
    komega : float
        Real kappa-omega angle in degrees.
    kappa : float
        Real kappa angle in degrees.
    kphi : float
        Real kappa-phi angle in degrees.
    alpha_deg : float, optional
        Kappa tilt angle in degrees.  Default 50°.

    Returns
    -------
    omega, chi, phi : tuple of float
        Virtual Eulerian pseudoangles in degrees.

    Raises
    ------
    ValueError
        If the combination of kappa and alpha_deg places chi outside
        the reachable range (``|kappa/2|`` · sin(alpha_0) > 1).

    Examples
    --------
    >>> kappa_to_eulerian(57.045, 134.756, 57.045, alpha_deg=50.0)
    (-0.0, 90.0, -0.0)
    """
    a0 = math.radians(alpha_deg)
    k = math.radians(kappa)

    sin_chi_half = math.sin(k / 2.0) * math.sin(a0)
    if abs(sin_chi_half) > 1.0:  # pragma: no cover
        raise ValueError(
            f"kappa_to_eulerian: kappa={kappa:.4f}°, alpha_deg={alpha_deg:.4f}° "
            f"gives |sin(chi/2)| = {abs(sin_chi_half):.6f} > 1 — unreachable chi."
        )

    chi_half = math.asin(sin_chi_half)
    chi = math.degrees(2.0 * chi_half)

    cos_chi_half = math.cos(chi_half)
    if abs(cos_chi_half) < 1e-14:  # pragma: no cover
        raise ValueError(
            f"kappa_to_eulerian: chi approaches ±180° (chi_half={math.degrees(chi_half):.4f}°) "
            f"— singularity in the offset calculation."
        )

    cos_ratio = math.cos(k / 2.0) / cos_chi_half
    cos_ratio = max(-1.0, min(1.0, cos_ratio))  # guard against float rounding
    offset = math.degrees(math.acos(cos_ratio))
    # The sign of offset matches the sign of kappa
    if kappa < 0.0:
        offset = -offset

    omega = komega - offset
    phi = kphi - offset

    return omega, chi, phi


def eulerian_to_kappa(
    omega: float,
    chi: float,
    phi: float,
    alpha_deg: float = 50.0,
    branch: int = +1,
) -> tuple[float, float, float]:
    """
    Convert virtual Eulerian pseudoangles to real kappa angles.

    Inverse of Walko (2016), eq. [16].

    Parameters
    ----------
    omega : float
        Virtual Eulerian omega in degrees.
    chi : float
        Virtual Eulerian chi in degrees.
    phi : float
        Virtual Eulerian phi in degrees.
    alpha_deg : float, optional
        Kappa tilt angle in degrees.  Default 50°.
    branch : {+1, -1}, optional
        Branch selection.  ``+1`` gives kappa ≥ 0 (default);
        ``-1`` gives kappa ≤ 0.

    Returns
    -------
    komega, kappa, kphi : tuple of float
        Real kappa motor angles in degrees.

    Raises
    ------
    ValueError
        If ``|chi|`` ≥ 2·alpha_deg (chi outside the reachable range).
    ValueError
        If branch is not +1 or -1.

    Examples
    --------
    >>> eulerian_to_kappa(0.0, 90.0, 0.0, alpha_deg=50.0)
    (57.04519..., 134.75594..., 57.04519...)
    >>> eulerian_to_kappa(0.0, 90.0, 0.0, alpha_deg=50.0, branch=-1)
    (-57.04519..., -134.75594..., -57.04519...)
    """
    if branch not in (+1, -1):
        raise ValueError(f"eulerian_to_kappa: branch must be +1 or -1; got {branch!r}.")

    a0 = math.radians(alpha_deg)
    chi_half = math.radians(chi / 2.0)

    chi_limit = 2.0 * alpha_deg
    if abs(chi) >= chi_limit - 1e-9:
        raise ValueError(
            f"eulerian_to_kappa: |chi| = {abs(chi):.4f}° must be less than "
            f"2·alpha_deg = {chi_limit:.4f}° (kappa geometry limit).  "
            "The requested chi is outside the reachable range."
        )

    # Use |chi/2| to get the magnitude of kappa, then apply branch for sign.
    sin_kappa_half = abs(math.sin(chi_half)) / math.sin(a0)
    # Clamp for float safety (already guarded above, but defensive)
    sin_kappa_half = max(-1.0, min(1.0, sin_kappa_half))
    kappa_half = math.asin(sin_kappa_half)
    kappa = math.degrees(2.0 * kappa_half) * branch

    # Recompute chi_half from the chosen kappa branch for consistency
    k = math.radians(kappa)
    sin_ch = math.sin(k / 2.0) * math.sin(a0)
    ch = math.asin(max(-1.0, min(1.0, sin_ch)))

    cos_ch = math.cos(ch)
    if abs(cos_ch) < 1e-14:  # pragma: no cover
        raise ValueError("eulerian_to_kappa: chi approaches ±180° — singularity.")

    cos_ratio = math.cos(k / 2.0) / cos_ch
    cos_ratio = max(-1.0, min(1.0, cos_ratio))
    offset = math.degrees(math.acos(cos_ratio))
    if kappa < 0.0:
        offset = -offset

    komega = omega + offset
    kphi = phi + offset

    return komega, kappa, kphi
