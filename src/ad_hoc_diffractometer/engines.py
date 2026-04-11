# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
engines.py — Alternative calculation engines.

Provides standalone conversion functions between the four representations
of a diffraction position:

- **hkl** — Miller indices (h, k, l); dimensionless integers or reals
- **Q-vector** — scattering vector (Qx, Qy, Qz) in Å⁻¹ (reciprocal Cartesian)
- **d-spacing** — interplanar spacing d in Å; scalar
- **two-theta** — diffraction angle 2θ in degrees; requires wavelength

All conversions follow Busing & Levy (1967) and You (1999).

Functions
---------
hkl_to_Q(sample, h, k, l) -> np.ndarray
    Miller indices → Q-vector (Å⁻¹) via Q = UB @ hkl.

Q_to_hkl(sample, Qx, Qy, Qz) -> tuple[float, float, float]
    Q-vector → Miller indices via hkl = UB⁻¹ @ Q.

Q_to_d(Qx, Qy, Qz) -> float
    Q-vector magnitude → d-spacing: d = 2π / |Q|.

d_to_Q_mag(d) -> float
    d-spacing → |Q|: |Q| = 2π / d.

hkl_to_d(sample, h, k, l) -> float
    Miller indices → d-spacing via |Q| = |UB @ hkl|, d = 2π / |Q|.

d_to_two_theta(d, wavelength) -> float
    d-spacing → 2θ (degrees) via Bragg's law: λ = 2d sin(θ).

two_theta_to_d(two_theta_deg, wavelength) -> float
    2θ (degrees) → d-spacing via Bragg's law.

hkl_to_two_theta(sample, h, k, l, wavelength) -> float
    Miller indices → 2θ (degrees).

two_theta_to_Q_mag(two_theta_deg, wavelength) -> float
    2θ (degrees) → |Q| (Å⁻¹): |Q| = 4π sin(θ) / λ.

Q_mag_to_two_theta(Q_mag, wavelength) -> float
    |Q| (Å⁻¹) → 2θ (degrees).

References
----------
Busing & Levy, Acta Cryst. 22, 457-464 (1967) — B matrix, Q = UB @ hkl.
You, J. Appl. Cryst. 32, 614-623 (1999) — Bragg condition, 2θ definition.
"""

from __future__ import annotations

import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .sample import Sample

_TWO_PI = 2.0 * math.pi


# ---------------------------------------------------------------------------
# hkl ↔ Q-vector
# ---------------------------------------------------------------------------


def hkl_to_Q(
    sample: Sample,
    h: float,
    k: float,
    l: float,  # noqa: E741
) -> np.ndarray:
    """
    Convert Miller indices to the scattering vector Q in Å⁻¹.

    Uses the UB matrix: **Q** = UB · **hkl**.

    Parameters
    ----------
    sample : Sample
        Sample with a UB matrix set.
    h, k, l : float
        Miller indices.

    Returns
    -------
    numpy.ndarray, shape (3,)
        Scattering vector (Qx, Qy, Qz) in Å⁻¹.

    Raises
    ------
    ValueError
        If ``sample.UB`` is None.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.sample.lattice = ahd.Lattice(a=3.905)
    >>> ahd.ub_identity(g.sample)
    >>> Q = ahd.hkl_to_Q(g.sample, 1, 0, 0)
    >>> float(Q[0])  # doctest: +ELLIPSIS
    1.6...
    """
    if sample.UB is None:
        raise ValueError(
            "hkl_to_Q() requires sample.UB to be set. "
            "Call ub_identity() or ub_from_*() first."
        )
    return sample.UB @ np.array([float(h), float(k), float(l)], dtype=float)


def Q_to_hkl(
    sample: Sample,
    Qx: float,
    Qy: float,
    Qz: float,
) -> tuple[float, float, float]:
    """
    Convert a Q-vector (Å⁻¹) to Miller indices.

    Solves **hkl** = UB⁻¹ · **Q**.

    Parameters
    ----------
    sample : Sample
        Sample with a UB matrix set.
    Qx, Qy, Qz : float
        Scattering vector components in Å⁻¹.

    Returns
    -------
    tuple of float
        Miller indices (h, k, l).

    Raises
    ------
    ValueError
        If ``sample.UB`` is None or singular.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> import numpy as np
    >>> g = ahd.fourcv()
    >>> g.sample.lattice = ahd.Lattice(a=3.905)
    >>> ahd.ub_identity(g.sample)
    >>> Q = ahd.hkl_to_Q(g.sample, 1, 0, 0)
    >>> ahd.Q_to_hkl(g.sample, *Q)  # doctest: +ELLIPSIS
    (1.0, 0.0, 0.0)
    """
    if sample.UB is None:
        raise ValueError(
            "Q_to_hkl() requires sample.UB to be set. "
            "Call ub_identity() or ub_from_*() first."
        )
    Q = np.array([float(Qx), float(Qy), float(Qz)], dtype=float)
    try:
        hkl = np.linalg.solve(sample.UB, Q)
    except np.linalg.LinAlgError as exc:
        raise ValueError("Q_to_hkl(): UB matrix is singular; cannot invert.") from exc
    return (float(hkl[0]), float(hkl[1]), float(hkl[2]))


# ---------------------------------------------------------------------------
# Q-magnitude ↔ d-spacing
# ---------------------------------------------------------------------------


def Q_to_d(Qx: float, Qy: float, Qz: float) -> float:
    """
    Convert a Q-vector to the d-spacing.

    d = 2π / |**Q**|

    Parameters
    ----------
    Qx, Qy, Qz : float
        Scattering vector components in Å⁻¹.

    Returns
    -------
    float
        d-spacing in Å.

    Raises
    ------
    ValueError
        If |Q| = 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.Q_to_d(1.0, 0.0, 0.0), 6)
    6.283185
    """
    Q_mag = math.sqrt(Qx**2 + Qy**2 + Qz**2)
    if Q_mag < 1e-14:
        raise ValueError("Q_to_d(): |Q| is zero; d-spacing is undefined.")
    return _TWO_PI / Q_mag


def d_to_Q_mag(d: float) -> float:
    """
    Convert a d-spacing to the magnitude of the scattering vector.

    |**Q**| = 2π / d

    Parameters
    ----------
    d : float
        d-spacing in Å.

    Returns
    -------
    float
        |Q| in Å⁻¹.

    Raises
    ------
    ValueError
        If ``d`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.d_to_Q_mag(1.0), 6)
    6.283185
    """
    if d <= 0.0:
        raise ValueError(f"d_to_Q_mag(): d must be > 0 Å; got {d}.")
    return _TWO_PI / d


def hkl_to_d(
    sample: Sample,
    h: float,
    k: float,
    l: float,  # noqa: E741
) -> float:
    """
    Convert Miller indices to d-spacing.

    d = 2π / |UB · **hkl**|

    Parameters
    ----------
    sample : Sample
        Sample with a UB matrix set.
    h, k, l : float
        Miller indices.

    Returns
    -------
    float
        d-spacing in Å.

    Raises
    ------
    ValueError
        If ``sample.UB`` is None or if **hkl** = (0, 0, 0).

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.sample.lattice = ahd.Lattice(a=3.905)
    >>> ahd.ub_identity(g.sample)
    >>> round(ahd.hkl_to_d(g.sample, 1, 0, 0), 4)
    3.905
    """
    Q = hkl_to_Q(sample, h, k, l)
    Q_mag = float(np.linalg.norm(Q))
    if Q_mag < 1e-14:
        raise ValueError(
            f"hkl_to_d(): |UB @ ({h}, {k}, {l})| is zero; "
            "d-spacing is undefined for the (0, 0, 0) reflection."
        )
    return _TWO_PI / Q_mag


# ---------------------------------------------------------------------------
# d-spacing ↔ two-theta  (Bragg's law)
# ---------------------------------------------------------------------------


def d_to_two_theta(d: float, wavelength: float) -> float:
    """
    Convert d-spacing to 2θ using Bragg's law.

    λ = 2 d sin(θ)  →  2θ = 2 arcsin(λ / (2d))

    Parameters
    ----------
    d : float
        d-spacing in Å.
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        2θ in degrees.

    Raises
    ------
    ValueError
        If ``d`` ≤ 0, ``wavelength`` ≤ 0, or λ / (2d) > 1 (reflection
        unreachable at this wavelength).

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.d_to_two_theta(2.0, 1.5406), 4)
    22.6956
    """
    if d <= 0.0:
        raise ValueError(f"d_to_two_theta(): d must be > 0 Å; got {d}.")
    if wavelength <= 0.0:
        raise ValueError(
            f"d_to_two_theta(): wavelength must be > 0 Å; got {wavelength}."
        )
    sin_theta = wavelength / (2.0 * d)
    if sin_theta > 1.0:
        raise ValueError(
            f"d_to_two_theta(): reflection unreachable at wavelength {wavelength} Å "
            f"(λ/(2d) = {sin_theta:.4f} > 1). Use a shorter wavelength or larger d."
        )
    return 2.0 * math.degrees(math.asin(sin_theta))


def two_theta_to_d(two_theta_deg: float, wavelength: float) -> float:
    """
    Convert 2θ to d-spacing using Bragg's law.

    λ = 2 d sin(θ)  →  d = λ / (2 sin(θ))

    Parameters
    ----------
    two_theta_deg : float
        Diffraction angle 2θ in degrees.
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        d-spacing in Å.

    Raises
    ------
    ValueError
        If ``two_theta_deg`` is not in (0°, 180°) or ``wavelength`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.two_theta_to_d(22.6956, 1.5406), 4)
    2.0
    """
    if wavelength <= 0.0:
        raise ValueError(
            f"two_theta_to_d(): wavelength must be > 0 Å; got {wavelength}."
        )
    if not (0.0 < two_theta_deg < 180.0):
        raise ValueError(
            f"two_theta_to_d(): 2θ must be in (0°, 180°); got {two_theta_deg}."
        )
    sin_theta = math.sin(math.radians(two_theta_deg / 2.0))
    return wavelength / (2.0 * sin_theta)


def hkl_to_two_theta(
    sample: Sample,
    h: float,
    k: float,
    l: float,  # noqa: E741
    wavelength: float,
) -> float:
    """
    Convert Miller indices to 2θ using the UB matrix and Bragg's law.

    Parameters
    ----------
    sample : Sample
        Sample with a UB matrix set.
    h, k, l : float
        Miller indices.
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        2θ in degrees.

    Raises
    ------
    ValueError
        If ``sample.UB`` is None, ``wavelength`` ≤ 0, **hkl** = (0, 0, 0),
        or the reflection is unreachable at this wavelength.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.sample.lattice = ahd.Lattice(a=3.905)
    >>> ahd.ub_identity(g.sample)
    >>> round(ahd.hkl_to_two_theta(g.sample, 1, 0, 0, 1.5406), 4)
    23.2066
    """
    d = hkl_to_d(sample, h, k, l)
    return d_to_two_theta(d, wavelength)


# ---------------------------------------------------------------------------
# |Q| ↔ two-theta
# ---------------------------------------------------------------------------


def two_theta_to_Q_mag(two_theta_deg: float, wavelength: float) -> float:
    """
    Convert 2θ to |Q| using the Bragg relation.

    |**Q**| = 4π sin(θ) / λ

    Parameters
    ----------
    two_theta_deg : float
        Diffraction angle 2θ in degrees.
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        |Q| in Å⁻¹.

    Raises
    ------
    ValueError
        If ``two_theta_deg`` is not in (0°, 180°) or ``wavelength`` ≤ 0.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.two_theta_to_Q_mag(30.0, 1.5406), 4)
    2.0904
    """
    if wavelength <= 0.0:
        raise ValueError(
            f"two_theta_to_Q_mag(): wavelength must be > 0 Å; got {wavelength}."
        )
    if not (0.0 < two_theta_deg < 180.0):
        raise ValueError(
            f"two_theta_to_Q_mag(): 2θ must be in (0°, 180°); got {two_theta_deg}."
        )
    return 2.0 * _TWO_PI * math.sin(math.radians(two_theta_deg / 2.0)) / wavelength


def Q_mag_to_two_theta(Q_mag: float, wavelength: float) -> float:
    """
    Convert |Q| to 2θ.

    2θ = 2 arcsin(|Q| λ / (4π))

    Parameters
    ----------
    Q_mag : float
        Magnitude of the scattering vector in Å⁻¹.
    wavelength : float
        Wavelength in Å.

    Returns
    -------
    float
        2θ in degrees.

    Raises
    ------
    ValueError
        If ``Q_mag`` < 0, ``wavelength`` ≤ 0, or the reflection is
        unreachable (|Q| λ / (4π) > 1).

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> round(ahd.Q_mag_to_two_theta(2.0904, 1.5406), 2)
    30.0
    """
    if wavelength <= 0.0:
        raise ValueError(
            f"Q_mag_to_two_theta(): wavelength must be > 0 Å; got {wavelength}."
        )
    if Q_mag < 0.0:
        raise ValueError(f"Q_mag_to_two_theta(): |Q| must be ≥ 0; got {Q_mag}.")
    sin_theta = Q_mag * wavelength / (2.0 * _TWO_PI)
    if sin_theta > 1.0:
        raise ValueError(
            f"Q_mag_to_two_theta(): reflection unreachable (|Q|λ/(4π) = "
            f"{sin_theta:.4f} > 1). Use a longer wavelength or smaller |Q|."
        )
    return 2.0 * math.degrees(math.asin(sin_theta))
