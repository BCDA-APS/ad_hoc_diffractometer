# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
surface.py — Surface geometry calculations.

Provides incidence angle αᵢ, emergence angle αf, and Q-vector decomposition
into in-plane (Q‖) and out-of-plane (Q⊥) components relative to the sample
surface.

All functions operate on a geometry instance and a set of motor angles
(defaulting to the current stage angles when ``angles`` is ``None``).

Public API
----------
incidence(geometry, angles=None)
    Angle of incidence αᵢ in degrees — angle between the incoming beam
    and the sample surface.  Accepts the deprecated alias ``alpha_i``.

emergence(geometry, angles=None)
    Angle of emergence αf in degrees — angle between the diffracted beam
    and the sample surface.  Accepts the deprecated alias ``alpha_f``.

q_components(geometry, angles=None)
    Decompose Q into Q⊥ (perpendicular to sample surface) and Q‖
    (parallel to sample surface) in Å⁻¹.  Returns a named dict.

is_specular(geometry, angles=None, atol=0.01)
    Return True when αᵢ ≈ αf within ``atol`` degrees.

is_evanescent(geometry, angles=None, critical_angle_deg=None)
    Return True when αᵢ < critical_angle_deg (evanescent / surface-
    sensitive regime).

Physical convention
-------------------
The sample surface normal is specified as Miller indices ``(h, k, l)``
via the geometry's ``surface_normal`` property (added to
``AdHocDiffractometer``).  If ``surface_normal`` is ``None`` the
``azimuth`` is used as a fallback.

The surface normal in the lab frame is computed as::

    n_phi = UB @ n_hkl          (phi-frame normal, same units as Q)
    n_lab = Z @ n_phi            (Z = sample rotation matrix)
    n̂_lab = n_lab / |n_lab|

The incident beam direction in the lab frame is the ``'longitudinal'``
basis vector (``+ŷ``) — the beam travels along +y.

The diffracted beam direction in the lab frame is ``D @ ŷ`` where D is
the total detector rotation matrix.

Then::

    sin(αᵢ) = ŷ · n̂_lab    (ŷ is incident direction, n̂_lab is surface normal)
    sin(αf) = (D @ ŷ) · n̂_lab

Q decomposition::

    Q_lab = (2π/λ) (D @ ŷ − ŷ)
    Q_perp_vec = (Q_lab · n̂_lab) n̂_lab
    Q_par_vec  = Q_lab − Q_perp_vec
    Q_perp = |Q_perp_vec|   (signed: positive if Q points outward)
    Q_par  = |Q_par_vec|

References
----------
* Lohmeier & Vlieg, J. Appl. Cryst. 26, 706-716 (1993) §4.2
* You, J. Appl. Cryst. 32, 614-623 (1999), eqs. 10-11
* Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987)
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _surface_vectors(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """
    Compute the vectors needed for surface angle calculations.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    angles : dict or None

    Returns
    -------
    y_hat : np.ndarray (3,)
        Unit incident-beam direction in the lab frame (longitudinal).
    kf_hat : np.ndarray (3,)
        Unit diffracted-beam direction in the lab frame.
    n_hat : np.ndarray (3,)
        Unit surface normal in the lab frame.
    Q_lab : np.ndarray (3,)
        Scattering vector in the lab frame, in the same units as
        ``angles_to_phi_vector`` (Å⁻¹ without the 2π factor that the
        B matrix omits).

    Raises
    ------
    ValueError
        If ``geometry.wavelength`` is None.
    ValueError
        If ``geometry.surface_normal`` and ``geometry.azimuth``
        are both None (no surface normal available).
    ValueError
        If ``geometry.sample.UB`` is None.
    """
    if geometry.wavelength is None:
        raise ValueError("surface calculations require geometry.wavelength to be set.")

    n_hkl = geometry.surface_normal
    if n_hkl is None:
        n_hkl = geometry.azimuth
    if n_hkl is None:
        raise ValueError(
            "surface calculations require geometry.surface_normal "
            "(or geometry.azimuth as a fallback) to be set."
        )

    if geometry.sample.UB is None:
        raise ValueError(
            "surface calculations require a UB matrix on the active sample. "
            "Set one with ub_identity(), ub_from_one_reflection(), or "
            "assign sample.UB directly."
        )

    # Resolve current angles
    if angles is None:
        angles = {s.name: s.angle for s in geometry._stages.values()}  # noqa: SLF001

    # Validate stage names up-front (raises KeyError for unknown names)
    for name in angles:
        geometry.stage(name)  # raises KeyError if not found

    # Compute rotation matrices statelessly (no save/restore needed).
    # Composition follows the BL1967 standard convention (Busing & Levy
    # 1967): outermost (floor-most) stage leftmost in the product, so
    # Z = R_0 @ R_1 @ ... @ R_{n-1} and likewise for D.  Z maps phi-frame
    # vectors to lab-frame vectors (v_lab = Z @ v_phi).
    from .rotation import _rotation_matrix_normalized

    Z = np.eye(3)
    for s in geometry.sample_stages:
        angle = angles.get(s.name, s.angle)
        Z = Z @ _rotation_matrix_normalized(s._axis_hat, angle)  # noqa: SLF001

    D = np.eye(3)
    for s in geometry.detector_stages:
        angle = angles.get(s.name, s.angle)
        D = D @ _rotation_matrix_normalized(s._axis_hat, angle)  # noqa: SLF001

    # Incident-beam direction: longitudinal basis vector (normalized)
    y_raw = np.asarray(geometry.basis["longitudinal"], dtype=float)
    y_hat = y_raw / np.linalg.norm(y_raw)

    # Diffracted-beam direction
    kf_hat = D @ y_hat
    kf_hat = kf_hat / np.linalg.norm(kf_hat)

    # Surface normal in phi frame then rotate to lab frame
    n_hkl_arr = np.asarray(n_hkl, dtype=float)
    n_phi = geometry.sample.UB @ n_hkl_arr
    n_phi_norm = np.linalg.norm(n_phi)
    if n_phi_norm < 1e-14:
        raise ValueError(
            "Surface normal UB @ n_hkl is the zero vector. "
            "Check the surface_normal Miller indices and the UB matrix."
        )
    n_phi = n_phi / n_phi_norm

    # Rotate to lab frame:  n_lab = Z @ n_phi
    n_lab = Z @ n_phi
    n_lab = n_lab / np.linalg.norm(n_lab)

    # Scattering vector in lab frame: Q_lab = (2π/λ)(kf - ki)
    two_pi_over_lambda = 2.0 * math.pi / geometry.wavelength
    Q_lab = two_pi_over_lambda * (kf_hat - y_hat)

    return y_hat, kf_hat, n_lab, Q_lab


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def incidence(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the angle of incidence αᵢ (degrees).

    αᵢ is the angle between the incoming beam and the sample surface
    (complement of the angle between the beam and the surface normal)::

        sin(αᵢ) = |ŷ · n̂_lab|

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``wavelength``, ``sample.UB``, and ``surface_normal``
        (or ``azimuth``) set.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, current stage angles are used.

    Returns
    -------
    float
        αᵢ in degrees, in [0°, 90°].

    Raises
    ------
    ValueError
        If any required attribute is not set.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.zaxis()
    >>> g.wavelength = 1.5406
    >>> g.surface_normal = (0, 0, 1)
    >>> ahd.ub_identity(g.sample)
    >>> ai = g.incidence({"alpha": 5.0, "Z": 0.0, "delta": 20.0, "gamma": 0.0})
    >>> round(ai, 4)
    5.0

    References
    ----------
    Lohmeier & Vlieg (1993), §4.2, eq. 15.
    """
    y_hat, _kf_hat, n_hat, _Q_lab = _surface_vectors(geometry, angles)
    sin_ai = float(np.dot(y_hat, n_hat))
    # Clamp for numerical safety
    sin_ai = max(-1.0, min(1.0, sin_ai))
    return math.degrees(math.asin(abs(sin_ai)))


def emergence(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the angle of emergence αf (degrees).

    αf is the angle between the diffracted beam and the sample surface::

        sin(αf) = |k̂f · n̂_lab|

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``wavelength``, ``sample.UB``, and ``surface_normal``
        (or ``azimuth``) set.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, current stage angles are used.

    Returns
    -------
    float
        αf in degrees, in [0°, 90°].

    Raises
    ------
    ValueError
        If any required attribute is not set.

    References
    ----------
    Lohmeier & Vlieg (1993), §4.2, eq. 16.
    """
    _y_hat, kf_hat, n_hat, _Q_lab = _surface_vectors(geometry, angles)
    sin_af = float(np.dot(kf_hat, n_hat))
    sin_af = max(-1.0, min(1.0, sin_af))
    return math.degrees(math.asin(abs(sin_af)))


# REMOVE-AT-V1.0: deprecated aliases for the canonical incidence /
# emergence functions.  Emit DeprecationWarning on every call and
# delegate to the canonical implementation.
def alpha_i(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """Deprecated alias for :func:`incidence`.  REMOVE-AT-V1.0."""
    import warnings

    warnings.warn(
        "ad_hoc_diffractometer.surface.alpha_i() is deprecated; "
        "use ad_hoc_diffractometer.surface.incidence() instead.  "
        "(issue #299)",
        DeprecationWarning,
        stacklevel=2,
    )
    return incidence(geometry, angles)


def alpha_f(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """Deprecated alias for :func:`emergence`.  REMOVE-AT-V1.0."""
    import warnings

    warnings.warn(
        "ad_hoc_diffractometer.surface.alpha_f() is deprecated; "
        "use ad_hoc_diffractometer.surface.emergence() instead.  "
        "(issue #299)",
        DeprecationWarning,
        stacklevel=2,
    )
    return emergence(geometry, angles)


def q_components(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> dict[str, float]:
    r"""
    Decompose Q into components perpendicular and parallel to the sample surface.

    ``Q_perp`` is the component of Q along the surface normal (out-of-plane).
    ``Q_par`` is the magnitude of the component of Q in the surface plane
    (in-plane).

    The signed perpendicular component ``Q_perp_signed`` is positive when Q
    has a component in the same direction as the outward surface normal.

    Units
    -----
    The returned values are in the same units as ``angles_to_phi_vector``,
    i.e. Å⁻¹ using the ``(2π/λ)(kf − ki)`` definition of Q (with full 2π).

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``wavelength``, ``sample.UB``, and ``surface_normal``
        (or ``azimuth``) set.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, current stage angles are used.

    Returns
    -------
    dict with keys:
        ``"Q_perp"`` : float
            Magnitude of the out-of-plane component (Å⁻¹), always ≥ 0.
        ``"Q_par"`` : float
            Magnitude of the in-plane component (Å⁻¹), always ≥ 0.
        ``"Q_perp_signed"`` : float
            Signed out-of-plane component (Å⁻¹); positive when Q points
            along the outward surface normal.
        ``"Q_total"`` : float
            Total \|Q\| (Å⁻¹).

    Raises
    ------
    ValueError
        If any required attribute is not set.

    References
    ----------
    Lohmeier & Vlieg (1993), §4.2, eq. 17.
    """
    _y_hat, _kf_hat, n_hat, Q_lab = _surface_vectors(geometry, angles)

    Q_perp_signed = float(np.dot(Q_lab, n_hat))
    Q_perp_vec = Q_perp_signed * n_hat
    Q_par_vec = Q_lab - Q_perp_vec
    Q_par = float(np.linalg.norm(Q_par_vec))
    Q_total = float(np.linalg.norm(Q_lab))

    return {
        "Q_perp": abs(Q_perp_signed),
        "Q_par": Q_par,
        "Q_perp_signed": Q_perp_signed,
        "Q_total": Q_total,
    }


def is_specular(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
    atol: float = 0.01,
) -> bool:
    """
    Return True when αᵢ ≈ αf within ``atol`` degrees.

    The specular condition (αᵢ = αf) is the constraint used in
    reflectometry and specular surface diffraction.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    angles : dict[str, float] or None
        Motor angles in degrees.
    atol : float
        Tolerance in degrees.  Default is 0.01°.

    Returns
    -------
    bool

    Raises
    ------
    ValueError
        If any required attribute is not set.
    """
    ai = incidence(geometry, angles)
    af = emergence(geometry, angles)
    return abs(ai - af) <= atol


def is_evanescent(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
    critical_angle_deg: float | None = None,
) -> bool:
    """
    Return True when αᵢ is below the critical angle (evanescent / surface-sensitive).

    If ``critical_angle_deg`` is ``None`` the function raises ``ValueError``
    since the critical angle depends on material properties (electron density)
    that are not stored on the geometry.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    angles : dict[str, float] or None
        Motor angles in degrees.
    critical_angle_deg : float
        Critical angle for total external reflection in degrees.  This is
        a material property (typically 0.1°–0.5° for hard X-rays) that
        must be supplied by the caller.

    Returns
    -------
    bool
        True if αᵢ < critical_angle_deg.

    Raises
    ------
    ValueError
        If ``critical_angle_deg`` is None.

    Examples
    --------
    >>> g.is_evanescent(critical_angle_deg=0.2)
    True   # if αᵢ < 0.2°
    """
    if critical_angle_deg is None:
        raise ValueError(
            "is_evanescent() requires critical_angle_deg to be supplied. "
            "The critical angle depends on the material's electron density "
            "and is not stored on the geometry.  "
            "Typical values: 0.1°–0.5° for hard X-rays."
        )
    ai = incidence(geometry, angles)
    return ai < critical_angle_deg
