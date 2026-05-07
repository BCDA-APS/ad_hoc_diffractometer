# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
reference.py — Reference pseudo-angle computations for diffraction modes.

Provides standalone functions for computing the physical pseudo-angles that
appear in :class:`~mode.ReferenceConstraint` conditions: incidence angle,
exit angle, azimuthal angle ψ, lab-frame azimuthal angle naz, and the
SPEC ``OMEGA`` pseudo-angle (angle between Q and the chi-circle plane).

These functions require the geometry's :attr:`surface_normal` or
:attr:`azimuthal_reference` to be set before calling, **except** for
:func:`omega_pseudo`, which is a pure motor-frame quantity and needs
neither.

Functions
---------
:func:`incidence_angle`
    Angle of incidence α_i between the incident beam and the sample surface.

:func:`exit_angle`
    Angle of exit β_out between the diffracted beam and the sample surface.

:func:`psi_angle`
    Azimuthal angle ψ of the reference vector n̂ about Q (You 1999, eq. 23).

:func:`naz_angle`
    Azimuthal angle of n̂ projected onto the lab-frame horizontal plane.

:func:`omega_pseudo`
    SPEC ``OMEGA`` pseudo-angle — angle between Q and the plane of the
    chi circle (SPEC ``psic`` ``def OMEGA 'Q[6]'``; You 1999 §5).

Notes
-----
All functions accept a ``geometry`` instance and an optional ``angles`` dict
of motor angles (defaulting to the geometry's current stage angles when
``None``).  They raise :class:`ValueError` when the required reference vector
is not set on the geometry (where applicable).

References
----------
* You, *J. Appl. Cryst.* **32**, 614-623 (1999), eqs. 10-11, 23, §5.
* Lohmeier & Vlieg, *J. Appl. Cryst.* **26**, 706-716 (1993), §4.2.
* Certified Scientific Software, *SPEC psic help*,
  https://certif.com/spec_help/psic.html — ``def OMEGA 'Q[6]'``.
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

logger = logging.getLogger(__name__)


def _require_surface_normal(geometry: AdHocDiffractometer) -> None:
    """Raise ValueError if surface_normal is not set on the geometry."""
    sn = geometry.surface_normal
    if sn is None:
        raise ValueError(
            f"geometry '{geometry.name}': surface_normal must be set before "
            "calling reference pseudo-angle functions.  "
            "Set g.surface_normal = (h, k, l) with the surface normal "
            "expressed as Miller indices."
        )


def _require_azimuthal_reference(geometry: AdHocDiffractometer) -> None:
    """Raise ValueError if azimuthal_reference is not set on the geometry."""
    ar = geometry.azimuthal_reference
    if ar is None:
        raise ValueError(
            f"geometry '{geometry.name}': azimuthal_reference must be set before "
            "computing the ψ angle.  "
            "Set g.azimuthal_reference = (h, k, l) with the reference direction "
            "expressed as Miller indices."
        )


def incidence_angle(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the angle of incidence α_i in degrees.

    The angle between the incident beam and the sample surface (the
    complement of the angle between the incident beam and the surface
    normal).  Positive when the beam strikes the front face.

    Requires :attr:`~geometry.AdHocDiffractometer.surface_normal` to be set.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer instance.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, the geometry's current
        stage angles are used.

    Returns
    -------
    float
        Incidence angle α_i in degrees.

    Raises
    ------
    ValueError
        If ``geometry.surface_normal`` is ``None``.

    References
    ----------
    * You (1999), eq. 10.
    * Lohmeier & Vlieg (1993), §4.2.
    """
    _require_surface_normal(geometry)
    return geometry.alpha_i(angles=angles)


def exit_angle(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the angle of exit β_out in degrees.

    The angle between the diffracted beam and the sample surface.
    Positive when the diffracted beam exits through the front face.

    Requires :attr:`~geometry.AdHocDiffractometer.surface_normal` to be set.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer instance.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, the geometry's current
        stage angles are used.

    Returns
    -------
    float
        Exit angle β_out in degrees.

    Raises
    ------
    ValueError
        If ``geometry.surface_normal`` is ``None``.

    References
    ----------
    * You (1999), eq. 11.
    * Lohmeier & Vlieg (1993), §4.2.
    """
    _require_surface_normal(geometry)
    return geometry.alpha_f(angles=angles)


def psi_angle(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the azimuthal angle ψ of the reference vector n̂ about Q.

    ψ is the angle between the azimuthal reference direction (projected
    onto the plane perpendicular to Q) and the incident beam direction
    (also projected onto that plane).  ψ = 0 when the reference vector
    lies in the scattering plane on the same side as the incident beam.

    Requires :attr:`~geometry.AdHocDiffractometer.azimuthal_reference`
    to be set.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer instance.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, the geometry's current
        stage angles are used.

    Returns
    -------
    float
        Azimuthal angle ψ in degrees, in the range (−180°, +180°].

    Raises
    ------
    ValueError
        If ``geometry.azimuthal_reference`` is ``None``.
    ValueError
        If the reference vector is parallel to Q (ψ is undefined).

    References
    ----------
    * You (1999), eq. 23.
    """
    _require_azimuthal_reference(geometry)
    return geometry.psi(angles=angles)


def naz_angle(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the lab-frame azimuthal angle of n̂ (naz) in degrees.

    naz is the azimuthal angle of the surface normal n̂ projected onto
    the horizontal plane of the lab frame (the plane perpendicular to
    the vertical axis).  It describes the in-plane orientation of the
    sample surface relative to the lab coordinate system.

    Requires :attr:`~geometry.AdHocDiffractometer.surface_normal` to be set.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer instance.
    angles : dict[str, float] or None
        Motor angles in degrees.  If ``None``, the geometry's current
        stage angles are used.

    Returns
    -------
    float
        Azimuthal angle naz in degrees, in the range (−180°, +180°].

    Raises
    ------
    ValueError
        If ``geometry.surface_normal`` is ``None``.

    References
    ----------
    * You (1999).
    """
    _require_surface_normal(geometry)

    if angles is None:
        angles = {s.name: s.angle for s in geometry._stages.values()}  # noqa: SLF001

    # Surface normal in phi frame
    n_hkl = np.asarray(geometry.surface_normal, dtype=float)
    n_phi = geometry.sample.UB @ n_hkl
    n_mag = float(np.linalg.norm(n_phi))
    if n_mag < 1e-14:  # pragma: no cover
        raise ValueError(
            "naz_angle: surface normal maps to zero in the phi frame. "
            "Check that the UB matrix is non-singular."
        )
    n_phi_hat = n_phi / n_mag

    # Vertical axis in the lab frame
    vertical = np.asarray(
        geometry.basis.get("vertical", np.array([1.0, 0.0, 0.0])),
        dtype=float,
    )
    vertical_hat = vertical / np.linalg.norm(vertical)

    # Transverse axis = longitudinal × vertical (right-handed)
    longitudinal = np.asarray(
        geometry.basis.get("longitudinal", np.array([0.0, 1.0, 0.0])),
        dtype=float,
    )
    longitudinal_hat = longitudinal / np.linalg.norm(longitudinal)

    # Project n_phi onto the horizontal plane (perpendicular to vertical)
    n_horiz = n_phi_hat - np.dot(n_phi_hat, vertical_hat) * vertical_hat
    n_horiz_mag = float(np.linalg.norm(n_horiz))
    if n_horiz_mag < 1e-10:
        # n̂ is vertical — naz is undefined (return 0 by convention)
        return 0.0

    n_horiz_hat = n_horiz / n_horiz_mag

    # naz = angle of n_horiz from longitudinal axis, in the horizontal plane
    cos_naz = float(np.clip(np.dot(n_horiz_hat, longitudinal_hat), -1.0, 1.0))
    sin_naz = float(np.dot(vertical_hat, np.cross(longitudinal_hat, n_horiz_hat)))
    return math.degrees(math.atan2(sin_naz, cos_naz))


# ---------------------------------------------------------------------------
# SPEC OMEGA pseudo-angle (psic)
# ---------------------------------------------------------------------------


def _chi_stage_axis_lab(
    geometry: AdHocDiffractometer,
    angles: dict[str, float],
) -> np.ndarray:
    """
    Return the chi stage's rotation-axis unit vector expressed in the lab frame.

    The chi-circle plane is the plane perpendicular to this axis, so this
    vector is also the normal to the chi-circle plane.  All sample stages
    *outside* (i.e. lower in the stack than) chi rotate the chi axis;
    stages *inside* chi (phi, etc.) do not.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must contain a sample stage named ``"chi"``.
    angles : dict[str, float]
        Motor angles in degrees, keyed by stage name.

    Returns
    -------
    numpy.ndarray, shape (3,)
        Unit vector along the chi rotation axis in the lab frame.

    Raises
    ------
    KeyError
        If the geometry has no stage named ``"chi"``.
    """
    from .rotation import _rotation_matrix_normalized

    # Find the chi stage and the sample stages that come before it
    sample_stages = list(geometry.sample_stages)
    chi_index = None
    for i, s in enumerate(sample_stages):
        if s.name == "chi":
            chi_index = i
            break
    if chi_index is None:
        raise KeyError(
            f"omega_pseudo: geometry {geometry.name!r} has no sample stage "
            "named 'chi'.  The OMEGA pseudo-angle (SPEC Q[6]) is defined "
            "only for psic-family geometries with a chi circle."
        )

    chi_stage = sample_stages[chi_index]
    outer_stages = sample_stages[:chi_index]

    # Apply the rotations of all stages *before* chi to the chi axis
    R = np.eye(3)
    for s in outer_stages:
        angle = angles.get(s.name, s.angle)
        R = _rotation_matrix_normalized(s._axis_hat, angle) @ R  # noqa: SLF001

    chi_axis_lab = R @ np.asarray(chi_stage._axis_hat, dtype=float)  # noqa: SLF001
    n = float(np.linalg.norm(chi_axis_lab))
    if n < 1e-14:  # pragma: no cover
        raise ValueError(
            "omega_pseudo: chi stage axis vector has zero length after rotation."
        )
    return chi_axis_lab / n


def omega_pseudo(
    geometry: AdHocDiffractometer,
    angles: dict[str, float] | None = None,
) -> float:
    """
    Compute the SPEC ``OMEGA`` pseudo-angle in degrees.

    OMEGA is the angle between the scattering vector Q and the **plane of
    the chi circle**, taking values in ``[-90°, +90°]``.  When
    ``OMEGA = 0`` the scattering vector lies in the chi-circle plane and
    the diffractometer is in the bisecting condition; non-zero OMEGA
    means Q is tilted out of that plane.

    The chi-circle plane is the plane perpendicular to the chi rotation
    axis; that axis itself is rotated by every sample stage that sits
    *outside* (i.e. lower in the stack than) chi.  In psic the outer
    sample stages are ``mu`` and ``eta``.

    Algorithm
    ---------
    1. Compute the chi-axis unit vector in the lab frame, applying the
       rotations of all sample stages outside chi (e.g. ``mu`` and
       ``eta`` in psic).
    2. Compute the scattering vector ``Q_lab`` from the supplied motor
       angles.
    3. The OMEGA pseudo-angle is::

           sin(OMEGA) = (Q̂_lab · n̂_chi)

       where ``n̂_chi`` is the unit chi-axis vector in the lab frame.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Diffractometer geometry.  Must have ``wavelength`` set and must
        contain a sample stage named ``"chi"``.
    angles : dict[str, float] or None
        Motor angles in degrees, keyed by stage name.  If ``None``
        (default), the geometry's current stage angles are used.

    Returns
    -------
    float
        OMEGA pseudo-angle in degrees, in ``[-90°, +90°]``.

    Raises
    ------
    ValueError
        If ``geometry.wavelength`` is ``None``.
    KeyError
        If the geometry has no stage named ``"chi"``.

    Notes
    -----
    Unlike :func:`incidence_angle`, :func:`exit_angle`, :func:`psi_angle`,
    and :func:`naz_angle`, this function does **not** require any
    reference vector (``surface_normal`` / ``azimuthal_reference``) to be
    set on the geometry.  OMEGA is a pure motor-frame quantity defined
    by the diffractometer's internal geometry.

    The sign convention follows the right-hand rule about the chi axis:
    OMEGA is positive when Q has a positive component along the chi-axis
    direction defined by the geometry's stage configuration.

    OMEGA is **independent of the inner sample stages** (``chi`` and
    ``phi`` in psic): they rotate Q but they also rotate the
    chi-circle plane normal — the relative angle between them is
    preserved.  The pseudo-angle therefore depends only on the outer
    sample stages, the detector stages, and the wavelength.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.make_geometry("psic")
    >>> g.wavelength = 1.5406
    >>> # Bisecting position: mu = nu = 0, eta = delta/2 — Q in chi-plane
    >>> ahd.reference.omega_pseudo(
    ...     g,
    ...     {"mu": 0, "eta": 15.0, "chi": 90, "phi": 0,
    ...      "nu": 0, "delta": 30.0},
    ... )  # doctest: +SKIP
    0.0

    References
    ----------
    * Certified Scientific Software, *SPEC psic help*,
      https://certif.com/spec_help/psic.html —
      ``def OMEGA 'Q[6]'`` — "The angle between Q and the plane of
      the chi circle."
    * H. You, *J. Appl. Cryst.* **32**, 614-623 (1999), §5.
    """
    if geometry.wavelength is None:
        raise ValueError("omega_pseudo() requires geometry.wavelength to be set.")

    if angles is None:
        angles = {s.name: s.angle for s in geometry._stages.values()}  # noqa: SLF001

    # Lab-frame chi-circle axis (= normal to the chi-circle plane)
    n_chi_lab = _chi_stage_axis_lab(geometry, angles)

    # Lab-frame scattering vector
    from .rotation import _rotation_matrix_normalized

    D = np.eye(3)
    for s in geometry.detector_stages:
        angle = angles.get(s.name, s.angle)
        D = _rotation_matrix_normalized(s._axis_hat, angle) @ D  # noqa: SLF001

    y_raw = np.asarray(geometry.basis["longitudinal"], dtype=float)
    y_hat = y_raw / np.linalg.norm(y_raw)

    Q_lab = D @ y_hat - y_hat
    Q_mag = float(np.linalg.norm(Q_lab))
    if Q_mag < 1e-12:
        # No scattering — Q is undefined; OMEGA is undefined; return 0
        return 0.0
    Q_hat = Q_lab / Q_mag

    sin_omega = float(np.clip(np.dot(Q_hat, n_chi_lab), -1.0, 1.0))
    return math.degrees(math.asin(sin_omega))
