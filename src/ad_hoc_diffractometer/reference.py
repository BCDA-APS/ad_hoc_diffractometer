# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
reference.py — Reference pseudo-angle computations for diffraction modes.

Provides standalone functions for computing the physical pseudo-angles that
appear in :class:`~mode.ReferenceConstraint` conditions: incidence angle,
exit angle, azimuthal angle ψ, and lab-frame azimuthal angle naz.

These functions require the geometry's :attr:`surface_normal` or
:attr:`azimuthal_reference` to be set before calling.

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

Notes
-----
All functions accept a ``geometry`` instance and an optional ``angles`` dict
of motor angles (defaulting to the geometry's current stage angles when
``None``).  They raise :class:`ValueError` when the required reference vector
is not set on the geometry.

References
----------
* You, *J. Appl. Cryst.* **32**, 614-623 (1999), eqs. 10-11, 23.
* Lohmeier & Vlieg, *J. Appl. Cryst.* **26**, 706-716 (1993), §4.2.
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
