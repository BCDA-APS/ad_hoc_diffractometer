# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
rotation.py — rotation matrix calculation.

Provides the Rodrigues rotation formula for computing 3x3 rotation matrices
from a unit axis vector and an angle.
"""

import logging

import numpy as np

logger = logging.getLogger(__name__)

logger = logging.getLogger(__name__)


def rotation_matrix(axis: np.ndarray, angle_deg: float) -> np.ndarray:
    """
    Compute a right-handed rotation matrix about a unit axis vector.

    Uses the Rodrigues formula:

        R = I cos(theta) + (1 - cos(theta))(n ⊗ n) + sin(theta)[n]_x

    where [n]_x is the skew-symmetric cross-product matrix of n.

    For a left-handed rotation (e.g. eta, phi, delta in You 1999), pass the
    negated axis: rotation_matrix(-ZHAT, angle_deg).

    The sign of the axis vector encodes handedness:
    ``+n_axis`` is a right-handed rotation about ``n_axis``,
    ``-n_axis`` is left-handed (equivalent to right-handed with
    a negated angle).

    .. note::

       Here ``n_axis`` is the **per-stage rotation axis** — a property
       of how each stage is wired in its YAML definition.  Do not
       confuse it with ``n̂`` (``n_hat`` in mode ``extras``,
       :attr:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.surface_normal`,
       :attr:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.azimuth`),
       which is the user-facing **reference vector** consumed by
       :class:`~ad_hoc_diffractometer.mode.ReferenceConstraint` modes.
       See the :ref:`glossary` entries for "n̂ (reference vector)" and
       "Stage rotation axis" for the disambiguation.

    Parameters
    ----------
    axis : numpy.ndarray, shape (3,)
        Rotation axis vector.  Need not be normalized; it will be normalized
        internally.
    angle_deg : float
        Rotation angle in degrees (right-handed sense about the given axis).

    Returns
    -------
    R : numpy.ndarray, shape (3, 3)
        Orthogonal rotation matrix with det(R) = +1.

    Examples
    --------
    >>> import numpy as np
    >>> from ad_hoc_diffractometer.constants import XHAT, ZHAT
    >>> rotation_matrix(XHAT, 0)          # identity
    >>> rotation_matrix(-ZHAT, 30)        # left-handed 30 deg about z
    """
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    return _rotation_matrix_normalized(n, angle_deg)


def _rotation_matrix_normalized(n: np.ndarray, angle_deg: float) -> np.ndarray:
    """
    Compute a rotation matrix from a pre-normalized axis vector.

    Fast path: skips the normalization step.  The caller must guarantee
    that ``n`` is a unit vector (``|n| = 1``).  Produces bit-identical
    results to :func:`rotation_matrix` when the input axis is already
    normalized.

    Parameters
    ----------
    n : numpy.ndarray, shape (3,)
        **Unit** rotation axis vector.
    angle_deg : float
        Rotation angle in degrees.

    Returns
    -------
    R : numpy.ndarray, shape (3, 3)
    """
    theta = np.deg2rad(angle_deg)
    c = np.cos(theta)
    s = np.sin(theta)
    nx, ny, nz = n
    skew = np.array(
        [
            [0, -nz, ny],
            [nz, 0, -nx],
            [-ny, nx, 0],
        ]
    )
    return c * np.eye(3) + (1 - c) * np.outer(n, n) + s * skew


def _rotation_matrix_and_derivative_normalized(
    n: np.ndarray, angle_deg: float
) -> tuple[np.ndarray, np.ndarray]:
    """
    Compute both the rotation matrix and its derivative with respect to angle.

    Returns ``(R, dR/dtheta)`` where ``theta`` is in **radians** — the
    derivative is with respect to the radian-valued angle, not degrees.
    The caller is responsible for the deg-to-rad conversion factor
    (``pi/180``) when composing chain-rule Jacobians that take degree
    inputs.

    Uses the same Rodrigues decomposition as
    :func:`_rotation_matrix_normalized` and shares the trig values, so the
    marginal cost of the derivative is negligible.

    The derivative formula:

        dR/dtheta = -sin(theta) * (I - n (x) n) + cos(theta) * [n]_x

    Parameters
    ----------
    n : numpy.ndarray, shape (3,)
        **Unit** rotation axis vector.
    angle_deg : float
        Rotation angle in degrees.

    Returns
    -------
    R : numpy.ndarray, shape (3, 3)
        Rotation matrix (identical to ``_rotation_matrix_normalized``).
    dR_dtheta : numpy.ndarray, shape (3, 3)
        Derivative of R with respect to theta (radians).
    """
    theta = np.deg2rad(angle_deg)
    c = np.cos(theta)
    s = np.sin(theta)
    nx, ny, nz = n
    skew = np.array(
        [
            [0, -nz, ny],
            [nz, 0, -nx],
            [-ny, nx, 0],
        ]
    )
    nn = np.outer(n, n)
    R = c * np.eye(3) + (1 - c) * nn + s * skew
    # dR/dtheta = -sin(theta)(I - n⊗n) + cos(theta)[n]×
    dR = -s * (np.eye(3) - nn) + c * skew
    return R, dR
