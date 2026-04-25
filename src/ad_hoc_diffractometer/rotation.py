# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
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
        +nHat  =>  right-handed rotation about nHat
        -nHat  =>  left-handed rotation about nHat
                   (equivalent to right-handed with negated angle)

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
