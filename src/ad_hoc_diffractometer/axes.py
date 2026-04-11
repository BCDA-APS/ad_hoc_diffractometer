# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
axes.py — caller-facing axis notation and physical direction mapping.

Provides a human-readable string notation for signed rotation axes:

    "+x"  ->  +XHAT  (right-handed rotation about the x-axis)
    "-x"  ->  -XHAT  (left-handed rotation about the x-axis)
    "+y"  ->  +YHAT
    "-y"  ->  -YHAT
    "+z"  ->  +ZHAT
    "-z"  ->  -ZHAT

Physical direction names ("vertical", "longitudinal", "lateral") are also
accepted and resolved against a caller-supplied basis dict that maps each
name to one of the three signed axis vectors.

For diffractometer geometries with tilted axes (e.g. kappa), kappa_axis()
computes the axis vector from the kappa angle alpha and the basis dict.

Internal stage attributes always store the axis as a numpy array.
This module handles only the conversion between the caller-facing string
notation and the internal numpy representation.

References
----------
International Tables for Crystallography, Vol. C, Section 2.2.6 (2006).
    DOI: 10.1107/97809553602060000103
    Confirms kappa axis tilted 50° from the omega (vertical) axis.
"""

import logging

import numpy as np

from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Canonical signed-axis string notation
# ---------------------------------------------------------------------------

# Maps caller-facing label -> numpy array (internal representation)
_LABEL_TO_VECTOR: dict[str, np.ndarray] = {
    "+x": +XHAT,
    "-x": -XHAT,
    "+y": +YHAT,
    "-y": -YHAT,
    "+z": +ZHAT,
    "-z": -ZHAT,
}

# Maps internal numpy array (as a tuple key) -> canonical label
_VECTOR_TO_LABEL: dict[tuple, str] = {tuple(v): k for k, v in _LABEL_TO_VECTOR.items()}


def parse_axis(label: str, basis: dict | None = None) -> np.ndarray:
    """
    Convert a caller-facing axis label to an internal numpy array.

    Accepted forms:

    1. Signed Cartesian label:  "+x", "-x", "+y", "-y", "+z", "-z"
       These are resolved directly against the standard basis vectors
       XHAT, YHAT, ZHAT regardless of the supplied basis dict.

    2. Physical direction name with optional sign prefix:
       "vertical", "+vertical", "-vertical",
       "longitudinal", "+longitudinal", "-longitudinal",
       "lateral", "+lateral", "-lateral"
       These are resolved against the supplied basis dict.

    Parameters
    ----------
    label : str
        Axis label in one of the forms above.  Case-insensitive.
    basis : dict or None
        Mapping from physical direction names to numpy arrays.
        Required if label contains a physical direction name.
        Ignored if label is a signed Cartesian label.

    Returns
    -------
    axis : numpy.ndarray, shape (3,)
        Signed axis vector for internal use.

    Raises
    ------
    ValueError
        If the label is not recognised, or if a physical direction name is
        given but no basis dict is supplied, or if the direction name is
        not in the basis dict.

    Examples
    --------
    >>> parse_axis("+x")
    array([1., 0., 0.])
    >>> parse_axis("-z")
    array([ 0.,  0., -1.])
    >>> basis = {"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT}
    >>> parse_axis("vertical", basis)
    array([1., 0., 0.])
    >>> parse_axis("-lateral", basis)
    array([ 0.,  0., -1.])
    """
    label = label.strip().lower()

    # 1. Signed Cartesian label
    if label in _LABEL_TO_VECTOR:
        return _LABEL_TO_VECTOR[label].copy()

    # 2. Physical direction name (with optional sign prefix)
    sign = +1.0
    name = label
    if label.startswith(("+", "-")):
        sign = +1.0 if label[0] == "+" else -1.0
        name = label[1:]

    if basis is None:
        raise ValueError(
            f"Axis label {label!r} is not a signed Cartesian label "
            f"(+x, -x, +y, -y, +z, -z); a basis dict is required to "
            f"resolve physical direction names."
        )

    if name not in basis:
        raise ValueError(
            f"Physical direction {name!r} not found in basis dict. "
            f"Available directions: {list(basis.keys())}."
        )

    return sign * np.asarray(basis[name], dtype=float)


def axis_label(vector: np.ndarray) -> str:
    """
    Convert an internal axis numpy array to its caller-facing string label.

    Only the six standard signed basis vectors are recognised.  For other
    vectors (e.g. a kappa-style tilted axis) the array is formatted
    numerically.

    Parameters
    ----------
    vector : numpy.ndarray, shape (3,)
        Signed axis vector.

    Returns
    -------
    label : str
        Caller-facing label, e.g. "+x", "-z", or "[0.5, 0.5, 0.707]".

    Examples
    --------
    >>> axis_label(XHAT)
    '+x'
    >>> axis_label(-ZHAT)
    '-z'
    """
    key = tuple(np.asarray(vector, dtype=float))
    if key in _VECTOR_TO_LABEL:
        return _VECTOR_TO_LABEL[key]
    return f"[{vector[0]:.4g}, {vector[1]:.4g}, {vector[2]:.4g}]"


def axis_from_physical(direction: str, sign: str, basis: dict) -> np.ndarray:
    """
    Convenience wrapper: resolve a physical direction and sign to an axis vector.

    Parameters
    ----------
    direction : str
        Physical direction name, e.g. "vertical", "lateral", "longitudinal".
    sign : str
        "+" for right-handed, "-" for left-handed rotation.
    basis : dict
        Mapping from physical direction names to numpy arrays.

    Returns
    -------
    axis : numpy.ndarray, shape (3,)

    Examples
    --------
    >>> basis = {"vertical": XHAT, "longitudinal": YHAT, "lateral": ZHAT}
    >>> axis_from_physical("lateral", "-", basis)
    array([ 0.,  0., -1.])
    """
    return parse_axis(f"{sign}{direction}", basis=basis)


def kappa_axis(alpha_deg: float, basis: dict | None = None) -> np.ndarray:
    """
    Compute the kappa rotation axis vector for a kappa-geometry diffractometer.

    The kappa axis lies in the vertical-lateral plane, tilted at angle alpha
    from the vertical axis toward the lateral axis:


        kappa_axis = vertical * cos(alpha) + lateral * sin(alpha)

    In the You (1999) / default basis (xHat=vertical, zHat=lateral):

        kappa_axis = XHAT * cos(alpha) + ZHAT * sin(alpha)

    At alpha=0  the axis is purely vertical (degenerate with komega).
    At alpha=90 the axis is purely lateral.
    Typical value: alpha=50 degrees (Walko 2016; ITC Vol. C Sec. 2.2.6;
    originally Enraf-Nonius).

    The axis is a unit vector by construction (cos²α + sin²α = 1).

    Parameters
    ----------
    alpha_deg : float
        Kappa tilt angle in degrees, measured from the vertical axis toward
        the lateral axis in the vertical-lateral plane.  Must be in (0, 90).
    basis : dict or None
        Mapping from physical direction names to numpy arrays.  If None,
        the default You (1999) basis is used:
            vertical  -> XHAT
            lateral   -> ZHAT

    Returns
    -------
    axis : numpy.ndarray, shape (3,)
        Unit vector of the kappa rotation axis in the lab frame.

    Raises
    ------
    ValueError
        If alpha_deg is not in the open interval (0, 90) degrees, or if
        the basis dict is missing 'vertical' or 'lateral' keys.

    Examples
    --------
    >>> kappa_axis(50.0)                    # typical kappa, default basis
    >>> kappa_axis(50.0, basis=my_basis)    # custom basis
    """
    if not (0.0 < alpha_deg < 90.0):
        raise ValueError(
            f"kappa alpha must be in (0, 90) degrees; got {alpha_deg}. "
            f"At 0 the kappa axis is degenerate with komega (vertical); "
            f"at 90 it is purely lateral."
        )

    if basis is None:
        vertical = XHAT
        lateral = ZHAT
    else:
        for key in ("vertical", "lateral"):
            if key not in basis:
                raise ValueError(
                    f"basis dict must contain 'vertical' and 'lateral' keys; "
                    f"missing: {key!r}.  Available: {list(basis.keys())}."
                )
        vertical = np.asarray(basis["vertical"], dtype=float)
        lateral = np.asarray(basis["lateral"], dtype=float)

    alpha_r = np.deg2rad(alpha_deg)
    return np.cos(alpha_r) * vertical + np.sin(alpha_r) * lateral
