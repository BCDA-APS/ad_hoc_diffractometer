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

Internal stage attributes always store the axis as a numpy array.
This module handles only the conversion between the caller-facing string
notation and the internal numpy representation.
"""

import numpy as np

from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT

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
