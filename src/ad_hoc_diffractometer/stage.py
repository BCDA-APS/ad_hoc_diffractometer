"""
stage.py — rotary stage description.

A Stage represents one rotary axis of a diffractometer.  Its rotation axis
is stored internally as a signed numpy array (the internal representation).
The caller-facing +x/-z string notation is handled in axes.py.
"""

import numpy as np

from .rotation import rotation_matrix


class Stage:
    """
    One rotary stage of a diffractometer.

    Each stage is characterised by:
      - a name (e.g. 'mu', 'eta', 'chi')
      - a rotation axis expressed as a signed Cartesian vector in the lab
        frame (e.g. +XHAT for right-handed vertical, -ZHAT for left-handed
        lateral).  Use axes.parse_axis() to construct this from the
        caller-facing notation.
      - a parent stage name (the stage on which this one sits), or None if
        it sits on the floor / fixed lab frame.

    The sign of the axis vector encodes handedness:
        +nHat  =>  right-handed rotation about nHat
        -nHat  =>  left-handed rotation about nHat

    Parameters
    ----------
    name : str
        Human-readable name for this stage (e.g. 'mu', 'S2-1').
    axis : numpy.ndarray, shape (3,)
        Signed rotation axis vector in the lab frame (internal representation).
        Construct from caller-facing notation using axes.parse_axis().
    parent : str or None
        Name of the stage on which this stage is mounted, or None if it
        sits directly on the lab frame (floor).
    role : str, optional
        'sample' or 'detector', for bookkeeping.  Default is 'sample'.
    angle : float, optional
        Current angle setting in degrees.  Default is 0.0.
    """

    def __init__(
        self,
        name: str,
        axis: np.ndarray,
        parent: str | None = None,
        role: str = "sample",
        angle: float = 0.0,
    ):
        self.name = name
        self.axis = np.asarray(axis, dtype=float)
        self.parent = parent
        self.role = role
        self.angle = angle

    def rotation_matrix(self) -> np.ndarray:
        """Return the 3x3 rotation matrix for the current angle setting."""
        return rotation_matrix(self.axis, self.angle)

    def __repr__(self) -> str:
        return (
            f"Stage(name={self.name!r}, axis={self.axis}, "
            f"parent={self.parent!r}, role={self.role!r}, angle={self.angle})"
        )
