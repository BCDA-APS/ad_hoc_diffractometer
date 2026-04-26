# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
stage.py — rotary stage description.

A Stage represents one rotary axis of a diffractometer.  Its rotation axis
is stored internally as a signed numpy array (the internal representation).
The caller-facing +x/-z string notation is handled in axes.py.
"""

import logging

import numpy as np

from .rotation import rotation_matrix

logger = logging.getLogger(__name__)


class Stage:
    """
    One rotary stage of a diffractometer.

    Each stage is characterized by:
      - a name (e.g. 'mu', 'eta', 'chi')
      - a rotation axis expressed as a signed Cartesian vector in the lab
        frame (e.g. +XHAT for right-handed vertical, -ZHAT for left-handed
        transverse).  Use axes.parse_axis() to construct this from the
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
        An arbitrary string label for bookkeeping.  The two conventional
        values are ``"sample"`` and ``"detector"``; any other string is
        accepted and can be used to model additional components such as
        analyzers, polarizers, slits, or azimuthal spinners.
        Geometry methods ``sample_stages`` and ``detector_stages`` filter
        by these conventional values; use
        :meth:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.stages_by_role`
        to query stages with any other role.  Default is ``"sample"``.
    angle : float, optional
        Current angle setting in degrees.  Default is 0.0.
    limits : tuple of (float, float), optional
        (min_angle, max_angle) in degrees.  Default is (-180.0, 180.0).
        Must satisfy min_angle < max_angle.

    Raises
    ------
    ValueError
        If limits[0] >= limits[1].
    """

    def __init__(
        self,
        name: str,
        axis: np.ndarray,
        parent: str | None = None,
        role: str = "sample",
        angle: float = 0.0,
        limits: tuple[float, float] = (-180.0, 180.0),
    ):
        self.name = name
        self.axis = np.asarray(axis, dtype=float)
        # Pre-normalized axis for fast rotation computation (avoids repeated
        # normalization in the inner loop of the forward solver).
        ax = self.axis
        ax_norm = np.linalg.norm(ax)
        self._axis_hat: np.ndarray = ax / ax_norm if ax_norm > 0 else ax.copy()
        self.parent = parent
        self.role = role
        self.angle = angle
        self.limits = limits  # validated via property setter

    @property
    def limits(self) -> tuple[float, float]:
        """Motor angle limits (min_angle, max_angle) in degrees."""
        return self._limits

    @limits.setter
    def limits(self, value: tuple[float, float]) -> None:
        lo, hi = float(value[0]), float(value[1])
        if lo >= hi:
            raise ValueError(
                f"Stage {self.name!r}: limits min ({lo}) must be less than max ({hi})."
            )
        self._limits = (lo, hi)

    def in_limits(self, angle_deg: float) -> bool:
        """
        Return True if angle_deg is within [min_angle, max_angle] (inclusive).

        Parameters
        ----------
        angle_deg : float
            Angle to test, in degrees.

        Returns
        -------
        bool
        """
        lo, hi = self._limits
        return lo <= angle_deg <= hi

    def rotation_matrix(self) -> np.ndarray:
        """Return the 3x3 rotation matrix for the current angle setting."""
        return rotation_matrix(self.axis, self.angle)

    def __repr__(self) -> str:
        return (
            f"Stage(name={self.name!r}, axis={self.axis}, "
            f"parent={self.parent!r}, role={self.role!r}, angle={self.angle}, "
            f"limits={self.limits})"
        )

    def to_dict(self) -> dict:
        """
        Return a JSON-serialisable ``dict`` representing this stage.

        Returns
        -------
        dict
            Keys: ``"name"`` (str), ``"axis"`` (list of 3 float),
            ``"role"`` (str), ``"parent"`` (str or None),
            ``"angle"`` (float), ``"limits"`` (list of 2 float).
        """
        return {
            "name": self.name,
            "axis": [float(x) for x in self.axis],
            "role": self.role,
            "parent": self.parent,
            "angle": float(self.angle),
            "limits": [float(self.limits[0]), float(self.limits[1])],
        }

    @classmethod
    def from_dict(cls, d: dict) -> "Stage":
        """
        Reconstruct a :class:`Stage` from a dict produced by :meth:`to_dict`.

        Parameters
        ----------
        d : dict
            Must contain ``"name"``, ``"axis"``, ``"role"``.  ``"parent"``,
            ``"angle"``, and ``"limits"`` are optional (defaulting to
            ``None``, ``0.0``, and ``(-180.0, 180.0)`` respectively).

        Returns
        -------
        Stage
        """
        import numpy as np

        s = cls(
            name=d["name"],
            axis=np.array(d["axis"], dtype=float),
            role=d["role"],
            parent=d.get("parent"),
            limits=tuple(d.get("limits", (-180.0, 180.0))),
        )
        s.angle = float(d.get("angle", 0.0))
        return s
