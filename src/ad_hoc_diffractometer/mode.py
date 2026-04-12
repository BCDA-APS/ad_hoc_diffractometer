# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
mode.py — Diffraction mode / operating constraints.

A diffraction mode specifies which degrees of freedom are constrained
during a diffraction calculation and which are free.  Modes are a
first-class concern of the geometry description.

Classes
-------
DiffractionMode
    Abstract base class for all modes.  Declares which stages are
    constrained (frozen) and which remain free, plus optional cut-points
    for branch selection.

FixedAngleMode
    Constrains one named stage to a fixed angle value.

BisectingMode
    Bisecting geometry: omega = ttheta/2 (half the detector angle).
    For psic-family geometries this means eta = delta/2.

ModeDict
    An ordered, guarded dict of named modes.  Validates that all entries
    are DiffractionMode instances; prevents inserting non-Mode values.

Notes
-----
Cut-points control which branch of a multi-valued angle solution is
returned (SPEC #G4 convention).  A cut-point for stage X with value C
means the returned angle is chosen in the interval [C, C + 360°).

Frozen-angle values mirror the SPEC #G0 convention: when a mode holds a
stage fixed, that stage's frozen angle is recorded on the mode object as
``frozen_angles``.

References
----------
Busing & Levy, Acta Cryst. 22, 457-464 (1967).
You, J. Appl. Cryst. 32, 614-623 (1999).
Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), section 3.3.
SPEC #G0, #G4 lines — references/2020-12-13-fourcc-alignment-7-id-c/spec_G_lines.md
"""

from __future__ import annotations

import logging
from abc import ABC
from abc import abstractmethod
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    pass

logger = logging.getLogger(__name__)


class DiffractionMode(ABC):
    """
    Abstract base class for diffraction modes.

    A mode specifies:
    - which stages are frozen (held at a fixed angle during the calculation)
    - which stages are free (computed by the forward calculation)
    - per-stage cut-points for branch selection (SPEC #G4)

    Subclasses must implement :meth:`constrain` and :meth:`constrained_stages`.

    Parameters
    ----------
    frozen_angles : dict[str, float] or None
        Mapping of stage name -> frozen angle (degrees).  Stages in this
        dict are held at the given angle during diffraction calculations.
        ``None`` is equivalent to an empty dict (no frozen angles).
    cut_points : dict[str, float] or None
        Mapping of stage name -> cut-point angle (degrees).  The returned
        solution angle for that stage lies in ``[cut, cut + 360)``.
        ``None`` is equivalent to an empty dict (no cut-points set, i.e.
        the default branch [-180, 180) is used).

    """

    def __init__(
        self,
        frozen_angles: dict[str, float] | None = None,
        cut_points: dict[str, float] | None = None,
    ) -> None:
        self.frozen_angles: dict[str, float] = dict(frozen_angles or {})
        self.cut_points: dict[str, float] = dict(cut_points or {})

    @property
    @abstractmethod
    def constrained_stages(self) -> list[str]:
        """
        Return the list of stage names that are constrained (not free) in this mode.

        Constrained stages include both:

        - frozen stages (held at a fixed angle)
        - stages whose angle is determined by a relationship to another stage
          (e.g. bisecting: omega = ttheta/2)

        Returns
        -------
        list of str
        """

    def apply_cut_point(self, stage_name: str, angle_deg: float) -> float:
        """
        Apply the cut-point for *stage_name* to *angle_deg*.

        If no cut-point is set for this stage, the angle is returned
        unchanged.

        A cut-point C means the returned angle lies in ``[C, C + 360)``.

        Parameters
        ----------
        stage_name : str
            Stage to apply the cut-point to.
        angle_deg : float
            Angle to normalise (degrees).

        Returns
        -------
        float
            Angle normalised to the cut-point interval.
        """
        if stage_name not in self.cut_points:
            return angle_deg
        cut = self.cut_points[stage_name]
        # Shift into [cut, cut + 360)
        remainder = (angle_deg - cut) % 360.0
        return cut + remainder

    def __repr__(self) -> str:  # pragma: no cover
        # Each concrete subclass provides its own __repr__; this fallback is
        # kept for completeness but is never reached in production use.
        return (
            f"{type(self).__name__}("
            f"frozen_angles={self.frozen_angles!r}, "
            f"cut_points={self.cut_points!r})"
        )

    def __eq__(self, other: object) -> bool:  # pragma: no cover
        # Each concrete subclass overrides __eq__ with type-specific checks;
        # this base implementation is never reached in production use.
        if type(self) is not type(other):
            return False
        return (
            self.frozen_angles == other.frozen_angles  # type: ignore[attr-defined]
            and self.cut_points == other.cut_points  # type: ignore[attr-defined]
        )


class FixedAngleMode(DiffractionMode):
    """
    Mode that constrains one stage to a fixed angle.

    The stage named ``stage`` is held at ``value`` degrees for all
    diffraction calculations performed in this mode.

    Parameters
    ----------
    stage : str
        Name of the stage to freeze.
    value : float
        Angle (degrees) at which to freeze the stage.
    cut_points : dict[str, float] or None
        Branch-selection cut-points (see :class:`DiffractionMode`).

    Examples
    --------
    >>> mode = FixedAngleMode(stage="chi", value=90.0)
    >>> mode.frozen_angles
    {'chi': 90.0}
    >>> mode.constrained_stages
    ['chi']
    """

    def __init__(
        self,
        stage: str,
        value: float,
        cut_points: dict[str, float] | None = None,
    ) -> None:
        super().__init__(
            frozen_angles={stage: float(value)},
            cut_points=cut_points,
        )
        self._stage = stage
        self._value = float(value)

    @property
    def constrained_stages(self) -> list[str]:
        """Stages constrained by this mode: just the one frozen stage."""
        return [self._stage]

    def __repr__(self) -> str:
        return (
            f"FixedAngleMode(stage={self._stage!r}, value={self._value!r}, "
            f"cut_points={self.cut_points!r})"
        )

    def __eq__(self, other: object) -> bool:
        if type(self) is not type(other):
            return False
        return (
            self._stage == other._stage  # type: ignore[attr-defined]
            and self._value == other._value  # type: ignore[attr-defined]
            and self.cut_points == other.cut_points  # type: ignore[attr-defined]
        )


class BisectingMode(DiffractionMode):
    """
    Bisecting geometry: the sample rotation is set to half the detector angle.

    In the four-circle / psic convention, "bisecting" means that the
    incident and diffracted beam make equal angles with the diffracting
    planes.  Concretely:

    - For ``fourcv`` / ``fourch``: ``omega = ttheta / 2``
    - For ``psic``: ``eta = delta / 2`` (with ``mu = nu = 0`` by convention)

    This mode records the convention via ``sample_stage`` (the stage that is
    driven to half the detector angle) and ``detector_stage`` (the reference
    detector stage).

    Parameters
    ----------
    sample_stage : str
        Name of the sample rotation stage driven to half the detector angle
        (e.g. ``"omega"`` for fourcv, ``"eta"`` for psic).
    detector_stage : str
        Name of the detector stage whose angle determines the constraint
        (e.g. ``"ttheta"`` for fourcv, ``"delta"`` for psic).
    frozen_angles : dict[str, float] or None
        Additional stages held fixed (e.g. ``{"mu": 0.0, "nu": 0.0}`` for psic).
    cut_points : dict[str, float] or None
        Branch-selection cut-points (see :class:`DiffractionMode`).

    Examples
    --------
    >>> mode = BisectingMode(sample_stage="omega", detector_stage="ttheta")
    >>> mode.constrained_stages
    ['omega']
    >>> mode_psic = BisectingMode(
    ...     sample_stage="eta",
    ...     detector_stage="delta",
    ...     frozen_angles={"mu": 0.0, "nu": 0.0},
    ... )
    >>> mode_psic.constrained_stages
    ['eta', 'mu', 'nu']
    """

    def __init__(
        self,
        sample_stage: str,
        detector_stage: str,
        frozen_angles: dict[str, float] | None = None,
        cut_points: dict[str, float] | None = None,
    ) -> None:
        super().__init__(frozen_angles=frozen_angles, cut_points=cut_points)
        self._sample_stage = sample_stage
        self._detector_stage = detector_stage

    @property
    def sample_stage(self) -> str:
        """The sample rotation stage driven to half the detector angle."""
        return self._sample_stage

    @property
    def detector_stage(self) -> str:
        """The detector stage that determines the bisecting constraint."""
        return self._detector_stage

    @property
    def constrained_stages(self) -> list[str]:
        """
        Stages constrained by the bisecting condition.

        Always includes ``sample_stage``, plus any additional stages in
        ``frozen_angles``.
        """
        constrained = [self._sample_stage]
        for name in self.frozen_angles:
            if name not in constrained:
                constrained.append(name)
        return constrained

    def __repr__(self) -> str:
        return (
            f"BisectingMode(sample_stage={self._sample_stage!r}, "
            f"detector_stage={self._detector_stage!r}, "
            f"frozen_angles={self.frozen_angles!r}, "
            f"cut_points={self.cut_points!r})"
        )

    def __eq__(self, other: object) -> bool:
        if type(self) is not type(other):
            return False
        return (
            self._sample_stage == other._sample_stage  # type: ignore[attr-defined]
            and self._detector_stage == other._detector_stage  # type: ignore[attr-defined]
            and self.frozen_angles == other.frozen_angles  # type: ignore[attr-defined]
            and self.cut_points == other.cut_points  # type: ignore[attr-defined]
        )


class ModeDict:
    """
    An ordered, guarded dict of named diffraction modes.

    Only :class:`DiffractionMode` instances may be stored.  Keys are
    mode names (str).  Iteration follows insertion order.

    Parameters
    ----------
    modes : dict[str, DiffractionMode] or None
        Initial modes.  If ``None``, an empty dict is created.

    Raises
    ------
    TypeError
        If any value is not a :class:`DiffractionMode` instance.

    Examples
    --------
    >>> md = ModeDict({"bisecting": BisectingMode("omega", "ttheta")})
    >>> "bisecting" in md
    True
    >>> len(md)
    1
    """

    def __init__(self, modes: dict[str, DiffractionMode] | None = None) -> None:
        self._data: dict[str, DiffractionMode] = {}
        if modes:
            for name, mode in modes.items():
                self[name] = mode  # validates type

    def __setitem__(self, name: str, mode: DiffractionMode) -> None:
        if not isinstance(mode, DiffractionMode):
            raise TypeError(
                f"ModeDict values must be DiffractionMode instances; "
                f"got {type(mode).__name__!r} for key {name!r}."
            )
        self._data[name] = mode

    def __getitem__(self, name: str) -> DiffractionMode:
        return self._data[name]

    def __delitem__(self, name: str) -> None:
        del self._data[name]

    def __contains__(self, name: object) -> bool:
        return name in self._data

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self) -> str:
        return f"ModeDict({self._data!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ModeDict):
            return False
        return self._data == other._data

    def keys(self):
        """Return mode names."""
        return self._data.keys()

    def values(self):
        """Return mode objects."""
        return self._data.values()

    def items(self):
        """Return (name, mode) pairs."""
        return self._data.items()
