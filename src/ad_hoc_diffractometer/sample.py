"""
sample.py — Sample class and SampleDict for crystallographic sample management.

A Sample bundles together everything that is specific to one crystal
mounted on the diffractometer:

  - a name (also used as the dict key in AdHocDiffractometer.samples)
  - a Lattice (crystal structure → B matrix)
  - a ReflectionList (geometry-specific orienting reflections)
  - U and UB matrices (orientation, computed later)

A SampleDict is a guarded ordered dict that:
  - only accepts Sample values (no None, no arbitrary objects)
  - prevents removing or replacing the currently active sample unless a
    different sample has first been selected as active

The diffractometer owns a SampleDict and keeps one sample active.
All reflection operations via AdHocDiffractometer.add_reflection()
target the active sample's ReflectionList.

Typical usage::

    g = psic()
    g.sample                        # -> Sample(name='test', ...)
    g.samples                       # -> SampleDict with 'test'
    g.add_sample("silicon", Lattice(a=5.4310))
    g.sample = "silicon"            # set active sample by name
    g.add_reflection("r1", hkl=(1,1,1), angles={...})  # -> silicon's reflections
"""

from __future__ import annotations

import numpy as np

from .lattice import Lattice
from .reflection import ReflectionList

# Default sample name and lattice used when a new geometry is constructed.
_DEFAULT_SAMPLE_NAME = "test"
_DEFAULT_LATTICE = Lattice(a=1.0)  # cubic, 1 Å


class SampleDict:
    """
    Guarded ordered dict of Sample objects.

    Enforces two invariants:
      1. Every value must be a Sample instance — no None, no arbitrary objects.
      2. The currently active sample cannot be removed or replaced.

    The active sample name is tracked via a reference to a one-element list
    ``[active_name]`` shared with the owning AdHocDiffractometer, so changes
    to the active selection are always visible here.

    Direct reassignment of the underlying container
    (``g._samples = something``) is blocked by making ``_samples`` a
    read-only property on AdHocDiffractometer.
    """

    def __init__(self, active_ref: list[str]) -> None:
        # active_ref is a one-element list [name] shared with the geometry;
        # mutating active_ref[0] updates the active selection everywhere.
        self._active_ref = active_ref
        self._data: dict[str, Sample] = {}

    # ------------------------------------------------------------------
    # Active-name helpers
    # ------------------------------------------------------------------

    @property
    def _active_name(self) -> str:
        return self._active_ref[0]

    def _guard_active(self, name: str, action: str) -> None:
        if name == self._active_name:
            raise ValueError(
                f"Cannot {action} the active sample {name!r}. "
                f"Select a different sample first."
            )

    # ------------------------------------------------------------------
    # Type guard
    # ------------------------------------------------------------------

    @staticmethod
    def _guard_type(value: object) -> None:
        if not isinstance(value, Sample):
            raise TypeError(
                f"SampleDict only accepts Sample instances; got {type(value).__name__!r}."
            )

    # ------------------------------------------------------------------
    # Dict-like interface
    # ------------------------------------------------------------------

    def __setitem__(self, name: str, value: object) -> None:
        self._guard_type(value)
        if name in self._data and name == self._active_name:
            self._guard_active(name, "replace")
        self._data[name] = value  # type: ignore[assignment]

    def __delitem__(self, name: str) -> None:
        if name not in self._data:
            raise KeyError(name)
        self._guard_active(name, "remove")
        del self._data[name]

    def __getitem__(self, name: str) -> Sample:
        return self._data[name]

    def __contains__(self, name: object) -> bool:
        return name in self._data

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self) -> str:
        return f"SampleDict({list(self._data)})"

    def clear(self) -> None:
        raise ValueError(
            "SampleDict.clear() is not permitted; the dict must always "
            "contain at least the active sample."
        )

    def pop(self, name: str, *args) -> Sample:
        if name not in self._data:
            if args:
                return args[0]
            raise KeyError(name)
        self._guard_active(name, "pop")
        return self._data.pop(name)

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()


class Sample:
    """
    One crystal sample mounted on a diffractometer.

    Parameters
    ----------
    name : str
        Unique label, matching the key in ``AdHocDiffractometer.samples``.
    lattice : Lattice
        Crystal lattice; provides the B matrix.
    reflections : ReflectionList
        Geometry-specific reflection collection.  Must be constructed
        with the owning geometry's stage names as ``valid_stages``.
    U : numpy.ndarray or None
        3×3 crystal orientation matrix.  None until computed.
    UB : numpy.ndarray or None
        3×3 UB = U @ B matrix.  None until computed.

    Notes
    -----
    U and UB are set by the U/UB computation routines (Priority 2 items
    #5 and #6).  They are stored here so that each sample carries its
    own orientation independently.
    """

    def __init__(
        self,
        name: str,
        lattice: Lattice,
        reflections: ReflectionList,
        U: np.ndarray | None = None,
        UB: np.ndarray | None = None,
        parent=None,
    ) -> None:
        self.name = name
        self.lattice = lattice
        self.reflections = reflections
        self.U = U
        self.UB = UB
        self.parent = parent  # owning AdHocDiffractometer, or None

    def __repr__(self) -> str:
        u_str = "set" if self.U is not None else "None"
        ub_str = "set" if self.UB is not None else "None"
        parent_str = self.parent.name if self.parent is not None else "(no parent)"
        return (
            f"Sample(name={self.name!r}, geometry={parent_str!r}, "
            f"lattice={self.lattice!r}, "
            f"reflections={len(self.reflections)} reflection(s), "
            f"U={u_str}, UB={ub_str})"
        )

    def __eq__(self, other: object) -> bool:
        """
        True if name, lattice, U, and UB all match.

        ``parent`` is excluded from equality — two samples with the same
        content owned by different geometries are considered equal.
        """
        if not isinstance(other, Sample):
            return NotImplemented
        u_eq = (self.U is None and other.U is None) or (
            self.U is not None
            and other.U is not None
            and np.array_equal(self.U, other.U)
        )
        ub_eq = (self.UB is None and other.UB is None) or (
            self.UB is not None
            and other.UB is not None
            and np.array_equal(self.UB, other.UB)
        )
        return (
            self.name == other.name and self.lattice == other.lattice and u_eq and ub_eq
        )

    def to_dict(self) -> dict:
        """
        Return a JSON-serialisable ``dict`` representing this sample.

        Returns
        -------
        dict
            Keys: ``"name"`` (str), ``"lattice"`` (lattice dict),
            ``"reflections"`` (reflection-list dict), ``"U"`` (3×3 list or
            None), ``"UB"`` (3×3 list or None).
        """
        return {
            "name": self.name,
            "lattice": self.lattice.to_dict(),
            "reflections": self.reflections.to_dict(),
            "U": self.U.tolist() if self.U is not None else None,
            "UB": self.UB.tolist() if self.UB is not None else None,
        }

    @classmethod
    def from_dict(cls, d: dict, parent=None) -> Sample:
        """
        Reconstruct a :class:`Sample` from a dict produced by :meth:`to_dict`.

        Parameters
        ----------
        d : dict
            Must contain ``"name"``, ``"lattice"``, ``"reflections"``.
            ``"U"`` and ``"UB"`` are optional (default ``None``).
        parent : AdHocDiffractometer or None
            Geometry to attach as the parent back-reference.

        Returns
        -------
        Sample
        """
        import numpy as np

        from .lattice import Lattice
        from .reflection import ReflectionList

        lattice = Lattice.from_dict(d["lattice"])
        reflections = ReflectionList.from_dict(d["reflections"])
        U = np.array(d["U"], dtype=float) if d.get("U") is not None else None
        UB = np.array(d["UB"], dtype=float) if d.get("UB") is not None else None
        return cls(
            name=d["name"],
            lattice=lattice,
            reflections=reflections,
            U=U,
            UB=UB,
            parent=parent,
        )
