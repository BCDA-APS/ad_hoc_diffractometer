"""
reflection.py — Reflection and ReflectionList for orienting reflections.

A Reflection stores the Miller indices (hkl), the diffractometer motor
angles at which the reflection was observed, the wavelength used, and
the name of the geometry on which it was recorded.

A ReflectionList owns an ordered dict of named Reflection objects and
manages the primary (or1) and secondary (or2) orienting designations.
It is held as ``AdHocDiffractometer.reflections``.

Typical usage::

    g = psic()
    rl = g.reflections                          # ReflectionList
    r1 = rl.add("r1", hkl=(1, 0, 0),
                 angles={"mu": 0, "eta": 20, ...},
                 valid_stages=set(g._stages))
    r2 = rl.add("r2", hkl=(0, 1, 0), ...)
    rl.setor1("r1")
    rl.setor2("r2")
    rl.orienting_reflections   # -> [r1, r2]
    rl["r1"]                   # -> r1
    del rl["r1"]
    rl.clear()
"""

from __future__ import annotations

from dataclasses import dataclass
from dataclasses import field


@dataclass
class Reflection:
    """
    One observed reflection used for diffractometer orientation.

    Parameters
    ----------
    name : str
        User-supplied label (e.g. ``"r1"``, ``"Si_111"``).  Set by
        ``ReflectionList.add()``.
    hkl : tuple of float
        Miller indices (h, k, l).  Need not be integers.
    angles : dict[str, float]
        Motor angles (degrees) keyed by stage name.  Keys are specific
        to the geometry named in ``geometry_name`` and are validated by
        ``ReflectionList.add()`` before the reflection is stored.
    wavelength : float or None
        Wavelength in Å.  If None the geometry's wavelength is assumed.
    geometry_name : str or None
        Name of the geometry on which this reflection was recorded.
        Angle keys are only meaningful for that geometry.

    Notes
    -----
    Reflections are geometry-specific: angle keys from a ``psic`` geometry
    cannot be applied to ``kappa6c`` because the stage names differ.
    ``ReflectionList.add()`` validates keys and records ``geometry_name``
    so cross-geometry misuse can be detected.
    """

    name: str
    hkl: tuple[float, float, float]
    angles: dict[str, float]
    wavelength: float | None = field(default=None)
    geometry_name: str | None = field(default=None)

    def __post_init__(self) -> None:
        h, k, l = self.hkl  # noqa: E741
        self.hkl = (float(h), float(k), float(l))
        self.angles = {n: float(v) for n, v in self.angles.items()}
        if self.wavelength is not None:
            self.wavelength = float(self.wavelength)
            if self.wavelength <= 0:
                raise ValueError(
                    f"Reflection wavelength must be > 0 Å; got {self.wavelength}."
                )

    def __eq__(self, other: object) -> bool:
        """
        True if all fields match.

        Two reflections from different geometries are never equal even if
        hkl and angle values coincide, because the angle keys carry
        different physical meanings.
        """
        if not isinstance(other, Reflection):
            return NotImplemented
        return (
            self.name == other.name
            and self.hkl == other.hkl
            and self.angles == other.angles
            and self.wavelength == other.wavelength
            and self.geometry_name == other.geometry_name
        )

    def __repr__(self) -> str:
        wl = f"{self.wavelength} Å" if self.wavelength is not None else "not set"
        return (
            f"Reflection(name={self.name!r}, hkl={self.hkl}, "
            f"angles={self.angles}, wavelength={wl}, "
            f"geometry_name={self.geometry_name!r})"
        )


class ReflectionList:
    """
    Ordered, named collection of Reflection objects with or1/or2 management.

    Held as ``AdHocDiffractometer.reflections``.  All mutation methods
    validate that angle keys belong to the owning geometry.

    Parameters
    ----------
    geometry_name : str
        Name of the owning ``AdHocDiffractometer``.  Stored on every
        reflection added through this list.
    valid_stages : set of str
        Stage names of the owning geometry.  Used to validate angle keys.

    Attributes
    ----------
    geometry_name : str
    valid_stages : set of str

    Examples
    --------
    >>> rl = ReflectionList(geometry_name="psic",
    ...                     valid_stages={"mu", "eta", "chi", "phi",
    ...                                   "nu", "delta"})
    >>> r = rl.add("r1", hkl=(1, 0, 0), angles={"mu": 0, "eta": 20})
    >>> rl.setor1("r1")
    >>> rl.orienting_reflections
    [Reflection(name='r1', ...)]
    """

    def __init__(self, geometry_name: str, valid_stages: set[str]) -> None:
        self.geometry_name = geometry_name
        self.valid_stages = set(valid_stages)
        self._data: dict[str, Reflection] = {}
        self._or1_name: str | None = None
        self._or2_name: str | None = None

    # ------------------------------------------------------------------
    # Dict-like interface
    # ------------------------------------------------------------------

    def __getitem__(self, name: str) -> Reflection:
        return self._data[name]

    def __delitem__(self, name: str) -> None:
        del self._data[name]
        if self._or1_name == name:
            self._or1_name = None
        if self._or2_name == name:
            self._or2_name = None

    def __contains__(self, name: object) -> bool:
        return name in self._data

    def __len__(self) -> int:
        return len(self._data)

    def __iter__(self):
        return iter(self._data)

    def __repr__(self) -> str:
        return (
            f"ReflectionList(geometry_name={self.geometry_name!r}, "
            f"reflections={list(self._data)})"
        )

    def keys(self):
        return self._data.keys()

    def values(self):
        return self._data.values()

    def items(self):
        return self._data.items()

    # ------------------------------------------------------------------
    # Mutation
    # ------------------------------------------------------------------

    def add(
        self,
        name: str,
        hkl: tuple[float, float, float],
        angles: dict[str, float],
        wavelength: float | None = None,
    ) -> Reflection:
        """
        Add a named reflection, validate angle keys, and return it.

        Parameters
        ----------
        name : str
            Unique label.
        hkl : tuple of float
            Miller indices.
        angles : dict[str, float]
            Motor angles; keys must be stage names of this geometry.
        wavelength : float or None
            Wavelength in Å; stored as-is (caller supplies the geometry's
            current wavelength if desired).

        Raises
        ------
        ValueError
            If ``name`` already exists or any angle key is not a valid stage.
        """
        if name in self._data:
            raise ValueError(
                f"A reflection named {name!r} already exists. Remove it first."
            )
        unknown = set(angles) - self.valid_stages
        if unknown:
            raise ValueError(
                f"Angle key(s) {sorted(unknown)} are not stage names in "
                f"geometry {self.geometry_name!r}. "
                f"Valid stages: {sorted(self.valid_stages)}."
            )
        r = Reflection(
            name=name,
            hkl=hkl,
            angles=dict(angles),
            wavelength=wavelength,
            geometry_name=self.geometry_name,
        )
        self._data[name] = r
        return r

    def remove(self, name: str) -> None:
        """
        Remove the named reflection and clear any or1/or2 designation it held.

        Raises
        ------
        KeyError
            If no reflection with that name exists.
        """
        del self[name]

    def clear(self) -> None:
        """Remove all reflections and clear or1/or2 designations."""
        self._data.clear()
        self._or1_name = None
        self._or2_name = None

    # ------------------------------------------------------------------
    # or1 / or2 management
    # ------------------------------------------------------------------

    def _resolve(self, reflection: str | Reflection) -> str:
        """Return the name of a reflection given a name or object."""
        name = reflection.name if isinstance(reflection, Reflection) else reflection
        if name not in self._data:
            raise KeyError(f"No reflection named {name!r}. Add it first.")
        return name

    def setor1(self, reflection: str | Reflection) -> None:
        """
        Designate a reflection as primary (or1).

        If it was previously the secondary (or2), that designation is cleared.
        The previous or1 remains in the list without any designation.
        """
        name = self._resolve(reflection)
        if name == self._or2_name:
            self._or2_name = None
        self._or1_name = name

    def setor2(self, reflection: str | Reflection) -> None:
        """
        Designate a reflection as secondary (or2).

        If it was previously the primary (or1), that designation is cleared.
        The previous or2 remains in the list without any designation.
        """
        name = self._resolve(reflection)
        if name == self._or1_name:
            self._or1_name = None
        self._or2_name = name

    @property
    def orienting_reflections(self) -> list[Reflection]:
        """
        Ordered list of designated orienting reflections.

        Returns ``[]``, ``[or1]``, or ``[or1, or2]``.  A secondary
        without a primary is returned as a single-element list containing
        the secondary.
        """
        result: list[Reflection] = []
        if self._or1_name is not None:
            result.append(self._data[self._or1_name])
        if self._or2_name is not None:
            result.append(self._data[self._or2_name])
        return result
