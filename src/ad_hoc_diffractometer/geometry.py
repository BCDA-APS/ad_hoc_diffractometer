"""
geometry.py — AdHocDiffractometer class.

Describes a diffractometer as an ordered collection of rotary stages.
The stacking order is encoded via the parent attribute of each Stage.
"""

import numpy as np

from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .reflection import Reflection
from .reflection import ReflectionList
from .stage import Stage


class AdHocDiffractometer:
    """
    Description of a diffractometer as an ordered collection of rotary stages.

    The geometry is defined by a list of Stage objects.  The stacking order
    (which stage sits on which) is encoded by the parent attribute of each
    Stage: a stage with parent=None sits on the fixed lab frame; all other
    stages sit on the stage named by their parent.

    The class supports arbitrary geometries including:
      - You (1999) psic  4S+2D  six-circle
      - Lohmeier & Vlieg (1993) sixc  six-circle
      - Busing & Levy (1967) fourc  four-circle Eulerian
      - Kappa four-circle

    The basis vectors for the lab frame are specified at construction time,
    allowing the same physical geometry to be described in different
    coordinate conventions.

    Parameters
    ----------
    name : str
        Name of the diffractometer geometry (e.g. 'psic', 'sixc', 'fourc').
    stages : list of Stage
        All stages in the geometry.  Order within the list does not determine
        the stacking; the parent attribute of each Stage does.
    basis : dict, optional
        Mapping from physical direction names to Cartesian basis vectors.
        Default is the You (1999) convention:
            {'vertical': XHAT, 'longitudinal': YHAT, 'lateral': ZHAT}
    description : str, optional
        Free-text description of the geometry.
    wavelength : float or None, optional
        X-ray or neutron wavelength in Å.  Must match the units used for
        unit cell edge lengths.  Default is None (unset).  Must be > 0
        if provided.
    kappa_alpha_deg : float or None, optional
        Kappa tilt angle in degrees for kappa geometries (kappa4cv,
        kappa4ch, kappa6c).  None for non-kappa geometries.  Set by
        the kappa factory functions; not intended to be changed after
        construction.

    Attributes
    ----------
    sample_stages : list of Stage
        Stages with role='sample', in stacking order (floor first).
    detector_stages : list of Stage
        Stages with role='detector', in stacking order (floor first).
    wavelength : float or None
        Wavelength in Å, or None if not set.
    kappa_alpha_deg : float or None
        Kappa tilt angle in degrees, or None for non-kappa geometries.
    """

    DEFAULT_BASIS = {
        "vertical": XHAT,
        "longitudinal": YHAT,
        "lateral": ZHAT,
    }

    def __init__(
        self,
        name: str,
        stages: list[Stage],
        basis: dict | None = None,
        description: str = "",
        wavelength: float | None = None,
        kappa_alpha_deg: float | None = None,
    ):
        self.name = name
        self.description = description
        self.basis = basis if basis is not None else dict(self.DEFAULT_BASIS)
        self.wavelength = wavelength  # validated via property setter
        self.kappa_alpha_deg = kappa_alpha_deg

        # Validate basis vectors
        self._check_basis()

        # Store stages in a dict keyed by name for fast lookup
        self._stages = {s.name: s for s in stages}
        if len(self._stages) != len(stages):
            raise ValueError("Stage names must be unique.")

        # Validate parent references
        for s in stages:
            if s.parent is not None and s.parent not in self._stages:
                raise ValueError(
                    f"Stage {s.name!r} has parent {s.parent!r} "
                    f"which is not in the stage list."
                )

        # Check for cycles in the parent graph
        self._check_no_cycles()

        # Build ordered lists per role
        self.sample_stages = self._ordered_stages("sample")
        self.detector_stages = self._ordered_stages("detector")

        # Ordered collection of named orienting reflections with or1/or2 management
        self.reflections: ReflectionList = ReflectionList(
            geometry_name=self.name,
            valid_stages=set(self._stages),
        )

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _check_basis(self) -> None:
        """
        Validate that the basis contains exactly three mutually orthogonal,
        non-zero, 3-dimensional vectors.

        The basis is a labelled dict; its keys have no defined ordering, so
        right-handedness cannot be checked here (the cross product of two
        vectors depends on which is 'first').  Right-handedness is the
        caller's responsibility and should be verified by the geometry
        factory functions or the user.

        Raises
        ------
        ValueError
            If the basis does not have exactly three entries, if any vector
            is not 3-dimensional or is zero, or if any two vectors are not
            mutually orthogonal (normalised dot product exceeds tolerance).
        """
        vecs = list(self.basis.values())
        names = list(self.basis.keys())

        if len(vecs) != 3:
            raise ValueError(f"Basis must contain exactly 3 vectors; got {len(vecs)}.")

        for name, v in zip(names, vecs, strict=False):
            v = np.asarray(v, dtype=float)
            if v.shape != (3,):
                raise ValueError(
                    f"Basis vector {name!r} must be 3-dimensional; got shape {v.shape}."
                )
            if np.linalg.norm(v) == 0.0:
                raise ValueError(f"Basis vector {name!r} must be non-zero.")

        atol = 1e-10
        for i, (n1, v1) in enumerate(zip(names, vecs, strict=False)):
            for n2, v2 in zip(names[i + 1 :], vecs[i + 1 :], strict=False):
                v1n = np.asarray(v1, dtype=float) / np.linalg.norm(v1)
                v2n = np.asarray(v2, dtype=float) / np.linalg.norm(v2)
                dot = np.dot(v1n, v2n)
                if abs(dot) > atol:
                    raise ValueError(
                        f"Basis vectors {n1!r} and {n2!r} are not orthogonal "
                        f"(normalised dot product = {dot:.6g}, tolerance {atol})."
                    )

    def _check_no_cycles(self) -> None:
        """Raise ValueError if the parent graph contains a cycle."""
        for start in self._stages:
            visited: set[str] = set()
            node = self._stages[start].parent
            while node is not None:
                if node in visited:
                    raise ValueError(
                        f"Cycle detected in parent chain starting from {start!r}."
                    )
                visited.add(node)
                node = self._stages[node].parent

    def _ordered_stages(self, role: str) -> list[Stage]:
        """
        Return stages of the given role in stacking order (floor-most first).

        Uses a topological sort on the parent graph restricted to stages of
        the requested role.
        """
        role_stages = [s for s in self._stages.values() if s.role == role]
        remaining = list(role_stages)
        ordered: list[Stage] = []
        role_names = {s.name for s in role_stages}

        max_iter = len(remaining) + 1
        while remaining:
            max_iter -= 1
            if max_iter < 0:
                raise RuntimeError(
                    "Could not determine stacking order; check parent chain."
                )
            for s in remaining:
                if s.parent is None or s.parent not in role_names:
                    ordered.append(s)
                    remaining.remove(s)
                    role_names.discard(s.name)
                    break

        return ordered

    # ------------------------------------------------------------------
    # Wavelength
    # ------------------------------------------------------------------

    @property
    def wavelength(self) -> float | None:
        """Wavelength in Å, or None if not set."""
        return self._wavelength

    @wavelength.setter
    def wavelength(self, value: float | None) -> None:
        if value is not None:
            value = float(value)
            if value <= 0:
                raise ValueError(f"wavelength must be > 0 Å; got {value}.")
        self._wavelength = value

    # ------------------------------------------------------------------
    # Kappa alpha
    # ------------------------------------------------------------------

    @property
    def kappa_alpha_deg(self) -> float | None:
        """
        Kappa tilt angle in degrees, or None for non-kappa geometries.

        This is the angle between the kappa rotation axis and the vertical
        axis (toward the lateral axis).  Typical value: 50 deg.
        Set by kappa factory functions (kappa4cv, kappa4ch, kappa6c).
        """
        return self._kappa_alpha_deg

    @kappa_alpha_deg.setter
    def kappa_alpha_deg(self, value: float | None) -> None:
        if value is not None:
            value = float(value)
        self._kappa_alpha_deg = value

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def stage(self, name: str) -> Stage:
        """Return the Stage with the given name."""
        return self._stages[name]

    def set_angle(self, name: str, angle_deg: float) -> None:
        """Set the angle of a named stage in degrees."""
        self._stages[name].angle = angle_deg

    def check_limits(self, **angles: float) -> None:
        """
        Verify that all supplied angles are within their stage limits.

        Parameters
        ----------
        **angles : float
            Keyword arguments mapping stage name to angle in degrees,
            e.g. check_limits(mu=10.0, eta=30.0).

        Raises
        ------
        KeyError
            If a supplied stage name does not exist in this geometry.
        ValueError
            If any supplied angle is outside its stage's limits.  The
            message names every out-of-range stage, its angle, and its
            limits, so all violations are reported in one call.
        """
        violations: list[str] = []
        for name, angle in angles.items():
            stage = self._stages[name]  # raises KeyError for unknown stages
            if not stage.in_limits(angle):
                lo, hi = stage.limits
                violations.append(
                    f"  {name!r}: angle {angle} deg is outside limits [{lo}, {hi}]"
                )
        if violations:
            raise ValueError(
                "The following stages have angles outside their limits:\n"
                + "\n".join(violations)
            )

    # ------------------------------------------------------------------
    # Reflections — convenience wrapper around self.reflections
    # ------------------------------------------------------------------

    def add_reflection(
        self,
        name: str,
        hkl: tuple[float, float, float],
        angles: dict[str, float],
        wavelength: float | None = None,
    ) -> Reflection:
        """
        Add a named orienting reflection recorded on this geometry.

        Delegates to ``self.reflections.add()``.  If ``wavelength`` is
        None the geometry's current wavelength is inherited.

        Parameters
        ----------
        name : str
            Unique label (e.g. ``"r1"``, ``"Si_111"``).
        hkl : tuple of float
            Miller indices (h, k, l).
        angles : dict[str, float]
            Motor angles in degrees keyed by stage name.  Keys must be
            stage names of this geometry.
        wavelength : float or None
            Wavelength in Å.  If None, inherits ``self.wavelength``.

        Returns
        -------
        Reflection
        """
        wl = wavelength if wavelength is not None else self._wavelength
        return self.reflections.add(name=name, hkl=hkl, angles=angles, wavelength=wl)

    def sample_rotation_matrix(self) -> np.ndarray:
        """
        Compute the total sample rotation matrix Z = R_floor * ... * R_top.

        Stages are applied in stacking order: the floor-most stage is applied
        first (leftmost in the product), the innermost stage last.

        Returns
        -------
        Z : numpy.ndarray, shape (3, 3)
        """
        Z = np.eye(3)
        for s in self.sample_stages:
            Z = s.rotation_matrix() @ Z
        return Z

    def detector_rotation_matrix(self) -> np.ndarray:
        """
        Compute the total detector rotation matrix in stacking order.

        Returns
        -------
        D : numpy.ndarray, shape (3, 3)
        """
        D = np.eye(3)
        for s in self.detector_stages:
            D = s.rotation_matrix() @ D
        return D

    def summary(self) -> None:
        """Print a human-readable summary of the geometry."""
        from .axes import axis_label
        from .display import fmt

        if self._wavelength is not None:
            wl_str = f"{fmt(self._wavelength)} Å"
        else:
            wl_str = "not set"

        lines = [
            f"AdHocDiffractometer: {self.name}",
            f"  {self.description}",
            f"  λ = {wl_str}",
            "",
            "  Basis vectors:",
        ]
        for direction, vec in self.basis.items():
            lines.append(f"    {direction:12s} -> {axis_label(vec)}")
        lines.append("")
        lines.append("  Sample stages (floor first):")
        for s in self.sample_stages:
            lines.append(
                f"    {s.name:8s}  axis={axis_label(s.axis):4s}  "
                f"parent={s.parent}  angle={s.angle} deg"
            )
        lines.append("")
        lines.append("  Detector stages (floor first):")
        for s in self.detector_stages:
            lines.append(
                f"    {s.name:8s}  axis={axis_label(s.axis):4s}  "
                f"parent={s.parent}  angle={s.angle} deg"
            )
        print("\n".join(lines))

    def __repr__(self) -> str:
        return (
            f"AdHocDiffractometer(name={self.name!r}, "
            f"stages={list(self._stages.keys())})"
        )
