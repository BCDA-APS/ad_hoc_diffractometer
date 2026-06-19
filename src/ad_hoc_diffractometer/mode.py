# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
"""
mode.py — Diffraction mode / operating constraints.

A diffraction mode specifies which degrees of freedom are constrained
during a diffraction calculation and which are free.  Modes are a
first-class concern of the geometry description.

Design
------
The constraint system is built around the physical insight that specifying
(h, k, l) provides exactly 3 equations on the motor angles.  A geometry
with N real motor axes therefore has N − 3 free parameters after the Bragg
condition is satisfied.  Each free parameter must be resolved by a
constraint.

Four types of constraint cover all physically meaningful cases:

**Bisect constraint** (:class:`BisectConstraint`) is a relational sample
constraint: one named sample stage is driven to half the angle of one named
detector stage (e.g. ``eta = delta / 2`` in psic, ``omega = ttheta / 2``
in fourcv).  Both stage names are declared explicitly in the mode
definition — the diffractometer designer knows their axis names and states
them directly.  :class:`BisectConstraint` belongs to the ``"sample"``
category because it constrains a sample stage.

**Sample constraints** (:class:`SampleConstraint`) fix one sample motor
angle to a given value.  The name is a real stage name or a virtual
Eulerian angle name (for kappa geometries: ``"omega"``, ``"chi"``,
``"phi"``).

**Detector constraints** (:class:`DetectorConstraint`) fix one detector
stage angle (or the pseudo-angle ``"qaz"`` from You 1999 eq. 18).  At
most one detector constraint is allowed — fixing two detector angles would
over-constrain the scattered beam direction.

**Reference constraints** (:class:`ReferenceConstraint`) express a
condition between Q and an external reference vector n̂ (surface normal,
polarization axis, etc.) stored on the geometry.  The named options are
physical pseudo-angles from You (1999) and Lohmeier & Vlieg (1993):
``"psi"``, ``"incidence"``, ``"emergence"``, ``"specular"``, ``"naz"``.
At most one reference constraint is allowed.

Classes
-------
:class:`BisectConstraint`
    Relational sample constraint: ``sample_stage = detector_stage / 2``.
    Stage names are explicit — declared by the mode designer.

:class:`SampleConstraint`
    Constrains one sample stage to a fixed angle.

:class:`DetectorConstraint`
    Constrains one detector stage to a fixed angle, or the pseudo-angle
    ``"qaz"``.

:class:`ReferenceConstraint`
    Constrains a physical pseudo-angle between Q and the reference vector.

:class:`ConstraintSet`
    An ordered, validated collection of constraints for a diffraction mode.
    Replaces the old ``DiffractionMode`` ABC.

:class:`ModeDict`
    An ordered, guarded dict of named modes.  Validates that all entries
    are :class:`ConstraintSet` instances.

Notes
-----
Cut-points control which branch of a multi-valued angle solution is
returned (SPEC #G4 convention).  They are stored on :class:`ConstraintSet`
and applied by the forward solver.

A few specialized modes are dispatched not by their constraint composition
but by the keys present in :attr:`ConstraintSet.extras`:

- **double_diffraction** modes carry ``h2``, ``k2``, ``l2`` extras and are
  routed to a 4-DOF simultaneous Bragg solver.
- **zone** modes (You 1999 §6, SPEC ``setmode 5``) carry ``z0`` and ``z1``
  extras (each a Miller-index 3-tuple) and confine the scattering vector
  Q to the plane spanned by ``UB @ z0`` and ``UB @ z1``.  The solver
  validates that the requested (h, k, l) lies in the plane and then
  returns the bisecting solutions; off-plane requests yield an empty
  list with a warning.

This extras-driven dispatch keeps the constraint taxonomy small while
allowing modes whose distinguishing input is a parameter (a secondary
reflection or a zone plane) rather than a fixed motor angle.

The DOF rule is geometry-dependent:

- 4-circle (N=4): 1 constraint needed
- 5-circle (N=5): 2 constraints needed
- 6-circle (N=6): 3 constraints needed

:meth:`ConstraintSet.is_fully_constrained` checks this rule against the
actual geometry at solve time.

References
----------

* Busing & Levy, Acta Cryst. 22, 457-464 (1967).
* You, J. Appl. Cryst. 32, 614-623 (1999).
* Lohmeier & Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
* Walko, Ref. Module Mater. Sci. Mater. Eng. (2016), section 3.3 and 5.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING
from typing import Any

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Custom exceptions
# ---------------------------------------------------------------------------


class EwaldSphereViolation(ValueError):
    """
    Raised when the requested reflection cannot be reached at the current wavelength.

    The magnitude of the scattering vector Q exceeds the diameter of the
    Ewald sphere (``|Q| > 4π/λ``), meaning the Bragg condition cannot be
    satisfied regardless of motor angles.

    Parameters
    ----------
    q_mag : float
        Requested scattering vector magnitude in Å⁻¹.
    q_max : float
        Maximum reachable ``|Q|`` at the current wavelength (``4π/λ``) in Å⁻¹.
    wavelength : float
        Current wavelength in Å.
    """

    def __init__(self, q_mag: float, q_max: float, wavelength: float) -> None:
        self.q_mag = q_mag
        self.q_max = q_max
        self.wavelength = wavelength
        super().__init__(
            f"|Q| = {q_mag:.6g} Å⁻¹ exceeds the Ewald sphere "
            f"(max = {q_max:.6g} Å⁻¹) at λ = {wavelength} Å.  "
            "The reflection cannot be reached."
        )


class ConstraintViolation(ValueError):
    """
    Raised when a forward-calculation solution violates a declared constraint.

    This indicates either a solver error or an unimplemented virtual-angle
    constraint (e.g. a kappa Eulerian pseudoangle that has not yet been
    inverted).

    Parameters
    ----------
    solution_index : int
        Zero-based index of the violating solution in the returned list.
    constraint_repr : str
        Human-readable representation of the violated constraint.
    residual : float
        The constraint residual in degrees (non-zero = violated).
    tolerance : float
        The tolerance used for the check (from :func:`~.display.precision_atol`).
    """

    def __init__(
        self,
        solution_index: int,
        constraint_repr: str,
        residual: float,
        tolerance: float,
    ) -> None:
        self.solution_index = solution_index
        self.constraint_repr = constraint_repr
        self.residual = residual
        self.tolerance = tolerance
        super().__init__(
            f"Solution {solution_index} violates {constraint_repr} "
            f"beyond tolerance ({tolerance:.2e}°): residual = {residual:.6g}°.  "
            "This indicates a solver error or an unimplemented virtual-angle "
            "constraint (e.g. kappa Eulerian pseudoangle)."
        )


# ---------------------------------------------------------------------------
# Vocabulary constants
# ---------------------------------------------------------------------------

REFERENCE_NAMES: frozenset[str] = frozenset(
    {
        "psi",
        "incidence",
        "emergence",
        "specular",
        "naz",
        "omega",
    }
)
"""Valid reference constraint names (physical pseudo-angles).

The set of strings accepted by :class:`ReferenceConstraint` as the ``name``
argument.  Most names correspond to a physically meaningful pseudo-angle
between the scattering vector Q and the reference vector n̂; ``"omega"`` is
the SPEC ``OMEGA`` pseudo-angle, defined entirely by the diffractometer
geometry (no reference vector required):

- ``"psi"`` — azimuthal angle of n̂ about Q (You 1999, eq. 23)
- ``"incidence"`` — angle of incidence (incident beam vs. surface plane)
- ``"emergence"`` — angle of emergence (diffracted beam vs. surface plane)
- ``"specular"`` — specular reflection (relational: incidence = emergence)
- ``"naz"`` — azimuthal angle of n̂ in the lab frame
- ``"omega"`` — SPEC ``OMEGA`` pseudo-angle (``Q[6]``): angle between
  Q and the plane of the chi circle (psic family); see
  :func:`~ad_hoc_diffractometer.reference.omega_pseudo`
"""


QAZ: str = "qaz"
"""Special name for the qaz detector pseudo-angle (You 1999, eq. 18).

When used as the ``name`` argument to :class:`DetectorConstraint`, this
constrains the azimuthal angle of Q in the plane spanned by the vertical
and transverse axes:  ``tan(qaz) = tan(delta) / sin(nu)``.
Setting ``qaz = 0`` constrains scattering to the horizontal plane;
``qaz = 90°`` constrains it to the vertical plane.
"""

# Sentinel objects for extras dict values.
REQUIRED: object = object()
"""Sentinel marking a required extra input in a :class:`ConstraintSet` extras dict.

Place this sentinel as the value for any key in ``ConstraintSet.extras`` that
*must* be supplied by the caller before :meth:`~.diffractometer.AdHocDiffractometer.forward`
is called::

    cs = ConstraintSet(
        [...],
        extras={"n_hat": REQUIRED},
    )

The forward solver should check for ``REQUIRED`` sentinels and raise a
descriptive error if they have not been replaced by the caller.
"""

OPTIONAL: object = object()
"""Sentinel marking an optional extra input in a :class:`ConstraintSet` extras dict.

Place this sentinel as the value for any key in ``ConstraintSet.extras`` that
*may* be supplied by the caller but has a sensible default::

    cs = ConstraintSet(
        [...],
        extras={"azimuth": OPTIONAL},
    )
"""


# ---------------------------------------------------------------------------
# Documentation-placeholder extras keys (issue #294)
# ---------------------------------------------------------------------------

_PLACEHOLDER_EXTRA_KEYS: frozenset[str] = frozenset({"n_hat"})
"""Extras-dict keys that are *documentation hints*, not settable input slots.

Issue #294.  Some ``extras`` entries on a
:class:`ConstraintSet` — currently only ``"n_hat"`` — name a quantity
that *looks* like a settable input but is actually a documentation
placeholder.  The corresponding value lives on the diffractometer under
:attr:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.surface_normal`
or
:attr:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.azimuth`,
chosen by the mode's
:class:`ReferenceConstraint`; the ``extras`` entry is just a hint that
documents the requirement.

A user who naively assigns ``cs.extras["n_hat"] = (0, 0, 1)`` sees no
effect on :meth:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward`
and is confused.  :class:`_ExtrasDict` (below) wraps the ``extras``
dict and emits a :class:`UserWarning` whenever a placeholder key in
this set is overwritten with a non-sentinel, non-``None`` value,
directing the caller at the correct geometry attribute.

To extend the warning to a new placeholder key, add the key name to
this frozenset.  The warning is intentionally minimal — no exception
is raised — because the assignment is harmless (silently ignored by
the solver).  The goal is to surface the footgun to a user who has
misread the ``extras`` dict as a settable input bag.
"""


class _ExtrasDict(dict):
    """A dict that warns when documentation-placeholder keys are overwritten.

    Issue #294.  Behaves identically to ``dict`` except that assigning a
    non-sentinel, non-``None`` value to a key in
    :data:`_PLACEHOLDER_EXTRA_KEYS` (currently ``{"n_hat"}``) emits a
    :class:`UserWarning` directing the caller at the correct geometry
    attribute.  Sentinel writes (``REQUIRED``, ``OPTIONAL``, ``None``)
    are silent: those are legitimate construction-time states produced
    by the YAML loader and :meth:`ConstraintSet.from_dict`.

    The check is intentionally minimal — no exception, no rejection of
    the write — because the value is harmless (it's silently ignored by
    the solver).  The goal is only to surface the footgun to a user who
    has misread the ``extras`` dict as a settable input bag.
    """

    def __setitem__(self, key: str, value: Any) -> None:
        if key in _PLACEHOLDER_EXTRA_KEYS and value not in (REQUIRED, OPTIONAL, None):
            import warnings

            warnings.warn(
                f"extras[{key!r}] is a documentation placeholder, not a "
                f"settable input slot.  Assigning a value here has no "
                f"effect on forward().  Set the reference vector on the "
                f"geometry instead: use 'g.surface_normal = (h, k, l)' "
                f"for surface-mode constraints "
                f"(incidence / emergence / specular), or "
                f"'g.azimuth = (h, k, l)' for "
                f"psi / naz constraints.  See "
                f"AdHocDiffractometer.required_reference_vector to "
                f"discover which attribute the active mode needs.  "
                f"(issue #294)",
                UserWarning,
                stacklevel=2,
            )
        super().__setitem__(key, value)


# ---------------------------------------------------------------------------
# Individual constraint classes
# ---------------------------------------------------------------------------


class BisectConstraint:
    """
    Relational sample constraint: ``sample_stage = detector_stage / 2``.

    This is the bisecting condition used in standard Eulerian diffractometer
    modes.  Both stage names are declared explicitly in the mode definition.

    :class:`BisectConstraint` belongs to the ``"sample"`` category because
    it drives a sample stage, even though it depends on the value of a
    detector stage.

    Parameters
    ----------
    sample_stage : str
        Name of the sample stage to be driven to half the detector angle
        (e.g. ``"eta"`` in psic, ``"omega"`` in fourcv).
    detector_stage : str
        Name of the detector stage whose angle is halved
        (e.g. ``"delta"`` in psic, ``"ttheta"`` in fourcv).

    Examples
    --------
    >>> BisectConstraint("eta", "delta")
    BisectConstraint('eta', 'delta')
    >>> BisectConstraint("omega", "ttheta")
    BisectConstraint('omega', 'ttheta')
    """

    category: str = "sample"
    """Constraint category identifier — always ``"sample"``."""

    is_bisect: bool = True
    """Always ``True`` — identifies this as the bisecting relational constraint."""

    name: str = "bisect"
    """Constraint name — always ``"bisect"``."""

    def __init__(self, sample_stage: str, detector_stage: str) -> None:
        self._sample_stage = str(sample_stage)
        self._detector_stage = str(detector_stage)

    @property
    def sample_stage(self) -> str:
        """Name of the sample stage driven to half the detector angle."""
        return self._sample_stage

    @property
    def detector_stage(self) -> str:
        """Name of the detector stage whose angle is halved."""
        return self._detector_stage

    def evaluate(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
    ) -> float:
        """
        Return the constraint residual (zero when satisfied).

        Residual = ``angles[sample_stage] - angles[detector_stage] / 2``.

        Parameters
        ----------
        angles : dict[str, float]
            Current motor angles in degrees.
        geometry : AdHocDiffractometer
            The diffractometer (not used directly; included for protocol
            consistency with other constraint types).

        Returns
        -------
        float
            Residual in degrees.  Zero means the constraint is satisfied.
        """
        return angles[self._sample_stage] - angles[self._detector_stage] / 2.0

    def is_satisfied(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
        tol: float = 1e-6,
    ) -> bool:
        """Return True when ``|evaluate(angles, geometry)| < tol``."""
        return abs(self.evaluate(angles, geometry)) < tol

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when both named stages exist in the geometry.

        Checks that ``sample_stage`` exists in ``geometry.sample_stages``
        and ``detector_stage`` exists in ``geometry.detector_stages``.
        """
        sample_names = {s.name for s in geometry.sample_stages}
        detector_names = {s.name for s in geometry.detector_stages}
        return (
            self._sample_stage in sample_names
            and self._detector_stage in detector_names
        )

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict."""
        return {
            "type": "BisectConstraint",
            "sample_stage": self._sample_stage,
            "detector_stage": self._detector_stage,
        }

    @classmethod
    def from_dict(cls, d: dict) -> BisectConstraint:
        """Reconstruct from a dict produced by :meth:`to_dict`."""
        return cls(sample_stage=d["sample_stage"], detector_stage=d["detector_stage"])

    def __repr__(self) -> str:
        return f"BisectConstraint({self._sample_stage!r}, {self._detector_stage!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, BisectConstraint):
            return False
        return (
            self._sample_stage == other._sample_stage
            and self._detector_stage == other._detector_stage
        )

    def __hash__(self) -> int:
        return hash((self._sample_stage, self._detector_stage))


class VirtualBisectConstraint(BisectConstraint):
    """
    True virtual bisecting constraint for kappa geometries.

    Identical in structure to :class:`BisectConstraint` but
    semantically distinct: the bisecting condition is enforced on the
    **virtual** Eulerian omega pseudoangle rather than on a real motor
    angle.  On a kappa diffractometer the virtual omega is computed
    from ``(komega, kappa, kphi)`` via the Walko (2016) eq. [16]
    relations (see :func:`~ad_hoc_diffractometer.kappa.kappa_to_eulerian`).

    The constraint residual is therefore::

        residual = omega_virtual(angles) - angles[detector_stage] / 2

    where ``omega_virtual`` depends nonlinearly on the kappa motor
    triple.  This is the physically correct bisecting condition for a
    kappa diffractometer; the literal ``komega = ttheta/2`` enforced by
    :class:`BisectConstraint` is only an approximation that coincides
    with true bisecting at ``kappa = 0``.

    The dispatcher in :mod:`~ad_hoc_diffractometer.forward` routes
    modes containing a :class:`VirtualBisectConstraint` to
    :func:`~ad_hoc_diffractometer.forward._solve_bisecting_kappa_virtual`,
    which solves a 2D Newton problem in virtual ``(chi, phi)`` space
    with ``omega_virtual = ttheta/2`` enforced exactly.

    Parameters
    ----------
    sample_stage : str
        Name of the *virtual* Eulerian angle to drive (typically
        ``"omega"``).  Stored for documentation and serialization;
        the solver always enforces ``omega_virtual = ttheta/2``
        regardless of this name.
    detector_stage : str
        Name of the detector stage whose angle is halved (e.g.
        ``"ttheta"`` on kappa4cv/kappa4ch, ``"delta"`` on kappa6c).

    Examples
    --------
    >>> VirtualBisectConstraint("omega", "ttheta")
    VirtualBisectConstraint('omega', 'ttheta')

    References
    ----------
    * D. A. Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016),
      eq. [16] — kappa pseudoangle relations.
    """

    name: str = "virtual_bisect"
    """Constraint name — ``"virtual_bisect"`` (overrides ``"bisect"``)."""

    def evaluate(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
    ) -> float:
        """
        Return the constraint residual ``omega_virtual − detector/2``.

        ``omega_virtual`` is computed from the kappa motor triple
        ``(komega, kappa, kphi)`` via the geometry-aware decomposition
        :func:`~ad_hoc_diffractometer.kappa.kappa_to_eulerian_axes`,
        using the per-geometry
        :class:`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention`
        attached to ``geometry``.

        Parameters
        ----------
        angles : dict[str, float]
            Current motor angles in degrees; must include ``komega``,
            ``kappa``, ``kphi``, and the detector stage named by
            :attr:`detector_stage`.
        geometry : AdHocDiffractometer
            The diffractometer, used to retrieve the
            ``kappa_pseudo_angle_convention``.

        Returns
        -------
        float
            Residual in degrees.  Zero means true virtual bisecting is
            satisfied.

        Raises
        ------
        KeyError
            If a required motor angle is missing from ``angles``.
        ValueError
            If ``geometry.kappa_pseudo_angle_convention`` is ``None``
            (the constraint is only meaningful on kappa geometries)
            or if the kappa triple lies outside the reachable virtual
            range.
        """
        from .kappa import kappa_to_eulerian_axes

        convention = geometry.kappa_pseudo_angle_convention
        if convention is None:
            raise ValueError(
                f"VirtualBisectConstraint requires "
                f"geometry.kappa_pseudo_angle_convention to be set; "
                f"geometry {geometry.name!r} is not a kappa geometry."
            )
        komega = angles["komega"]
        kappa = angles["kappa"]
        kphi = angles["kphi"]
        omega_v, _chi_v, _phi_v = kappa_to_eulerian_axes(
            komega, kappa, kphi, convention
        )
        return omega_v - angles[self._detector_stage] / 2.0

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when the geometry is a kappa diffractometer with
        the required ``komega``, ``kappa``, ``kphi`` stages and the
        named detector stage.
        """
        if geometry.kappa_alpha_deg is None:
            return False
        sample_names = {s.name for s in geometry.sample_stages}
        if not {"komega", "kappa", "kphi"}.issubset(sample_names):  # pragma: no cover
            return False
        detector_names = {s.name for s in geometry.detector_stages}
        return self._detector_stage in detector_names

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict."""
        return {
            "type": "VirtualBisectConstraint",
            "sample_stage": self._sample_stage,
            "detector_stage": self._detector_stage,
        }

    def __repr__(self) -> str:
        return (
            f"VirtualBisectConstraint({self._sample_stage!r}, {self._detector_stage!r})"
        )

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, VirtualBisectConstraint):
            return False
        return (
            self._sample_stage == other._sample_stage
            and self._detector_stage == other._detector_stage
        )

    def __hash__(self) -> int:
        return hash(
            ("VirtualBisectConstraint", self._sample_stage, self._detector_stage)
        )


class SampleConstraint:
    """
    Constrains one sample motor angle to a fixed value.

    Parameters
    ----------
    name : str
        Stage name, or a virtual Eulerian angle name for kappa geometries
        (e.g. ``"omega"``, ``"chi"``, ``"phi"``).
    value : float
        Fixed angle in degrees.

    Examples
    --------
    >>> SampleConstraint("chi", 90.0)
    SampleConstraint('chi', 90.0)
    >>> SampleConstraint("mu", 0.0)
    SampleConstraint('mu', 0.0)
    """

    category: str = "sample"
    """Constraint category identifier — always ``"sample"``."""

    is_bisect: bool = False
    """Always ``False`` — :class:`SampleConstraint` fixes a stage, does not bisect."""

    def __init__(self, name: str, value: float) -> None:
        self._name = str(name)
        self._value = float(value)

    @property
    def name(self) -> str:
        """Constraint name: a stage name (real or virtual Eulerian)."""
        return self._name

    @property
    def value(self) -> float:
        """Constraint value: a fixed angle in degrees."""
        return self._value

    def evaluate(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
    ) -> float:
        """
        Return the constraint residual (zero when satisfied).

        Residual = ``angles[name] - value``.

        Parameters
        ----------
        angles : dict[str, float]
            Current motor angles in degrees.
        geometry : AdHocDiffractometer
            Not used; included for protocol consistency.

        Returns
        -------
        float
            Residual in degrees.  Zero means the constraint is satisfied.
        """
        return angles[self._name] - self._value

    def is_satisfied(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
        tol: float = 1e-6,
    ) -> bool:
        """Return True when ``|evaluate(angles, geometry)| < tol``."""
        return abs(self.evaluate(angles, geometry)) < tol

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when the named stage exists in the geometry, or when the
        name is a virtual Eulerian kappa pseudoangle on a kappa geometry.

        Virtual kappa angles ``"omega"``, ``"chi"``, and ``"phi"`` are
        implemented when the geometry has a ``kappa_alpha_deg`` attribute set
        (i.e. it is a kappa diffractometer) and contains a stage named
        ``"kappa"``.  The kappa inversion solver (:func:`~kappa.eulerian_to_kappa`)
        converts the virtual angle constraint to real kappa motor angles.
        """
        stage_names = {s.name for s in geometry._stages.values()}  # noqa: SLF001
        if self._name in stage_names:
            return True
        # Virtual Eulerian pseudoangles on kappa geometries
        _KAPPA_VIRTUAL = {"omega", "chi", "phi"}
        if self._name in _KAPPA_VIRTUAL:
            return geometry.kappa_alpha_deg is not None and "kappa" in stage_names
        return False

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict."""
        return {
            "type": "SampleConstraint",
            "name": self._name,
            "value": self._value,
        }

    @classmethod
    def from_dict(cls, d: dict) -> SampleConstraint:
        """Reconstruct from a dict produced by :meth:`to_dict`."""
        return cls(name=d["name"], value=d["value"])

    def __repr__(self) -> str:
        return f"SampleConstraint({self._name!r}, {self._value!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, SampleConstraint):
            return False
        return self._name == other._name and self._value == other._value

    def __hash__(self) -> int:
        return hash((self._name, self._value))


class DetectorConstraint:
    """
    Constrains one detector stage angle (or the ``"qaz"`` pseudo-angle).

    Parameters
    ----------
    name : str
        Detector stage name, or ``"qaz"`` for the azimuthal angle of Q
        (You 1999, eq. 18: ``tan(qaz) = tan(delta) / sin(nu)``).
    value : float
        Fixed angle in degrees.

    Examples
    --------
    >>> DetectorConstraint("nu", 0.0)
    DetectorConstraint('nu', 0.0)
    >>> DetectorConstraint("qaz", 90.0)
    DetectorConstraint('qaz', 90.0)
    """

    category: str = "detector"
    """Constraint category identifier — always ``"detector"``."""

    def __init__(self, name: str, value: float) -> None:
        self._name = str(name)
        self._value = float(value)

    @property
    def name(self) -> str:
        """Constraint name: a detector stage name or ``'qaz'``."""
        return self._name

    @property
    def value(self) -> float:
        """Constraint value: a fixed angle in degrees."""
        return self._value

    @property
    def is_qaz(self) -> bool:
        """True when this constrains the qaz pseudo-angle."""
        return self._name == QAZ

    def evaluate(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
    ) -> float:
        """
        Return the constraint residual (zero when satisfied).

        For a named detector stage: ``angles[name] - value``.
        For ``"qaz"``: computed from the detector angles via You (1999) eq. 18.
        """
        if self.is_qaz:
            return _qaz_residual(angles, geometry, self._value)
        return angles[self._name] - self._value

    def is_satisfied(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
        tol: float = 1e-6,
    ) -> bool:
        """Return True when ``|evaluate(angles, geometry)| < tol``."""
        return abs(self.evaluate(angles, geometry)) < tol

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when this constraint has a solver implementation.

        Fixed detector-stage constraints are implemented when the stage
        exists in the geometry.  The ``"qaz"`` pseudo-angle constraint is
        implemented when the geometry has at least two detector stages
        (an outer ``nu``-like stage and an inner ``delta``-like stage),
        which is the case for all 6-circle geometries that support
        ``lifting_detector_*`` modes.
        """
        if self.is_qaz:
            return len(geometry.detector_stages) >= 2
        det_names = {s.name for s in geometry.detector_stages}
        return self._name in det_names

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict."""
        return {
            "type": "DetectorConstraint",
            "name": self._name,
            "value": self._value,
        }

    @classmethod
    def from_dict(cls, d: dict) -> DetectorConstraint:
        """Reconstruct from a dict produced by :meth:`to_dict`."""
        return cls(name=d["name"], value=d["value"])

    def __repr__(self) -> str:
        return f"DetectorConstraint({self._name!r}, {self._value!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, DetectorConstraint):
            return False
        return self._name == other._name and self._value == other._value

    def __hash__(self) -> int:
        return hash((self._name, self._value))


class ReferenceConstraint:
    """
    Constrains a physical pseudo-angle between Q and the reference vector n̂.

    The reference vector n̂ must be stored on the geometry as
    ``geometry.surface_normal`` (preferred) or ``geometry.azimuth``
    before calling ``forward()``.

    Valid names (from You 1999 and Lohmeier & Vlieg 1993):

    ``"psi"``
        Azimuthal angle of n̂ about Q (You 1999, eq. 23).
    ``"incidence"``
        Angle of incidence: angle between the incident beam and the surface
        plane (perpendicular to n̂).
    ``"emergence"``
        Angle of emergence: angle between the diffracted beam and the
        surface plane.
    ``"specular"``
        Specular reflection: relational condition ``incidence = emergence``.
        ``value`` must be ``True``.
    ``"naz"``
        Azimuthal angle of n̂ in the lab frame (You 1999).

    Parameters
    ----------
    name : str
        One of ``"psi"``, ``"incidence"``, ``"emergence"``, ``"specular"``,
        ``"naz"``, ``"omega"``.
    value : float or bool
        Target value in degrees (or ``True`` for ``"specular"``).

    Examples
    --------
    >>> ReferenceConstraint("psi", 90.0)
    ReferenceConstraint('psi', 90.0)
    >>> ReferenceConstraint("specular", True)
    ReferenceConstraint('specular', True)
    >>> ReferenceConstraint("incidence", 5.0)
    ReferenceConstraint('incidence', 5.0)
    """

    category: str = "reference"
    """Constraint category identifier — always ``"reference"``."""

    def __init__(self, name: str, value: float | bool) -> None:
        if name not in REFERENCE_NAMES:
            raise ValueError(
                f"ReferenceConstraint name must be one of "
                f"{sorted(REFERENCE_NAMES)}; got {name!r}."
            )
        if name == "specular" and value is not True:
            raise ValueError(
                "ReferenceConstraint('specular', value): value must be True; "
                f"got {value!r}."
            )
        self._name = name
        self._value: float | bool = True if name == "specular" else float(value)  # type: ignore[arg-type]

    @property
    def name(self) -> str:
        """Constraint name (one of the reference pseudo-angle names)."""
        return self._name

    @property
    def value(self) -> float | bool:
        """Target value in degrees, or ``True`` for ``"specular"``."""
        return self._value

    def evaluate(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
    ) -> float:
        """
        Return the constraint residual (zero when satisfied).

        Not yet implemented — all reference constraints require the
        reference vector infrastructure (Issue J).  Raises
        ``NotImplementedError`` when called.
        """
        raise NotImplementedError(
            f"ReferenceConstraint('{self._name}') solver is not yet implemented. "
            "Reference constraint solvers require the reference vector "
            "infrastructure (Issue J)."
        )

    def is_satisfied(
        self,
        angles: dict[str, float],
        geometry: AdHocDiffractometer,
        tol: float = 1e-6,
    ) -> bool:
        """Return True when the constraint is satisfied (not yet implemented)."""
        raise NotImplementedError(
            f"ReferenceConstraint('{self._name}') is not yet implemented."
        )

    def has_reference_vector(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when the required reference vector is set on the geometry.

        For ``"psi"`` and ``"naz"``: requires
        :attr:`~geometry.AdHocDiffractometer.azimuth` to be set.

        For ``"incidence"``, ``"emergence"``, and ``"specular"``: requires
        :attr:`~geometry.AdHocDiffractometer.surface_normal` to be set.

        For ``"omega"``: no reference vector is required (always returns
        ``True``); the SPEC ``OMEGA`` pseudo-angle is a pure motor-frame
        quantity.

        This is a prerequisite check, separate from :meth:`is_implemented`.
        A reference constraint can have its vector set but still lack a forward
        solver — in that case ``has_reference_vector`` returns ``True`` but
        ``is_implemented`` returns ``False``.
        """
        if self._name == "omega":
            return True
        if self._name in {"psi", "naz"}:
            return geometry.azimuth is not None
        return geometry.surface_normal is not None

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when the required reference vector is set on the geometry
        and a forward solver is available for this constraint name.

        Implemented constraints:

        - ``"incidence"`` — requires :attr:`~geometry.AdHocDiffractometer.surface_normal`
        - ``"emergence"`` — requires :attr:`~geometry.AdHocDiffractometer.surface_normal`
        - ``"specular"`` — requires :attr:`~geometry.AdHocDiffractometer.surface_normal`
        - ``"psi"`` — requires :attr:`~geometry.AdHocDiffractometer.azimuth`.
          The forward solver treats ψ as a **validation filter**: for a given
          (h,k,l) and UB, ψ is a pure phi-frame quantity that is the same for
          every Bragg solution.  The solver computes the natural ψ from UB and
          the reference direction; if it matches the stored target the bisecting
          solutions are returned, otherwise an empty list is returned.  See
          issue #176 for the full analysis.
        - ``"omega"`` — SPEC ``OMEGA`` pseudo-angle (``Q[6]``).  Implemented
          for psic-family geometries (those with a sample stage named
          ``"chi"``); does not require any reference vector.  See
          :func:`~ad_hoc_diffractometer.reference.omega_pseudo`.

        Not yet implemented:

        - ``"naz"`` — no forward solver yet.
        """
        if self._name == "naz":
            return False
        if self._name == "psi":
            return geometry.azimuth is not None
        if self._name == "omega":
            # Implemented for any geometry with a chi sample stage.
            return any(s.name == "chi" for s in geometry.sample_stages)
        # incidence, emergence, specular — implemented when surface_normal is set
        return geometry.surface_normal is not None

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict."""
        return {
            "type": "ReferenceConstraint",
            "name": self._name,
            "value": self._value,
        }

    @classmethod
    def from_dict(cls, d: dict) -> ReferenceConstraint:
        """Reconstruct from a dict produced by :meth:`to_dict`."""
        return cls(name=d["name"], value=d["value"])

    def __repr__(self) -> str:
        return f"ReferenceConstraint({self._name!r}, {self._value!r})"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ReferenceConstraint):
            return False
        return self._name == other._name and self._value == other._value

    def __hash__(self) -> int:
        v = self._value
        return hash((self._name, v if isinstance(v, bool) else float(v)))


# Union type for all constraint kinds.
AnyConstraint = (
    BisectConstraint | SampleConstraint | DetectorConstraint | ReferenceConstraint
)


# ---------------------------------------------------------------------------
# ConstraintSet
# ---------------------------------------------------------------------------


class ConstraintSet:
    """
    An ordered, validated collection of constraints for a diffraction mode.

    A ``ConstraintSet`` describes one named diffraction mode by listing the
    constraints that resolve the geometry's free degrees of freedom.  It
    replaces the old ``BisectingMode`` / ``FixedAngleMode`` classes.

    **Taxonomy rules** (enforced at construction):

    - At most **one** :class:`BisectConstraint`.
    - At most **one** :class:`DetectorConstraint`.
    - At most **one** :class:`ReferenceConstraint`.
    - All remaining constraints must be :class:`SampleConstraint`.

    The total number of constraints needed (N − 3) depends on the geometry
    and is checked against the actual geometry at solve time via
    :meth:`is_fully_constrained`.

    Parameters
    ----------
    constraints : list of BisectConstraint | SampleConstraint | DetectorConstraint | ReferenceConstraint
        The constraints defining this mode, in definition order.
    computed : list of str or None
        Stage names that the solver computes (the "writable" stages in
        hklpy2 terminology).  Informational; used by documentation and
        introspection.  ``None`` means "all stages not in ``constant_stages``".
    extras : dict[str, Any] or None
        Additional parameters the mode needs beyond (h, k, l).  Use the
        :data:`REQUIRED` sentinel for mandatory inputs and :data:`OPTIONAL`
        for optional ones.  Solver output quantities are stored here as
        ``None`` placeholders and populated after ``forward()`` runs.
        Example::

            extras={
                "n_hat": REQUIRED,   # surface normal, must be supplied
                "psi": None,         # output: populated by solver
            }
    cut_points : dict[str, float] or None
        Per-stage SPEC #G4 cut-points.  A cut-point C for stage X means
        the returned angle lies in ``[C, C + 360°)``.

    Raises
    ------
    ValueError
        If more than one ``BisectConstraint`` is supplied.
    ValueError
        If more than one ``DetectorConstraint`` is supplied.
    ValueError
        If more than one ``ReferenceConstraint`` is supplied.
    ValueError
        If any constraint is not a recognized constraint type.

    Examples
    --------
    >>> # psic bisecting_vertical: BisectConstraint + mu=0 + nu=0
    >>> cs = ConstraintSet([
    ...     BisectConstraint("eta", "delta"),
    ...     SampleConstraint("mu", 0.0),
    ...     DetectorConstraint("nu", 0.0),
    ... ])
    >>> cs.has_bisect
    True
    >>> cs.bisect_stages()
    ('eta', 'delta')
    >>> cs.detector_constraint.name
    'nu'
    >>> cs.constant_stages
    ['eta', 'mu', 'nu']

    >>> # fourcv bisecting: BisectConstraint only (1 constraint for N=4)
    >>> cs4 = ConstraintSet([BisectConstraint("omega", "ttheta")])
    >>> cs4.has_bisect
    True
    >>> cs4.bisect_stages()
    ('omega', 'ttheta')
    >>> cs4.constant_stages
    ['omega']
    """

    def __init__(
        self,
        constraints: list[AnyConstraint],
        computed: list[str] | None = None,
        extras: dict[str, Any] | None = None,
        cut_points: dict[str, float] | None = None,
    ) -> None:
        # Validate taxonomy
        bisect_count = sum(1 for c in constraints if isinstance(c, BisectConstraint))
        det_count = sum(1 for c in constraints if isinstance(c, DetectorConstraint))
        ref_count = sum(1 for c in constraints if isinstance(c, ReferenceConstraint))

        if bisect_count > 1:
            raise ValueError(
                f"ConstraintSet: at most one BisectConstraint is allowed; "
                f"got {bisect_count}."
            )
        if det_count > 1:
            raise ValueError(
                f"ConstraintSet: at most one DetectorConstraint is allowed; "
                f"got {det_count}."
            )
        if ref_count > 1:
            raise ValueError(
                f"ConstraintSet: at most one ReferenceConstraint is allowed; "
                f"got {ref_count}."
            )
        for c in constraints:
            if not isinstance(
                c,
                BisectConstraint
                | SampleConstraint
                | DetectorConstraint
                | ReferenceConstraint,
            ):
                raise ValueError(
                    f"ConstraintSet: all constraints must be BisectConstraint, "
                    f"SampleConstraint, DetectorConstraint, or ReferenceConstraint; "
                    f"got {type(c).__name__!r}."
                )

        self._constraints: list[AnyConstraint] = list(constraints)
        self.computed: list[str] | None = (
            list(computed) if computed is not None else None
        )
        # Use _ExtrasDict so a later assignment to a placeholder key
        # (e.g. ``cs.extras["n_hat"] = (0, 0, 1)``) warns the caller —
        # see issue #294.  Construction-time values bypass the warning
        # path because we populate via super().__setitem__ semantics
        # (the dict() copy uses internal C-level assignment).
        self.extras: dict[str, Any] = _ExtrasDict(extras or {})
        self.cut_points: dict[str, float] = dict(cut_points or {})

    # ------------------------------------------------------------------
    # Constraint access
    # ------------------------------------------------------------------

    @property
    def constraints(self) -> list[AnyConstraint]:
        """All constraints in definition order."""
        return list(self._constraints)

    @property
    def bisect_constraint(self) -> BisectConstraint | None:
        """The :class:`BisectConstraint`, or ``None`` if not present."""
        for c in self._constraints:
            if isinstance(c, BisectConstraint):
                return c
        return None

    @property
    def sample_constraints(self) -> list[SampleConstraint]:
        """All :class:`SampleConstraint` entries (fixed-angle, not bisect)."""
        return [c for c in self._constraints if isinstance(c, SampleConstraint)]

    @property
    def detector_constraint(self) -> DetectorConstraint | None:
        """The :class:`DetectorConstraint`, or ``None`` if not present."""
        for c in self._constraints:
            if isinstance(c, DetectorConstraint):
                return c
        return None

    @property
    def reference_constraint(self) -> ReferenceConstraint | None:
        """The :class:`ReferenceConstraint`, or ``None`` if not present."""
        for c in self._constraints:
            if isinstance(c, ReferenceConstraint):
                return c
        return None

    @property
    def has_bisect(self) -> bool:
        """True when a :class:`BisectConstraint` is present."""
        return any(isinstance(c, BisectConstraint) for c in self._constraints)

    @property
    def fixed_sample_constraints(self) -> list[SampleConstraint]:
        """All :class:`SampleConstraint` entries (alias for :attr:`sample_constraints`)."""
        return self.sample_constraints

    @property
    def constant_stages(self) -> list[str]:
        """
        Stage names held constant by this constraint set, derived from the
        constraints themselves.

        Includes:

        - The sample stage from any :class:`BisectConstraint` (it is driven
          to a fixed fraction of the detector stage, so it is not free).
        - The name of every :class:`SampleConstraint` (fixed value).
        - The name of every non-qaz :class:`DetectorConstraint` (frozen at
          its declared value).

        :class:`ReferenceConstraint` entries do not map to a single stage
        name and are therefore not included.

        Returns
        -------
        list of str
            Stage names in definition order, without duplicates.
        """
        return self.constrained_stages()

    # ------------------------------------------------------------------
    # Validation
    # ------------------------------------------------------------------

    def is_fully_constrained(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when the number of constraints equals the geometry's
        free DOF count (N − 3).

        Parameters
        ----------
        geometry : AdHocDiffractometer
            The diffractometer whose ``free_dof_after_bragg`` property
            is used.

        Returns
        -------
        bool
        """
        return len(self._constraints) == geometry.free_dof_after_bragg

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return True when all constraints in this set have solver
        implementations for the given geometry.

        Parameters
        ----------
        geometry : AdHocDiffractometer

        Returns
        -------
        bool
        """
        return all(c.is_implemented(geometry) for c in self._constraints)

    def constrained_stages(
        self, geometry: AdHocDiffractometer | None = None
    ) -> list[str]:
        """
        Return the list of stage names that are constrained (fixed or
        relationally determined) by this constraint set.

        The ``geometry`` argument is accepted for API consistency but is not
        used — all stage names come directly from the constraint objects.

        Parameters
        ----------
        geometry : AdHocDiffractometer or None
            Accepted but unused.

        Returns
        -------
        list of str
        """
        names: list[str] = []
        for c in self._constraints:
            if isinstance(c, BisectConstraint):
                if c.sample_stage not in names:  # pragma: no branch
                    names.append(c.sample_stage)
            elif isinstance(c, SampleConstraint):
                if c.name not in names:  # pragma: no branch
                    names.append(c.name)
            elif isinstance(c, DetectorConstraint) and not c.is_qaz:
                if c.name not in names:  # pragma: no branch
                    names.append(c.name)
        return names

    def apply_cut_point(self, stage_name: str, angle_deg: float) -> float:
        """
        Apply the cut-point for *stage_name* to *angle_deg*.

        A cut-point C means the returned angle lies in ``[C, C + 360°)``.
        Returns *angle_deg* unchanged if no cut-point is set for this stage.
        """
        if stage_name not in self.cut_points:
            return angle_deg
        cut = self.cut_points[stage_name]
        remainder = (angle_deg - cut) % 360.0
        return cut + remainder

    # ------------------------------------------------------------------
    # Bisect stage access
    # ------------------------------------------------------------------

    def bisect_stages(self) -> tuple[str, str] | None:
        """
        Return the ``(sample_stage_name, detector_stage_name)`` pair from
        the :class:`BisectConstraint`, or ``None`` if no bisect constraint
        is present.

        The stage names are those declared explicitly in the mode definition.
        No geometry lookup or axis-geometry heuristic is performed.

        Returns
        -------
        tuple of (str, str) or None
        """
        bc = self.bisect_constraint
        if bc is None:
            return None
        return (bc.sample_stage, bc.detector_stage)

    # ------------------------------------------------------------------
    # Functional update (issue #293)
    # ------------------------------------------------------------------

    def with_constraint_values(self, **updates: float | bool) -> ConstraintSet:
        """Return a new :class:`ConstraintSet` with named constraint values replaced.

        The receiver is left unchanged.  Each keyword argument names a
        constraint by its ``.name`` attribute and supplies the new value.
        Constraint order, :attr:`computed`, :attr:`extras`, and
        :attr:`cut_points` are preserved.  The method always returns a
        fresh instance, even when ``updates`` is empty.

        Parameters
        ----------
        **updates : float or bool
            Mapping of constraint name → new value.  Each key must
            exactly match the ``.name`` attribute of an existing
            :class:`SampleConstraint`, :class:`DetectorConstraint`, or
            :class:`ReferenceConstraint` in the set.  Values are floats
            (or ``bool`` for the ``"specular"`` :class:`ReferenceConstraint`).

        Returns
        -------
        ConstraintSet
            A new instance with the named constraints replaced.

        Raises
        ------
        KeyError
            If any key in ``updates`` does not match a constraint name in
            the set.  The message lists every unknown key at once so a
            user can fix them all in a single edit.
        ValueError
            If two constraints in the receiver share the same
            ``.name``.  This indicates a malformed
            :class:`ConstraintSet` (the declarative YAML loader does not
            produce these; a user who hand-builds one with duplicate
            sample-stage names will get this error).

        Examples
        --------
        Single sample-stage value:

        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.make_geometry("psic")
        >>> g.modes["fixed_chi_vertical"] = (
        ...     g.modes["fixed_chi_vertical"].with_constraint_values(chi=45.0)
        ... )

        Multiple values at once (psic ``fixed_incidence_fixed_chi_fixed_phi``):

        >>> g.modes["fixed_incidence_fixed_chi_fixed_phi"] = (
        ...     g.modes["fixed_incidence_fixed_chi_fixed_phi"]
        ...     .with_constraint_values(chi=15.0, phi=30.0, incidence=5.0)
        ... )

        Notes
        -----
        :class:`BisectConstraint` (and :class:`VirtualBisectConstraint`)
        carry no scalar value — they are *relational* constraints
        between two stages.  This method intentionally ignores them: a
        kwarg whose key happens to equal a bisect's class-level
        ``name`` identifier (``"bisect"`` / ``"virtual_bisect"``) will
        not match and raises :class:`KeyError` like any other unknown
        key.  A bisect cannot be overridden by changing a value
        because it has none.
        """
        # Build a {name: index} map of constraints that have a settable
        # scalar value (SampleConstraint, DetectorConstraint,
        # ReferenceConstraint).  Bisect constraints are deliberately
        # excluded — they are relational and have no `.value` to
        # override (see Notes).
        settable = (SampleConstraint, DetectorConstraint, ReferenceConstraint)
        name_to_index: dict[str, int] = {}
        for idx, c in enumerate(self._constraints):
            if not isinstance(c, settable):
                continue
            cname = c.name
            if cname in name_to_index:
                raise ValueError(
                    f"with_constraint_values: this ConstraintSet contains "
                    f"two constraints both named {cname!r}; cannot resolve "
                    f"an override unambiguously."
                )
            name_to_index[cname] = idx

        unknown = sorted(k for k in updates if k not in name_to_index)
        if unknown:
            available = sorted(name_to_index)
            raise KeyError(
                f"with_constraint_values: no constraint(s) named "
                f"{unknown!r} in this ConstraintSet; available names "
                f"are {available!r}."
            )

        new_constraints: list[AnyConstraint] = list(self._constraints)
        for name, new_value in updates.items():
            idx = name_to_index[name]
            original = new_constraints[idx]
            # Each settable-value constraint is constructed with the same
            # (name, value) signature, so a single replacement pattern
            # suffices across SampleConstraint, DetectorConstraint, and
            # ReferenceConstraint.
            new_constraints[idx] = type(original)(name, new_value)

        # extras and cut_points are shallow-copied; sentinel identity
        # (REQUIRED / OPTIONAL) is preserved by the copy.
        return ConstraintSet(
            constraints=new_constraints,
            computed=list(self.computed) if self.computed is not None else None,
            extras=dict(self.extras),
            cut_points=dict(self.cut_points),
        )

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """Return a JSON-serialisable dict."""
        # extras: replace sentinel objects with string markers
        extras_serial: dict[str, Any] = {}
        for k, v in self.extras.items():
            if v is REQUIRED:
                extras_serial[k] = "__REQUIRED__"
            elif v is OPTIONAL:
                extras_serial[k] = "__OPTIONAL__"
            elif v is None:
                extras_serial[k] = None
            else:
                extras_serial[k] = v

        return {
            "type": "ConstraintSet",
            "constraints": [c.to_dict() for c in self._constraints],
            "computed": self.computed,
            "extras": extras_serial,
            "cut_points": self.cut_points,
        }

    @classmethod
    def from_dict(cls, d: dict) -> ConstraintSet:
        """Reconstruct from a dict produced by :meth:`to_dict`."""
        constraints: list[AnyConstraint] = []
        for cd in d.get("constraints", []):
            t = cd.get("type", "")
            if t == "BisectConstraint":
                constraints.append(BisectConstraint.from_dict(cd))
            elif t == "VirtualBisectConstraint":
                constraints.append(VirtualBisectConstraint.from_dict(cd))
            elif t == "SampleConstraint":
                constraints.append(SampleConstraint.from_dict(cd))
            elif t == "DetectorConstraint":
                constraints.append(DetectorConstraint.from_dict(cd))
            elif t == "ReferenceConstraint":
                constraints.append(ReferenceConstraint.from_dict(cd))
            else:
                raise ValueError(
                    f"ConstraintSet.from_dict: unknown constraint type {t!r}."
                )

        # Restore extras sentinels
        extras_raw = d.get("extras", {})
        extras: dict[str, Any] = {}
        for k, v in extras_raw.items():
            if v == "__REQUIRED__":
                extras[k] = REQUIRED
            elif v == "__OPTIONAL__":
                extras[k] = OPTIONAL
            else:
                extras[k] = v

        return cls(
            constraints=constraints,
            computed=d.get("computed"),
            extras=extras,
            cut_points=d.get("cut_points", {}),
        )

    def __repr__(self) -> str:
        parts = ", ".join(repr(c) for c in self._constraints)
        return f"ConstraintSet([{parts}])"

    def __eq__(self, other: object) -> bool:
        if not isinstance(other, ConstraintSet):
            return False
        return (
            self._constraints == other._constraints
            and self.cut_points == other.cut_points
        )

    def __len__(self) -> int:
        return len(self._constraints)


# ---------------------------------------------------------------------------
# ModeDict
# ---------------------------------------------------------------------------


class ModeDict:
    """
    An ordered, guarded dict of named diffraction modes.

    Only :class:`ConstraintSet` instances may be stored.  Keys are mode
    names (str).  Iteration follows insertion order.

    Parameters
    ----------
    modes : dict[str, ConstraintSet] or None
        Initial modes.  If ``None``, an empty dict is created.

    Raises
    ------
    TypeError
        If any value is not a :class:`ConstraintSet` instance.

    Examples
    --------
    >>> cs = ConstraintSet([BisectConstraint("omega", "ttheta")])
    >>> md = ModeDict({"bisecting": cs})
    >>> "bisecting" in md
    True
    >>> len(md)
    1
    """

    def __init__(self, modes: dict[str, ConstraintSet] | None = None) -> None:
        self._data: dict[str, ConstraintSet] = {}
        if modes:
            for name, mode in modes.items():
                self[name] = mode  # validates type

    def __setitem__(self, name: str, mode: ConstraintSet) -> None:
        if not isinstance(mode, ConstraintSet):
            raise TypeError(
                f"ModeDict values must be ConstraintSet instances; "
                f"got {type(mode).__name__!r} for key {name!r}."
            )
        self._data[name] = mode

    def __getitem__(self, name: str) -> ConstraintSet:
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


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------


def _qaz_residual(
    angles: dict[str, float],
    geometry: AdHocDiffractometer,
    target_qaz_deg: float,
) -> float:
    """
    Compute the residual for a qaz detector constraint (You 1999, eq. 18).

    ``qaz`` is the azimuthal angle of Q in the plane spanned by the
    vertical and transverse axes::

        tan(qaz) = tan(delta) / sin(nu)

    where ``nu`` is the outer (vertical-axis) detector stage angle and
    ``delta`` is the inner (transverse-axis) detector stage angle.  The
    two detector stages are identified from ``geometry.detector_stages``
    by position: the outermost stage (index 0) plays the role of ``nu``
    and the innermost stage (index -1) plays the role of ``delta``.

    Parameters
    ----------
    angles : dict[str, float]
        Current motor angles in degrees.  Must contain entries for both
        detector stage names.
    geometry : AdHocDiffractometer
        The diffractometer.  ``geometry.detector_stages`` must contain
        exactly two stages.
    target_qaz_deg : float
        Target qaz angle in degrees.

    Returns
    -------
    float
        Residual in degrees: ``qaz_computed - target_qaz_deg``.
        Zero means the constraint is satisfied.

    Raises
    ------
    ValueError
        If the geometry has fewer than two detector stages.
    """
    import math

    det_stages = geometry.detector_stages
    if len(det_stages) < 2:
        raise ValueError(
            f"_qaz_residual requires at least 2 detector stages; "
            f"geometry {geometry.name!r} has {len(det_stages)}."
        )
    # Outer detector stage (nu-like): det_stages[0]
    # Inner detector stage (delta-like): det_stages[-1]
    nu_name = det_stages[0].name
    delta_name = det_stages[-1].name

    nu_deg = angles[nu_name]
    delta_deg = angles[delta_name]

    nu_rad = math.radians(nu_deg)
    delta_rad = math.radians(delta_deg)

    sin_nu = math.sin(nu_rad)
    tan_delta = math.tan(delta_rad)

    # You (1999) eq. 18: tan(qaz) = tan(delta) / sin(nu)
    # atan2 gives the correct quadrant and handles sin(nu) = 0 edge cases.
    qaz_rad = math.atan2(tan_delta, sin_nu)
    qaz_deg = math.degrees(qaz_rad)

    return qaz_deg - target_qaz_deg
