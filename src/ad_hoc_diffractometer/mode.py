# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
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
polarisation axis, etc.) stored on the geometry.  The named options are
physical pseudo-angles from You (1999) and Lohmeier & Vlieg (1993):
``"psi"``, ``"alpha_i"``, ``"beta_out"``, ``"a_eq_b"``, ``"naz"``.  At
most one reference constraint is allowed.

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
    from .geometry import AdHocDiffractometer

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
    {"psi", "alpha_i", "beta_out", "a_eq_b", "naz"}
)
"""Valid reference constraint names (physical pseudo-angles).

The set of strings accepted by :class:`ReferenceConstraint` as the ``name``
argument.  Each name corresponds to a physically meaningful pseudo-angle
between the scattering vector Q and the reference vector n̂:

- ``"psi"`` — azimuthal angle of n̂ about Q (You 1999, eq. 23)
- ``"alpha_i"`` — angle of incidence (incident beam vs. surface plane)
- ``"beta_out"`` — angle of exit (diffracted beam vs. surface plane)
- ``"a_eq_b"`` — relational: alpha_i = beta_out (symmetric reflection)
- ``"naz"`` — azimuthal angle of n̂ in the lab frame
"""

QAZ: str = "qaz"
"""Special name for the qaz detector pseudo-angle (You 1999, eq. 18).

When used as the ``name`` argument to :class:`DetectorConstraint`, this
constrains the azimuthal angle of Q in the plane spanned by the vertical
and lateral axes:  ``tan(qaz) = tan(delta) / sin(nu)``.
Setting ``qaz = 0`` constrains scattering to the horizontal plane;
``qaz = 90°`` constrains it to the vertical plane.

.. note::
   The qaz solver is not yet implemented.
   :meth:`DetectorConstraint.is_implemented` returns ``False`` for this name.
"""

# Sentinel objects for extras dict values.
REQUIRED: object = object()
"""Sentinel marking a required extra input in a :class:`ConstraintSet` extras dict.

Place this sentinel as the value for any key in ``ConstraintSet.extras`` that
*must* be supplied by the caller before :meth:`~.geometry.AdHocDiffractometer.forward`
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
        not yet implemented (requires the Q-azimuth solver).
        """
        if self.is_qaz:
            return False  # qaz solver not yet implemented
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
    ``geometry.surface_normal`` (preferred) or ``geometry.azimuthal_reference``
    before calling ``forward()``.

    Valid names (from You 1999 and Lohmeier & Vlieg 1993):

    ``"psi"``
        Azimuthal angle of n̂ about Q (You 1999, eq. 23).
    ``"alpha_i"``
        Angle of incidence: angle between the incident beam and the surface
        plane (perpendicular to n̂).
    ``"beta_out"``
        Angle of exit: angle between the diffracted beam and the surface
        plane.
    ``"a_eq_b"``
        Relational: alpha_i = beta_out (symmetric reflection).
        ``value`` must be ``True``.
    ``"naz"``
        Azimuthal angle of n̂ in the lab frame (You 1999).

    Parameters
    ----------
    name : str
        One of ``"psi"``, ``"alpha_i"``, ``"beta_out"``, ``"a_eq_b"``,
        ``"naz"``.
    value : float or bool
        Target value in degrees (or ``True`` for ``"a_eq_b"``).

    Examples
    --------
    >>> ReferenceConstraint("psi", 90.0)
    ReferenceConstraint('psi', 90.0)
    >>> ReferenceConstraint("a_eq_b", True)
    ReferenceConstraint('a_eq_b', True)
    """

    category: str = "reference"
    """Constraint category identifier — always ``"reference"``."""

    def __init__(self, name: str, value: float | bool) -> None:
        if name not in REFERENCE_NAMES:
            raise ValueError(
                f"ReferenceConstraint name must be one of {sorted(REFERENCE_NAMES)}; "
                f"got {name!r}."
            )
        if name == "a_eq_b" and value is not True:
            raise ValueError(
                "ReferenceConstraint('a_eq_b', value): value must be True; "
                f"got {value!r}."
            )
        self._name = name
        self._value: float | bool = True if name == "a_eq_b" else float(value)  # type: ignore[arg-type]

    @property
    def name(self) -> str:
        """Constraint name (one of the reference pseudo-angle names)."""
        return self._name

    @property
    def value(self) -> float | bool:
        """Target value in degrees, or ``True`` for ``"a_eq_b"``."""
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
        :attr:`~geometry.AdHocDiffractometer.azimuthal_reference` to be set.

        For ``"alpha_i"``, ``"beta_out"``, and ``"a_eq_b"``: requires
        :attr:`~geometry.AdHocDiffractometer.surface_normal` to be set.

        This is a prerequisite check, separate from :meth:`is_implemented`.
        A reference constraint can have its vector set but still lack a forward
        solver — in that case ``has_reference_vector`` returns ``True`` but
        ``is_implemented`` returns ``False``.
        """
        if self._name in {"psi", "naz"}:
            return geometry.azimuthal_reference is not None
        return geometry.surface_normal is not None

    def is_implemented(self, geometry: AdHocDiffractometer) -> bool:
        """
        Return False — reference constraint solvers are not yet implemented.

        Will return True once the forward solvers for reference-constraint
        modes are registered (Issue J).  The reference vector infrastructure
        is available (``has_reference_vector``), but the solver that uses
        it to compute motor angles has not been written yet.
        """
        return False

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
        If any constraint is not a recognised constraint type.

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
        self.extras: dict[str, Any] = dict(extras or {})
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
    Compute the residual for a qaz detector constraint.

    ``qaz`` is the azimuthal angle of Q in the plane spanned by the
    vertical and lateral axes (You 1999, eq. 18):
    ``tan(qaz) = tan(delta) / sin(nu)``

    Not yet implemented — raises ``NotImplementedError``.
    """
    raise NotImplementedError("qaz detector constraint solver is not yet implemented.")
