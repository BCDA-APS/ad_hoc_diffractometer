# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
geometry.py — AdHocDiffractometer class.

Describes a diffractometer as an ordered collection of rotary stages.
The stacking order is encoded via the parent attribute of each Stage.
"""

import builtins
import logging

import numpy as np

from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .mode import ConstraintSet
from .mode import ModeDict
from .reflection import Reflection
from .reflection import ReflectionList
from .sample import _DEFAULT_LATTICE
from .sample import _DEFAULT_SAMPLE_NAME
from .sample import Sample
from .sample import SampleDict
from .stage import Stage

logger = logging.getLogger(__name__)


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
        Default is the You (1999) convention
        ``{'vertical': XHAT, 'longitudinal': YHAT, 'lateral': ZHAT}``.
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
    azimuthal_reference : tuple of float or None, optional
        Azimuthal reference direction as Miller indices (h, k, l).  Used
        by :meth:`psi` to compute the azimuthal angle ψ.  ``None`` (default)
        means no reference is set.  Must be a non-zero vector.
    modes : dict[str, ConstraintSet] or ModeDict or None, optional
        Named diffraction modes available for this geometry.  Keys are
        mode names (str); values are :class:`~mode.DiffractionMode`
        instances.  ``None`` (default) means no modes are declared.
    default_mode : str or None, optional
        Name of the mode that is active at construction time.  Must be a
        key of ``modes`` if both are supplied.  ``None`` (default) means no
        active mode is set (``mode_name`` returns ``None``).
    cut_points : dict[str, float] or None, optional
        Geometry-level SPEC #G4 cut-points per stage name.  These are the
        default cut-points when no mode-level cut-point overrides them.
        A cut-point C for stage X means the returned angle lies in
        ``[C, C + 360°)``.  ``None`` is equivalent to an empty dict.

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
        azimuthal_reference: tuple[float, float, float] | None = None,
        modes: dict | ModeDict | None = None,
        default_mode: str | None = None,
        cut_points: dict[str, float] | None = None,
    ):
        self.name = name
        self.description = description
        self.basis = basis if basis is not None else dict(self.DEFAULT_BASIS)
        self._source_type: str = "xray"  # default; validated via property setter
        self.wavelength = wavelength  # validated via property setter
        self.kappa_alpha_deg = kappa_alpha_deg
        self.azimuthal_reference = azimuthal_reference  # validated via property setter
        self._surface_normal: tuple[float, float, float] | None = None
        self._detector_distance: float | None = None
        self._detector_tilt: float | None = None
        self._detector_offset: tuple[float, float] | None = None
        self._inclination_matrix: np.ndarray = np.eye(3)

        # Diffraction modes
        if isinstance(modes, ModeDict):
            self._modes = modes
        else:
            self._modes = ModeDict(modes)  # accepts None or plain dict

        # Validate and set the active mode name
        if default_mode is not None and default_mode not in self._modes:
            available = sorted(self._modes.keys())
            raise ValueError(
                f"default_mode {default_mode!r} is not in the modes dict. "
                f"Available modes: {available}."
            )
        self._mode_name: str | None = default_mode

        # Geometry-level cut-points (SPEC #G4)
        self.cut_points: dict[str, float] = dict(cut_points or {})

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

        # Samples: a SampleDict guarded against removing the active sample
        # or inserting non-Sample values.  The active name is shared via a
        # one-element list so SampleDict always sees the current selection.
        self._active_ref: list[str] = [_DEFAULT_SAMPLE_NAME]
        self.__samples = SampleDict(active_ref=self._active_ref)
        _default = Sample(
            name=_DEFAULT_SAMPLE_NAME,
            lattice=_DEFAULT_LATTICE,
            reflections=ReflectionList(
                geometry_name=self.name,
                valid_stages=set(self._stages),
            ),
            parent=self,
        )
        self.__samples._data[_DEFAULT_SAMPLE_NAME] = _default

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
            if max_iter < 0:  # pragma: no cover
                raise RuntimeError(
                    "Could not determine stacking order; check parent chain."
                )
            for s in remaining:  # pragma: no branch
                if s.parent is None or s.parent not in role_names:
                    ordered.append(s)
                    remaining.remove(s)
                    role_names.discard(s.name)
                    break

        return ordered

    def stages_by_role(self, role: str) -> list[Stage]:
        """
        Return all stages with the given role, in stacking order (base first).

        The two conventional roles are ``"sample"`` and ``"detector"``,
        which are also available as the :attr:`sample_stages` and
        :attr:`detector_stages` attributes.  This method accepts any
        arbitrary role string, making it useful for non-standard components
        such as analysers, polarisers, slits, or azimuthal spinners.

        Parameters
        ----------
        role : str
            The role string to filter on (case-sensitive).

        Returns
        -------
        list of Stage
            Stages with the requested role, in stacking order (base-most
            first).  An empty list if no stages have the requested role.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> import numpy as np
        >>> from ad_hoc_diffractometer import Stage
        >>> g = ahd.fourcv()
        >>> analyser = Stage("analyser", np.array([1., 0., 0.]), role="analyser")
        >>> g._stages["analyser"] = analyser
        >>> [s.name for s in g.stages_by_role("analyser")]
        ['analyser']
        >>> g.stages_by_role("unknown")
        []
        """
        return self._ordered_stages(role)

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

    @property
    def source_type(self) -> str:
        """
        Radiation source type: ``"xray"`` (default) or ``"neutron"``.

        Controls which energy formula is used by :attr:`energy` and
        :attr:`energy_units`:

        ``"xray"``
            Photon energy: E (keV) = hc/λ = 12.39842 / λ (Å).

        ``"neutron"``
            Fixed-wavelength **reactor** neutron kinetic energy via the de
            Broglie relation: E (meV) = h²/(2 m_n λ²) = 81.8042 / λ² (Å).
            Spallation (time-of-flight) sources are out of scope.

        Raises
        ------
        ValueError
            If set to a value not in ``("xray", "neutron")``.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.source_type
        'xray'
        >>> g.source_type = "neutron"
        >>> g.source_type
        'neutron'
        """
        return self._source_type

    @source_type.setter
    def source_type(self, value: str) -> None:
        from .radiation import SOURCE_TYPES

        if value not in SOURCE_TYPES:
            raise ValueError(
                f"source_type must be one of {SOURCE_TYPES!r}; got {value!r}."
            )
        self._source_type = value

    @property
    def energy_units(self) -> str:
        """
        Units of the :attr:`energy` property.

        ``"keV"`` for X-ray geometries; ``"meV"`` for neutron geometries.

        Returns
        -------
        str
        """
        return "keV" if self._source_type == "xray" else "meV"

    @property
    def energy(self) -> float | None:
        """
        Radiation energy in units of :attr:`energy_units`, derived from wavelength.

        - X-ray: E (keV) = hc / λ = 12.39842 / λ (Å)
        - Neutron (reactor): E (meV) = h²/(2 m_n λ²) = 81.8042 / λ² (Å)

        Returns ``None`` when ``wavelength`` is not set.  Setting ``energy``
        updates ``wavelength`` using the inverse formula for the current
        :attr:`source_type`.  The value must be in the units given by
        :attr:`energy_units`.

        Raises
        ------
        ValueError
            If set to a non-positive value.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> round(g.energy, 4)    # keV for xray
        8.0479
        >>> g.source_type = "neutron"
        >>> round(g.energy, 4)    # meV for neutron
        34.4634
        >>> g.energy = 25.0       # sets wavelength via de Broglie
        >>> round(g.wavelength, 4)
        1.809
        """
        if self._wavelength is None:
            return None
        if self._source_type == "xray":
            from .radiation import wavelength_to_energy

            return wavelength_to_energy(self._wavelength)
        else:
            from .radiation import neutron_wavelength_to_energy

            return neutron_wavelength_to_energy(self._wavelength)

    @energy.setter
    def energy(self, value: float | None) -> None:
        if value is None:
            self._wavelength = None
            return
        if self._source_type == "xray":
            from .radiation import energy_to_wavelength

            self._wavelength = energy_to_wavelength(float(value))
        else:
            from .radiation import neutron_energy_to_wavelength

            self._wavelength = neutron_energy_to_wavelength(float(value))

    @property
    def wavenumber(self) -> float | None:
        """
        Wave number k in Å⁻¹, derived from wavelength.

        k (Å⁻¹) = 2π / λ (Å)

        Returns ``None`` when ``wavelength`` is not set.  Setting
        ``wavenumber`` updates ``wavelength`` accordingly.

        Raises
        ------
        ValueError
            If set to a non-positive value.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> round(g.wavenumber, 4)
        4.0774
        >>> g.wavenumber = 4.0
        >>> round(g.wavelength, 4)
        1.5708
        """
        if self._wavelength is None:
            return None
        from .radiation import wavelength_to_wavenumber

        return wavelength_to_wavenumber(self._wavelength)

    @wavenumber.setter
    def wavenumber(self, value: float | None) -> None:
        if value is None:
            self._wavelength = None
            return
        from .radiation import wavenumber_to_wavelength

        self._wavelength = wavenumber_to_wavelength(float(value))

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

    @property
    def azimuthal_reference(self) -> tuple[float, float, float] | None:
        """
        Azimuthal reference vector as Miller indices (h, k, l), or ``None``.

        The azimuthal reference defines the direction in reciprocal space
        used to set the azimuthal angle ψ = 0.  ψ is zero when this vector
        lies in the scattering plane (the plane containing the incident beam
        and the scattering vector Q).

        Setting to ``None`` clears the reference (ψ becomes undefined).

        Setting to a tuple ``(h, k, l)`` stores it as a tuple of three
        floats.  Any three-element sequence is accepted.

        Parameters
        ----------
        value : tuple of float or None
            Miller indices of the azimuthal reference direction, e.g.
            ``(0, 0, 1)`` for the surface normal of a (001)-cut crystal.

        Raises
        ------
        ValueError
            If the value is not ``None`` and cannot be interpreted as a
            length-3 sequence of numbers, or if the vector is the zero
            vector.

        Examples
        --------
        >>> g = psic()
        >>> g.azimuthal_reference = (0, 0, 1)   # c-axis surface normal
        >>> g.azimuthal_reference
        (0.0, 0.0, 1.0)
        >>> g.azimuthal_reference = None          # clear reference
        """
        return self._azimuthal_reference

    @azimuthal_reference.setter
    def azimuthal_reference(self, value: tuple[float, float, float] | None) -> None:
        if value is None:
            self._azimuthal_reference = None
            return
        try:
            h, k, l = (float(x) for x in value)  # noqa: E741
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "azimuthal_reference must be a length-3 sequence of numbers "
                f"or None; got {value!r}."
            ) from exc
        if h == 0.0 and k == 0.0 and l == 0.0:
            raise ValueError(
                "azimuthal_reference must be a non-zero vector; (0, 0, 0) is not allowed."
            )
        self._azimuthal_reference = (h, k, l)

    # ------------------------------------------------------------------
    # Diffraction modes
    # ------------------------------------------------------------------

    @property
    def free_dof_after_bragg(self) -> int:
        """
        Number of free degrees of freedom after satisfying the Bragg
        condition.

        Equals ``N − 3`` where N is the total number of real motor axes
        (sample stages + detector stages).  The value determines how many
        constraints a :class:`~mode.ConstraintSet` must contain for this
        geometry.

        Returns
        -------
        int

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> ahd.fourcv().free_dof_after_bragg
        1
        >>> ahd.psic().free_dof_after_bragg
        3
        """
        return len(self.sample_stages) + len(self.detector_stages) - 3

    @property
    def modes(self) -> ModeDict:
        """
        The collection of named diffraction modes for this geometry.

        A ``ModeDict`` mapping mode name (str) to
        :class:`~mode.ConstraintSet` instance.  Empty when no modes
        have been declared.

        Returns
        -------
        ModeDict
        """
        return self._modes

    @property
    def mode_name(self) -> str | None:
        """
        Name of the currently active diffraction mode, or ``None``.

        Setting this property switches the active mode.  The name must
        be a key of :attr:`modes`.

        Returns
        -------
        str or None

        Raises
        ------
        ValueError
            If the supplied name is not in :attr:`modes`.
        """
        return self._mode_name

    @mode_name.setter
    def mode_name(self, name: str | None) -> None:
        if name is not None and name not in self._modes:
            available = sorted(self._modes.keys())
            raise ValueError(
                f"Mode {name!r} is not available in this geometry. "
                f"Available modes: {available}."
            )
        self._mode_name = name

    @property
    def mode(self) -> ConstraintSet | None:
        """
        The currently active :class:`~mode.ConstraintSet`, or ``None``.

        Equivalent to ``self.modes[self.mode_name]`` when a mode is
        active; ``None`` when :attr:`mode_name` is ``None``.

        Returns
        -------
        ConstraintSet or None
        """
        if self._mode_name is None:
            return None
        return self._modes[self._mode_name]

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
    # Samples
    # ------------------------------------------------------------------

    @property
    def _samples(self) -> SampleDict:  # pragma: no cover
        """
        The guarded sample dict.  Read-only property prevents replacement
        of the entire dict (``g._samples = something`` raises AttributeError).
        """
        return self.__samples

    @property
    def samples(self) -> SampleDict:
        """
        The guarded SampleDict of all samples.

        Supports read-only iteration (``for name in g.samples``, ``g.samples[name]``).
        Mutation goes through ``add_sample()`` / ``remove_sample()`` to keep
        invariants intact.
        """
        return self.__samples

    @property
    def sample(self) -> Sample:
        """The currently active Sample."""
        return self.__samples[self._active_ref[0]]

    @sample.setter
    def sample(self, value: str | Sample) -> None:
        """
        Set the active sample by name or Sample object.

        Raises
        ------
        KeyError
            If the name is not in the sample dict.
        """
        name = value.name if isinstance(value, Sample) else value
        if name not in self.__samples:
            raise KeyError(f"No sample named {name!r}. Add it first with add_sample().")
        self._active_ref[0] = name

    def add_sample(
        self,
        name: str,
        lattice=None,
    ) -> Sample:
        """
        Add a new sample and return it.

        Parameters
        ----------
        name : str
            Unique label for this sample.
        lattice : Lattice or None
            Crystal lattice.  If None, defaults to cubic a = 1 Å.

        Returns
        -------
        Sample

        Raises
        ------
        ValueError
            If a sample with that name already exists.
        """
        if name in self.__samples:
            raise ValueError(f"A sample named {name!r} already exists.")
        from .lattice import Lattice as _Lattice

        lat = lattice if lattice is not None else _Lattice(a=1.0)
        s = Sample(
            name=name,
            lattice=lat,
            reflections=ReflectionList(
                geometry_name=self.name,
                valid_stages=set(self._stages),
            ),
            parent=self,
        )
        self.__samples._data[name] = s  # bypass guard: add is always safe
        return s

    def remove_sample(self, name: str) -> None:
        """
        Remove a sample from the dict.

        The active sample cannot be removed; select a different sample
        first.  Clears ``sample.parent`` on the removed sample so that any
        external reference to it knows it is detached.  Delegates to
        SampleDict which enforces the active-sample invariant.

        Raises
        ------
        KeyError
            If the sample does not exist.
        ValueError
            If the sample is the currently active sample.
        """
        s = self.__samples[name]  # raises KeyError if not found
        del self.__samples[name]  # raises ValueError if active
        s.parent = None  # break the back-reference

    # ------------------------------------------------------------------
    # Reflections — convenience wrapper targeting the active sample
    # ------------------------------------------------------------------

    def add_reflection(
        self,
        name: str,
        hkl: tuple[float, float, float],
        angles: dict[str, float],
        wavelength: float | None = None,
    ) -> Reflection:
        """
        Add a named orienting reflection to the active sample.

        Delegates to ``self.sample.reflections.add()``.  If ``wavelength``
        is None the geometry's current wavelength is inherited.

        Parameters
        ----------
        name : str
            Unique label (e.g. ``"r1"``, ``"Si_111"``).
        hkl : tuple of float
            Miller indices (h, k, l).
        angles : dict[str, float]
            Motor angles in degrees keyed by stage name.
        wavelength : float or None
            Wavelength in Å.  If None, inherits ``self.wavelength``.

        Returns
        -------
        Reflection
        """
        wl = wavelength if wavelength is not None else self._wavelength
        return self.sample.reflections.add(
            name=name,
            hkl=hkl,
            angles=angles,
            wavelength=wl,
        )

    def forward(
        self,
        h: float,
        k: float,
        l: float,  # noqa: E741
    ) -> list[dict[str, float]]:
        r"""
        Compute all valid motor-angle solutions for the reflection (h, k, l).

        This is the *forward* diffraction calculation: given a reciprocal-
        lattice point, compute the motor angles that satisfy the Bragg
        condition under the active diffraction mode.

        Algorithm (delegated to :mod:`.forward`):

        1. Compute ``Q_phi = UB @ (h, k, l)`` — the target scattering vector
           in the phi frame.
        2. Apply Bragg's law to find the detector angle (2θ).
        3. Apply mode constraints (frozen angles, bisecting condition) to
           determine the remaining free stage angles.
        4. Return all valid solutions, filtered by stage limits.

        Parameters
        ----------
        h, k, l : float
            Miller indices of the target reflection.

        Returns
        -------
        list of dict[str, float]
            Each element is a complete set of motor angles (all stage names
            as keys, values in degrees) that satisfies the Bragg condition
            under the active mode.  May be empty if no valid solution exists
            within limits.

        Raises
        ------
        ValueError
            If ``self.wavelength`` is None.
        ValueError
            If ``self.sample.UB`` is None.
        ValueError
            If (h, k, l) == (0, 0, 0).
        ValueError
            If the requested \|Q\| exceeds the Ewald sphere.
        NotImplementedError
            If no active mode is set or the mode type is not supported.

        See Also
        --------
        inverse : Convert motor angles to (h, k, l).

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> ahd.ub_identity(g.sample)
        >>> solutions = g.forward(1, 0, 0)
        >>> len(solutions) > 0
        True
        >>> sorted(solutions[0].keys())
        ['chi', 'omega', 'phi', 'ttheta']
        """
        from .forward import compute_forward

        return compute_forward(self, h, k, l)

    def inverse(self, angles: dict[str, float]) -> tuple[float, float, float]:
        """
        Convert a set of motor angles to Miller indices (h, k, l).

        This is the *inverse* diffraction calculation: given where the
        diffractometer motors are pointing, compute the reciprocal-lattice
        point being measured.

        Algorithm:

        1. Compute ``Q_phi`` from the motor angles via
           :func:`angles_to_phi_vector`.
        2. Solve ``UB @ hkl = Q_phi`` for ``hkl``, i.e.
           ``hkl = UB⁻¹ @ Q_phi``.

        Parameters
        ----------
        angles : dict[str, float]
            Motor angles in degrees, keyed by stage name.  All stages
            present in the geometry may be supplied; stages not supplied
            keep their current ``angle`` attribute.

        Returns
        -------
        hkl : tuple of float
            Miller indices ``(h, k, l)`` corresponding to the supplied
            motor angles.

        Raises
        ------
        KeyError
            If a supplied stage name does not exist in the geometry.
        ValueError
            If ``self.wavelength`` is None.
        ValueError
            If the active sample has no UB matrix set
            (``sample.UB is None``).
        numpy.linalg.LinAlgError
            If the UB matrix is singular and cannot be inverted.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.psic()
        >>> g.wavelength = 1.5406
        >>> ahd.ub_identity(g.sample)          # set UB = B (cubic a=1)
        >>> hkl = g.inverse(
        ...     {"mu": 0, "eta": 20.97, "chi": 90, "phi": 0,
        ...      "nu": 0, "delta": 41.94}
        ... )
        """
        from .orientation import angles_to_phi_vector

        sample = self.sample
        if sample.UB is None:
            raise ValueError(
                f"Sample {sample.name!r} has no UB matrix. "
                "Set one with ub_identity(), ub_from_one_reflection(), or "
                "assign sample.UB directly."
            )

        Q_phi = angles_to_phi_vector(self, **angles)
        hkl = np.linalg.solve(sample.UB, Q_phi)
        return (float(hkl[0]), float(hkl[1]), float(hkl[2]))

    def psi(self, angles: dict[str, float] | None = None) -> float:
        """
        Compute the azimuthal angle ψ (psi) for the given motor angles.

        ψ is the rotation of the azimuthal reference vector about the
        scattering vector Q, measured from the zero-ψ position.  ψ = 0
        when the reference vector lies in the scattering plane (the plane
        containing the incident beam direction and Q).

        Algorithm (You 1999, eqs. 10-11):

        1. Compute ``Q_phi`` from the motor angles via
           :func:`angles_to_phi_vector`.
        2. Express the reference direction in the phi frame:
           ``n_phi = UB @ n_hkl``  (maps reciprocal-lattice direction to
           phi frame in the same units as Q).
        3. Project both ``n_phi`` and the incident-beam direction
           ``y_hat`` (the ``'longitudinal'`` basis vector) onto the plane
           perpendicular to ``Q_phi``.
        4. ψ = ``atan2(Q_hat · (y_perp × n_perp), y_perp · n_perp)``
           where all vectors are unit vectors in the perpendicular plane.

        Parameters
        ----------
        angles : dict[str, float] or None
            Motor angles in degrees, keyed by stage name.  If ``None``
            (default), the current stage angles are used.

        Returns
        -------
        psi : float
            Azimuthal angle in degrees, in the range (-180, 180].

        Raises
        ------
        ValueError
            If ``self.wavelength`` is None.
        ValueError
            If ``self.sample.UB`` is None.
        ValueError
            If ``self.azimuthal_reference`` is None (no reference set).
        ValueError
            If the reference vector is parallel to Q (ψ is undefined).

        Notes
        -----
        The zero-ψ reference direction is the component of the incident
        beam (``'longitudinal'`` basis vector) perpendicular to Q.  This
        matches the SPEC and You (1999) convention: ψ = 0 when the
        reference vector lies in the scattering plane.

        The incident beam direction is taken from the geometry's basis
        dict under the key ``'longitudinal'``.  For You (1999) geometries
        this is YHAT; for Busing & Levy geometries it is also YHAT.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> g.azimuthal_reference = (0, 0, 1)
        >>> # ... set UB, move motors ...
        >>> psi = g.psi()

        References
        ----------
        You, J. Appl. Cryst. 32, 614-623 (1999), eqs. 10-11.
        """
        import math

        if self.wavelength is None:
            raise ValueError("psi() requires geometry.wavelength to be set.")
        if self.sample.UB is None:
            raise ValueError(
                "psi() requires a UB matrix on the active sample. "
                "Call ub_identity() or ub_from_two_reflections_bl1967() first."
            )
        if self.azimuthal_reference is None:
            raise ValueError(
                "psi() requires azimuthal_reference to be set on the geometry. "
                "Set geometry.azimuthal_reference = (h, k, l) first."
            )

        if angles is None:
            angles = {s.name: s.angle for s in self._stages.values()}

        from .orientation import angles_to_phi_vector

        Q_phi = angles_to_phi_vector(self, **angles)
        Q_mag = np.linalg.norm(Q_phi)
        if Q_mag < 1e-14:
            raise ValueError("psi() is undefined when Q = 0 (all motors at zero).")
        Q_hat = Q_phi / Q_mag

        # Reference direction in phi frame
        n_hkl = np.asarray(self.azimuthal_reference, dtype=float)
        n_phi = self.sample.UB @ n_hkl
        n_mag = np.linalg.norm(n_phi)
        if n_mag < 1e-14:
            raise ValueError(
                "Azimuthal reference vector maps to zero in the phi frame. "
                "Check that the UB matrix is non-singular."
            )
        n_hat = n_phi / n_mag

        # Incident beam direction from basis ('longitudinal' key)
        y_hat = np.asarray(
            self.basis.get("longitudinal", np.array([0.0, 1.0, 0.0])),
            dtype=float,
        )
        y_hat = y_hat / np.linalg.norm(y_hat)

        # Project n and y onto the plane perpendicular to Q
        n_perp = n_hat - np.dot(n_hat, Q_hat) * Q_hat
        y_perp = y_hat - np.dot(y_hat, Q_hat) * Q_hat

        n_perp_mag = np.linalg.norm(n_perp)
        y_perp_mag = np.linalg.norm(y_perp)

        if n_perp_mag < 1e-10:
            raise ValueError(
                "The azimuthal reference vector is parallel to Q; "
                "psi is undefined at this reflection."
            )
        if y_perp_mag < 1e-10:
            raise ValueError(
                "The incident beam direction is parallel to Q; "
                "psi is undefined at this motor position."
            )

        n_perp_hat = n_perp / n_perp_mag
        y_perp_hat = y_perp / y_perp_mag

        cos_psi = float(np.clip(np.dot(y_perp_hat, n_perp_hat), -1.0, 1.0))
        sin_psi = float(np.dot(Q_hat, np.cross(y_perp_hat, n_perp_hat)))
        return math.degrees(math.atan2(sin_psi, cos_psi))

    # ------------------------------------------------------------------
    # Surface geometry: surface_normal property
    # ------------------------------------------------------------------

    @property
    def surface_normal(self) -> tuple[float, float, float] | None:
        """
        Surface normal direction as Miller indices (h, k, l), or ``None``.

        The surface normal defines the direction perpendicular to the sample
        surface in reciprocal space.  It is used by the surface geometry
        calculations (:meth:`alpha_i`, :meth:`alpha_f`, :meth:`q_components`,
        :meth:`is_specular`, :meth:`is_evanescent`).

        When ``None``, the surface calculations fall back to
        :attr:`azimuthal_reference` if that is set.

        Setting this to ``None`` clears the surface normal.
        Setting to a non-zero (h, k, l) 3-tuple stores it.

        Raises
        ------
        ValueError
            If the supplied value is not a length-3 sequence of numbers,
            or is the zero vector (0, 0, 0).

        Examples
        --------
        >>> g = psic()
        >>> g.surface_normal = (0, 0, 1)   # c-axis surface normal
        >>> g.surface_normal
        (0.0, 0.0, 1.0)
        >>> g.surface_normal = None        # clear
        """
        return self._surface_normal

    @surface_normal.setter
    def surface_normal(self, value: tuple[float, float, float] | None) -> None:
        if value is None:
            self._surface_normal = None
            return
        try:
            h, k, l = (float(x) for x in value)  # noqa: E741
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "surface_normal must be a length-3 sequence of numbers "
                f"or None; got {value!r}."
            ) from exc
        if h == 0.0 and k == 0.0 and l == 0.0:
            raise ValueError(
                "surface_normal must be a non-zero vector; (0, 0, 0) is not allowed."
            )
        self._surface_normal = (h, k, l)

    # ------------------------------------------------------------------
    # Detector geometry parameters
    # ------------------------------------------------------------------

    @property
    def detector_distance(self) -> float | None:
        """
        Sample-to-detector distance in millimetres, or ``None`` if not set.

        Relevant primarily for area detectors, where the distance is needed
        to convert pixel positions to scattering angles (and thence to Q).
        Point detectors do not require this parameter.

        Must be strictly positive when set.

        Raises
        ------
        ValueError
            If the supplied value is not positive.

        Examples
        --------
        >>> g = psic()
        >>> g.detector_distance = 500.0   # 500 mm
        >>> g.detector_distance
        500.0
        >>> g.detector_distance = None    # clear
        """
        return self._detector_distance

    @detector_distance.setter
    def detector_distance(self, value: float | None) -> None:
        if value is None:
            self._detector_distance = None
            return
        value = float(value)
        if value <= 0.0:
            raise ValueError(f"detector_distance must be positive (got {value!r} mm).")
        self._detector_distance = value

    @property
    def detector_tilt(self) -> float | None:
        """
        Detector tilt angle in degrees, or ``None`` if not set.

        A non-zero tilt indicates that the detector plane is not perfectly
        perpendicular to the diffracted beam at the nominal detector position.
        The tilt is a rotation of the detector face about an axis parallel
        to the detector plane (a pitch/roll correction for non-ideal
        alignment).

        Any real number is accepted (positive or negative; modulo 360° is
        not applied).

        Examples
        --------
        >>> g = psic()
        >>> g.detector_tilt = 0.3    # 0.3° tilt
        >>> g.detector_tilt
        0.3
        >>> g.detector_tilt = None   # clear
        """
        return self._detector_tilt

    @detector_tilt.setter
    def detector_tilt(self, value: float | None) -> None:
        if value is None:
            self._detector_tilt = None
            return
        self._detector_tilt = float(value)

    @property
    def detector_offset(self) -> tuple[float, float] | None:
        """
        In-plane detector offset (dx, dy) in millimetres, or ``None`` if not set.

        Describes the displacement of the beam centre from the detector
        centre in the detector plane.  ``dx`` is the horizontal offset and
        ``dy`` is the vertical offset, both in millimetres.  A zero offset
        means the beam strikes the geometric centre of the detector.

        Must be a two-element sequence of real numbers.

        Raises
        ------
        ValueError
            If the supplied value is not a length-2 sequence of numbers.

        Examples
        --------
        >>> g = psic()
        >>> g.detector_offset = (2.5, -1.0)   # 2.5 mm right, 1.0 mm down
        >>> g.detector_offset
        (2.5, -1.0)
        >>> g.detector_offset = None           # clear
        """
        return self._detector_offset

    @detector_offset.setter
    def detector_offset(self, value: tuple[float, float] | None) -> None:
        if value is None:
            self._detector_offset = None
            return
        try:
            dx, dy = (float(x) for x in value)
        except (TypeError, ValueError) as exc:
            raise ValueError(
                "detector_offset must be a length-2 sequence of numbers (dx, dy) "
                f"in mm, or None; got {value!r}."
            ) from exc
        self._detector_offset = (dx, dy)

    # ------------------------------------------------------------------
    # Surface geometry: incidence / emergence / Q-decomposition methods
    # ------------------------------------------------------------------

    def alpha_i(self, angles: dict[str, float] | None = None) -> float:
        """
        Angle of incidence αᵢ (degrees).

        αᵢ is the angle between the incoming beam and the sample surface.
        Requires ``wavelength``, ``sample.UB``, and ``surface_normal``
        (or ``azimuthal_reference``) to be set.

        Parameters
        ----------
        angles : dict[str, float] or None
            Motor angles.  If ``None``, current stage angles are used.

        Returns
        -------
        float
            αᵢ in degrees, in [0°, 90°].

        See Also
        --------
        alpha_f : Angle of emergence.

        Examples
        --------
        >>> g = zaxis()
        >>> g.wavelength = 1.5406
        >>> g.surface_normal = (0, 0, 1)
        >>> ub_identity(g.sample)
        >>> g.alpha_i({"alpha": 5.0, "Z": 0.0, "delta": 20.0, "gamma": 0.0})
        5.0
        """
        from .surface import alpha_i as _alpha_i

        return _alpha_i(self, angles)

    def alpha_f(self, angles: dict[str, float] | None = None) -> float:
        """
        Angle of emergence αf (degrees).

        αf is the angle between the diffracted beam and the sample surface.
        Requires ``wavelength``, ``sample.UB``, and ``surface_normal``
        (or ``azimuthal_reference``) to be set.

        Parameters
        ----------
        angles : dict[str, float] or None
            Motor angles.  If ``None``, current stage angles are used.

        Returns
        -------
        float
            αf in degrees, in [0°, 90°].

        See Also
        --------
        alpha_i : Angle of incidence.
        """
        from .surface import alpha_f as _alpha_f

        return _alpha_f(self, angles)

    def q_components(self, angles: dict[str, float] | None = None) -> dict[str, float]:
        """
        Decompose Q into components parallel and perpendicular to the surface.

        Parameters
        ----------
        angles : dict[str, float] or None
            Motor angles.  If ``None``, current stage angles are used.

        Returns
        -------
        dict with keys ``"Q_perp"``, ``"Q_par"``, ``"Q_perp_signed"``,
        ``"Q_total"`` — all in Å⁻¹.

        See Also
        --------
        alpha_i, alpha_f : Incidence and emergence angles.
        """
        from .surface import q_components as _q_components

        return _q_components(self, angles)

    def is_specular(
        self,
        angles: dict[str, float] | None = None,
        atol: float = 0.01,
    ) -> bool:
        """
        Return True when αᵢ ≈ αf within ``atol`` degrees.

        Parameters
        ----------
        angles : dict[str, float] or None
        atol : float
            Tolerance in degrees (default 0.01°).

        Returns
        -------
        bool
        """
        from .surface import is_specular as _is_specular

        return _is_specular(self, angles, atol)

    def is_evanescent(
        self,
        angles: dict[str, float] | None = None,
        critical_angle_deg: float | None = None,
    ) -> bool:
        """
        Return True when αᵢ < ``critical_angle_deg`` (evanescent regime).

        Parameters
        ----------
        angles : dict[str, float] or None
        critical_angle_deg : float
            Critical angle for total external reflection in degrees.
            Must be supplied (material-dependent; not stored on geometry).

        Returns
        -------
        bool

        Raises
        ------
        ValueError
            If ``critical_angle_deg`` is None.
        """
        from .surface import is_evanescent as _is_evanescent

        return _is_evanescent(self, angles, critical_angle_deg)

    # ------------------------------------------------------------------
    # Diffractometer inclination
    # ------------------------------------------------------------------

    @property
    def inclination_matrix(self) -> np.ndarray:
        """
        3×3 rotation matrix describing the diffractometer inclination
        relative to the incident beam direction.

        When the entire diffractometer is mounted at a non-zero angle with
        respect to the incident beam (e.g. tilted for grazing-incidence
        geometry or non-standard instrument mounting), this matrix encodes
        that orientation.

        The matrix **R** is applied as follows: the effective incident-beam
        direction in the diffractometer coordinate frame is
        ``R.T @ ŷ`` rather than ``ŷ``.  A zero inclination (identity matrix)
        reproduces the standard geometry exactly.

        The matrix must be a proper rotation (orthonormal, det = +1).

        Default: ``numpy.eye(3)`` (no inclination).

        Raises
        ------
        ValueError
            If the supplied array is not shape (3, 3), not orthonormal
            within 1e-8, or has det ≠ +1.

        See Also
        --------
        set_inclination : Set the inclination from an axis and angle.

        Examples
        --------
        >>> import numpy as np
        >>> g = psic()
        >>> g.inclination_matrix  # default: identity
        array([[1., 0., 0.],
               [0., 1., 0.],
               [0., 0., 1.]])
        >>> g.set_inclination(axis=[0, 0, 1], angle_deg=2.0)  # 2° about z
        >>> np.round(g.inclination_matrix, 4)
        array([[ 0.9994,  0.0349,  0.    ],
               [-0.0349,  0.9994,  0.    ],
               [ 0.    ,  0.    ,  1.    ]])
        """
        return self._inclination_matrix

    @inclination_matrix.setter
    def inclination_matrix(self, value: np.ndarray) -> None:
        arr = np.asarray(value, dtype=float)
        if arr.shape != (3, 3):
            raise ValueError(
                f"inclination_matrix must be a (3, 3) array; got shape {arr.shape}."
            )
        # Check orthonormality
        err = np.max(np.abs(arr.T @ arr - np.eye(3)))
        if err > 1e-8:
            raise ValueError(
                f"inclination_matrix must be orthonormal (R.T @ R = I); "
                f"max deviation = {err:.2e}."
            )
        det = float(np.linalg.det(arr))
        if abs(det - 1.0) > 1e-8:
            raise ValueError(
                f"inclination_matrix must be a proper rotation (det = +1); "
                f"got det = {det:.8f}."
            )
        self._inclination_matrix = arr

    def set_inclination(
        self,
        axis: np.ndarray | list[float],
        angle_deg: float,
    ) -> None:
        """
        Set the inclination matrix from a rotation axis and angle.

        Convenience wrapper around :attr:`inclination_matrix` that uses the
        :func:`~.rotation.rotation_matrix` function to build R from an axis
        vector and an angle in degrees.

        A zero angle resets the inclination to the identity matrix.

        Parameters
        ----------
        axis : array-like, shape (3,)
            Rotation axis (need not be a unit vector; it is normalised
            internally).
        angle_deg : float
            Rotation angle in degrees.

        Raises
        ------
        ValueError
            If ``axis`` is the zero vector.

        Examples
        --------
        >>> g = psic()
        >>> g.set_inclination(axis=[0, 0, 1], angle_deg=2.0)
        >>> g.set_inclination(axis=[1, 0, 0], angle_deg=0.0)  # reset
        """
        from .rotation import rotation_matrix

        ax = np.asarray(axis, dtype=float)
        ax_norm = np.linalg.norm(ax)
        if ax_norm < 1e-14:
            raise ValueError("set_inclination(): axis must be a non-zero vector.")
        ax = ax / ax_norm
        self.inclination_matrix = rotation_matrix(ax, float(angle_deg))

    def wh(self, print: bool = True) -> str:
        """
        Terse one-screen status string (the SPEC ``wh`` command).

        Shows the current reciprocal-space position (H K L), azimuthal
        angle ψ, wavelength, and a motor-angle table with SPEC-style
        column names.  Graceful fallback text is shown for any field
        that cannot be computed (e.g. no UB matrix or no wavelength).

        Parameters
        ----------
        print : bool, optional
            If ``True`` (default), the string is printed to stdout before
            being returned.  Pass ``print=False`` to suppress output and
            capture the string instead (useful in tests and logging).

        Returns
        -------
        str
            Multi-line status string.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> g.wh()                          # prints and returns
        H K L = not available
        Psi = not available
        Lambda = 1.5406
        <BLANKLINE>
          TwoTheta     Theta       Chi       Phi
             0.000     0.000     0.000     0.000
        >>> out = g.wh(print=False)         # capture without printing
        >>> "Lambda" in out
        True

        References
        ----------
        Align4Pete.log — ``wh`` command outputs, 7-ID-C fourc session,
        Dec 2020.
        """
        _print = builtins.print
        lines: list[str] = []

        # HKL position
        hkl_str = "not available"
        try:
            current_angles = {s.name: s.angle for s in self._stages.values()}
            hkl = self.inverse(current_angles)
            hkl_str = "  {:g}  {:g}  {:g}".format(*[self._clean_zero(v) for v in hkl])
        except Exception as exc:  # noqa: BLE001
            logger.debug("wh: could not compute HKL: %s", exc)
        lines.append(f"H K L = {hkl_str}")

        # Azimuthal angle ψ
        psi_str = "not available"
        try:
            psi_str = f"{self.psi():.4g}"
        except Exception as exc:  # noqa: BLE001
            logger.debug("wh: could not compute psi: %s", exc)
        lines.append(f"Psi = {psi_str}")

        # Wavelength
        lam = self.wavelength
        lines.append(f"Lambda = {lam:g}" if lam is not None else "Lambda = not set")

        # Motor angle table
        lines.append("")
        stage_names = list(self._stages.keys())
        lines.append("".join(f"{self._spec_motor_name(n):>10s}" for n in stage_names))
        lines.append("".join(f"{self._stages[n].angle:>10.3f}" for n in stage_names))

        result = "\n".join(lines)
        if print:
            _print(result)
        return result

    def pa(self, print: bool = True) -> str:
        """
        Verbose parameter listing (the SPEC ``pa`` command).

        Shows the geometry name, the two designated orienting reflections
        (if set), lattice constants in real and reciprocal space, the
        azimuthal reference vector, and the wavelength.

        Parameters
        ----------
        print : bool, optional
            If ``True`` (default), the string is printed to stdout before
            being returned.  Pass ``print=False`` to suppress output and
            capture the string instead (useful in tests and logging).

        Returns
        -------
        str
            Multi-line parameter string.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.pa()                          # prints and returns
        Geometry: fourcv
        ...
        >>> out = g.pa(print=False)         # capture without printing
        >>> "Geometry" in out
        True

        References
        ----------
        Align4Pete.log — ``pa`` command outputs, 7-ID-C fourc session,
        Dec 2020.
        """
        _print = builtins.print
        lines: list[str] = []
        lines.append(f"Geometry: {self.name}")
        lines.append("")

        # Orienting reflections
        sample = self.sample
        ors = sample.reflections.orienting_reflections

        def _refl_block(label: str, refl) -> list[str]:
            if refl is None:
                return [f"  {label}: not set"]
            lam_r = refl.wavelength if refl.wavelength is not None else self.wavelength
            lam_str = f"{lam_r:g}" if lam_r is not None else "not set"
            ang_str = "  ".join(f"{v:g}" for v in refl.angles.values())
            ang_keys = "  ".join(refl.angles.keys())
            h, k, l = refl.hkl  # noqa: E741
            return [
                f"  {label} (at lambda {lam_str}):",
                f"    {ang_keys} = {ang_str}",
                f"    H K L = {self._clean_zero(h):g}"
                f"  {self._clean_zero(k):g}  {self._clean_zero(l):g}",
            ]

        or1 = ors[0] if len(ors) >= 1 else None
        or2 = ors[1] if len(ors) >= 2 else None
        lines.extend(_refl_block("Primary Reflection", or1))
        lines.append("")
        lines.extend(_refl_block("Secondary Reflection", or2))
        lines.append("")

        # Lattice constants
        lat = sample.lattice
        lines.append("  Lattice Constants (lengths / angles):")
        lines.append(
            f"      real space = {lat.a:g} {lat.b:g} {lat.c:g}"
            f" / {lat.alpha:g} {lat.beta:g} {lat.gamma:g}"
        )

        rvecs = lat.reciprocal_lattice_vectors
        r1, r2, r3 = (
            np.asarray(rvecs[0]),
            np.asarray(rvecs[1]),
            np.asarray(rvecs[2]),
        )
        a_star = float(np.linalg.norm(r1))
        b_star = float(np.linalg.norm(r2))
        c_star = float(np.linalg.norm(r3))

        def _ang(u, v):
            import math as _math

            cos_a = float(
                np.clip(
                    np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1.0, 1.0
                )
            )
            return _math.degrees(_math.acos(cos_a))

        lines.append(
            f"    reciprocal space = {a_star:.4g} {b_star:.4g} {c_star:.4g}"
            f" / {_ang(r2, r3):.4g} {_ang(r1, r3):.4g} {_ang(r1, r2):.4g}"
        )
        lines.append("")

        # Azimuthal reference
        az_ref = self.azimuthal_reference
        if az_ref is not None:
            h, k, l = az_ref  # noqa: E741
            lines.append(
                f"  Azimuthal Reference:  H K L = "
                f"{self._clean_zero(h):g}  {self._clean_zero(k):g}"
                f"  {self._clean_zero(l):g}"
            )
        else:
            lines.append("  Azimuthal Reference:  not set")
        lines.append("")

        # Wavelength
        lam = self.wavelength
        lines.append(f"  Lambda = {lam:g}" if lam is not None else "  Lambda = not set")

        result = "\n".join(lines)
        if print:
            _print(result)
        return result

    # ------------------------------------------------------------------
    # Private status helpers (used by wh and pa methods)
    # ------------------------------------------------------------------

    @staticmethod
    def _clean_zero(v: float, atol: float = 1e-10) -> float:
        """Return 0.0 if ``abs(v) < atol``, else ``v``.

        Prevents ``-0.0`` or floating-point noise like ``1e-16`` from
        appearing in status output.
        """
        return 0.0 if abs(v) < atol else v

    @staticmethod
    def _spec_motor_name(name: str) -> str:
        """Map internal stage names to SPEC-style column headers.

        Unknown names pass through unchanged.
        """
        _MAP = {
            "ttheta": "TwoTheta",
            "omega": "Theta",
            "chi": "Chi",
            "phi": "Phi",
            "mu": "Mu",
            "eta": "Eta",
            "nu": "Nu",
            "delta": "Delta",
        }
        return _MAP.get(name, name)

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
            wl_str = (
                f"{fmt(self._wavelength)} Å"
                f"  (E = {fmt(self.energy)} {self.energy_units})"
            )
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

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    @staticmethod
    def _mode_to_dict(mode: "ConstraintSet") -> dict:
        """
        Serialise a :class:`~mode.ConstraintSet` to a JSON-serialisable dict.
        """
        return mode.to_dict()

    @staticmethod
    def _mode_from_dict(d: dict) -> "ConstraintSet":
        """
        Reconstruct a :class:`~mode.ConstraintSet` from a dict produced by
        :meth:`_mode_to_dict`.

        Raises
        ------
        ValueError
            If the dict does not have ``"type": "ConstraintSet"``.
        """
        from .mode import ConstraintSet

        type_name = d.get("type", "")
        if type_name != "ConstraintSet":
            raise ValueError(
                f"_mode_from_dict: expected type 'ConstraintSet'; got {type_name!r}. "
                "The old BisectingMode / FixedAngleMode serialisation format is no "
                "longer supported.  Re-save the geometry with the current version."
            )
        return ConstraintSet.from_dict(d)

    def to_dict(self) -> dict:
        """
        Export the complete diffractometer configuration as a
        JSON-serialisable ``dict``.

        The returned dict captures all geometry state: stage definitions,
        current motor angles, all samples (including their lattices,
        reflections, and U/UB matrices), the active sample name,
        wavelength, azimuthal reference, and identifying metadata.

        Returns
        -------
        dict
            JSON-serialisable mapping with the following top-level keys:

            ``"_meta"``
                Sub-dict with ``"software"`` (``"ad_hoc_diffractometer"``),
                ``"version"`` (package version string), and
                ``"created"`` (ISO-8601 UTC timestamp).
            ``"name"``
                Geometry name (str).
            ``"description"``
                Free-text description (str).
            ``"wavelength"``
                Wavelength in Å (float or None).
            ``"kappa_alpha_deg"``
                Kappa tilt angle in degrees (float or None).
            ``"azimuthal_reference"``
                Miller indices [h, k, l] (list of 3 float, or None).
            ``"basis"``
                Dict mapping physical direction names to [x, y, z] lists.
            ``"stages"``
                List of stage dicts (name, axis, role, parent, angle,
                limits).
            ``"active_sample"``
                Name of the currently active sample (str).
            ``"samples"``
                Dict mapping sample names to sample dicts.

        Examples
        --------
        >>> import json, ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> d = g.to_dict()
        >>> assert json.dumps(d)          # must be JSON-serialisable
        >>> d["name"]
        'fourcv'
        >>> d["wavelength"]
        1.5406
        """
        import datetime

        try:
            from importlib.metadata import version as _version

            _pkg_version = _version("ad_hoc_diffractometer")
        except Exception as exc:  # noqa: BLE001
            logger.debug("Could not read package version: %s", exc)
            _pkg_version = "unknown"

        stages = [s.to_dict() for s in self._stages.values()]

        samples = {
            name: sample.to_dict() for name, sample in self.__samples._data.items()
        }

        # Serialise modes: each mode records its type name and constructor fields.
        modes_dict = {}
        for mode_name, mode_obj in self._modes.items():
            modes_dict[mode_name] = self._mode_to_dict(mode_obj)

        return {
            "_meta": {
                "software": "ad_hoc_diffractometer",
                "version": _pkg_version,
                "created": datetime.datetime.now(datetime.timezone.utc).isoformat(),
            },
            "name": self.name,
            "description": self.description,
            "wavelength": self._wavelength,
            "kappa_alpha_deg": self._kappa_alpha_deg,
            "azimuthal_reference": (
                list(self._azimuthal_reference)
                if self._azimuthal_reference is not None
                else None
            ),
            "surface_normal": (
                list(self._surface_normal) if self._surface_normal is not None else None
            ),
            "detector_distance": self._detector_distance,
            "detector_tilt": self._detector_tilt,
            "detector_offset": (
                list(self._detector_offset)
                if self._detector_offset is not None
                else None
            ),
            "inclination_matrix": self._inclination_matrix.tolist(),
            "basis": {k: [float(x) for x in v] for k, v in self.basis.items()},
            "stages": stages,
            "active_sample": self._active_ref[0],
            "samples": samples,
            "modes": modes_dict,
            "mode_name": self._mode_name,
            "cut_points": dict(self.cut_points),
        }

    @classmethod
    def from_dict(cls, d: dict) -> "AdHocDiffractometer":
        """
        Reconstruct an :class:`AdHocDiffractometer` from a dict produced by
        :meth:`to_dict`.

        Parameters
        ----------
        d : dict
            Dict as returned by :meth:`to_dict`.  The ``"_meta"`` key is
            read but not validated; version mismatches do not raise.

        Returns
        -------
        AdHocDiffractometer
            A fully-configured geometry with all stages, samples, and
            settings restored.

        Notes
        -----
        Stages are reconstructed from the stored axis, role, parent, limits,
        and angle values.  The stage parent graph is rebuilt from the stored
        ``"parent"`` field (a name string or None).

        Samples are restored via :meth:`Sample.from_dict`.  The active
        sample is restored by name.

        The factory-function identity (e.g. ``psic``, ``fourcv``) is **not**
        stored; ``from_dict`` always returns a plain
        :class:`AdHocDiffractometer` instance.  If you need the factory
        function, look it up via ``get_geometry(d["name"])``.

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> g2 = ahd.AdHocDiffractometer.from_dict(g.to_dict())
        >>> g2.name
        'fourcv'
        >>> g2.wavelength
        1.5406
        """
        import numpy as np

        from .sample import Sample
        from .stage import Stage

        basis = {k: np.array(v, dtype=float) for k, v in d["basis"].items()}

        # Build Stage objects via Stage.from_dict()
        stage_objects: dict[str, Stage] = {}
        for sd in d["stages"]:
            s = Stage.from_dict(sd)
            stage_objects[sd["name"]] = s

        stages_list = list(stage_objects.values())

        # Restore modes
        raw_modes = d.get("modes", {})
        restored_modes: dict[str, ConstraintSet] = {}
        for mname, mdict in raw_modes.items():
            restored_modes[mname] = cls._mode_from_dict(mdict)

        geom = cls(
            name=d["name"],
            stages=stages_list,
            basis=basis,
            description=d.get("description", ""),
            wavelength=d.get("wavelength"),
            kappa_alpha_deg=d.get("kappa_alpha_deg"),
            azimuthal_reference=d.get("azimuthal_reference"),
            modes=restored_modes if restored_modes else None,
            default_mode=d.get("mode_name"),
            cut_points=d.get("cut_points"),
        )

        # Restore surface_normal
        sn = d.get("surface_normal")
        if sn is not None:
            geom.surface_normal = tuple(sn)

        # Restore detector geometry parameters
        if d.get("detector_distance") is not None:
            geom.detector_distance = d["detector_distance"]
        if d.get("detector_tilt") is not None:
            geom.detector_tilt = d["detector_tilt"]
        do = d.get("detector_offset")
        if do is not None:
            geom.detector_offset = tuple(do)

        # Restore inclination_matrix (always present in to_dict output; the
        # 'is None' branch guards against hand-crafted or legacy dicts)
        inc = d.get("inclination_matrix")
        if inc is not None:  # pragma: no branch
            geom.inclination_matrix = np.array(inc, dtype=float)

        # Restore samples.
        #
        # Order of operations is important to avoid the SampleDict guard
        # blocking deletion of the default "test" sample:
        #
        #   1. Write all restored samples into _data directly (bypass the
        #      guard — we are building a consistent state from scratch).
        #   2. Update the active-sample pointer to the restored active name
        #      BEFORE removing any samples, so the guard sees the correct
        #      active name when we delete stale entries.
        #   3. Remove any samples that were not in the exported dict (e.g.
        #      the default "test" sample created by __init__ that was not
        #      present in the original geometry at export time).
        active_name = d.get("active_sample", "test")
        saved_samples: dict[str, Sample] = {}
        for sample_name, sd in d.get("samples", {}).items():
            sample = Sample.from_dict(sd, parent=geom)
            sample.reflections.geometry_name = geom.name
            saved_samples[sample_name] = sample

        # Write restored samples directly into _data (bypass guard)
        geom.samples._data.update(saved_samples)

        # Switch active pointer before removing stale samples (step 2)
        if active_name in geom.samples._data:  # pragma: no branch
            geom._active_ref[0] = active_name

        # Remove samples that were not in the exported dict (step 3).
        # The active-sample guard now protects the correct (restored) active
        # sample, so any stale default sample that is not the active one
        # can be safely removed.
        stale = [n for n in list(geom.samples._data) if n not in saved_samples]
        for n in stale:
            if n != geom._active_ref[0]:  # pragma: no branch
                del geom.samples._data[n]

        return geom


# ---------------------------------------------------------------------------
# Public top-level convenience functions (SPEC-familiar status commands)
# ---------------------------------------------------------------------------


def wh(geometry: AdHocDiffractometer, print: bool = True) -> str:
    """
    Terse one-screen status for a diffractometer (the SPEC ``wh`` command).

    Shows the current reciprocal-space position (H K L), azimuthal angle ψ,
    wavelength, and a motor-angle table with SPEC-style column names.
    Graceful fallback text is shown for any field that cannot be computed
    (e.g. no UB matrix or no wavelength set).

    Equivalent to ``geometry.wh(print=print)``.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer to report on.
    print : bool, optional
        If ``True`` (default), the string is printed to stdout before
        being returned.  Pass ``print=False`` to suppress output and
        capture the string instead.

    Returns
    -------
    str
        Multi-line status string.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.wavelength = 1.5406
    >>> out = ahd.wh(g, print=False)
    >>> "Lambda" in out
    True

    See Also
    --------
    :func:`pa` : Verbose parameter listing.
    :meth:`AdHocDiffractometer.wh` : Equivalent method on the geometry object.
    """
    return geometry.wh(print=print)


def pa(geometry: AdHocDiffractometer, print: bool = True) -> str:
    """
    Verbose parameter listing for a diffractometer (the SPEC ``pa`` command).

    Shows the geometry name, the two designated orienting reflections (if
    set), lattice constants in real and reciprocal space, the azimuthal
    reference vector, and the wavelength.

    Equivalent to ``geometry.pa(print=print)``.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer to report on.
    print : bool, optional
        If ``True`` (default), the string is printed to stdout before
        being returned.  Pass ``print=False`` to suppress output and
        capture the string instead.

    Returns
    -------
    str
        Multi-line parameter string.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> out = ahd.pa(g, print=False)
    >>> "Geometry" in out
    True

    See Also
    --------
    :func:`wh` : Terse one-screen status.
    :meth:`AdHocDiffractometer.pa` : Equivalent method on the geometry object.
    """
    return geometry.pa(print=print)
