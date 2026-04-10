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
from .sample import _DEFAULT_LATTICE
from .sample import _DEFAULT_SAMPLE_NAME
from .sample import Sample
from .sample import SampleDict
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
    azimuthal_reference : tuple of float or None, optional
        Azimuthal reference direction as Miller indices (h, k, l).  Used
        by :meth:`psi` to compute the azimuthal angle ψ.  ``None`` (default)
        means no reference is set.  Must be a non-zero vector.

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
    azimuthal_reference : tuple of float or None
        Azimuthal reference vector (h, k, l), or None if not set.
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
    ):
        self.name = name
        self.description = description
        self.basis = basis if basis is not None else dict(self.DEFAULT_BASIS)
        self.wavelength = wavelength  # validated via property setter
        self.kappa_alpha_deg = kappa_alpha_deg
        self.azimuthal_reference = azimuthal_reference  # validated via property setter

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
    def _samples(self) -> SampleDict:
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

    @property
    def wh(self) -> str:
        """
        Terse one-screen status string (the SPEC ``wh`` command).

        Returns the same string as ``status.wh(self)``: current H K L
        position, azimuthal angle ψ, wavelength, and a motor-angle table
        with SPEC-style column names.

        Read as a property so the call site reads naturally::

            print(g.wh)

        Returns
        -------
        str
            Multi-line status string, ready to ``print()``.

        See Also
        --------
        ad_hoc_diffractometer.status.wh : Module-level function (thin wrapper).

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> g.wavelength = 1.5406
        >>> print(g.wh)
        H K L = not available
        Psi = not available
        Lambda = 1.5406
        <BLANKLINE>
          TwoTheta     Theta       Chi       Phi
             0.000     0.000     0.000     0.000
        """
        from .status import wh as _wh

        return _wh(self)

    @property
    def pa(self) -> str:
        """
        Verbose parameter listing (the SPEC ``pa`` command).

        Returns the same string as ``status.pa(self)``: geometry name,
        primary and secondary orienting reflections, lattice constants
        in real and reciprocal space, azimuthal reference, and wavelength.

        Read as a property so the call site reads naturally::

            print(g.pa)

        Returns
        -------
        str
            Multi-line parameter string, ready to ``print()``.

        See Also
        --------
        ad_hoc_diffractometer.status.pa : Module-level function (thin wrapper).

        Examples
        --------
        >>> import ad_hoc_diffractometer as ahd
        >>> g = ahd.fourcv()
        >>> print(g.pa)                     # doctest: +ELLIPSIS
        Geometry: fourcv
        ...
        """
        from .status import pa as _pa

        return _pa(self)

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
