# Copyright (c) 2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-ANL-Open-Source-License
"""
lattice.py — crystallographic lattice calculations.

Provides:
  - Lattice class: describes a crystal lattice with lazy-computed properties
    for Cartesian lattice vectors, reciprocal lattice vectors, and the B matrix.
  - lattice_vectors(): standalone function (kept for backward compatibility)
  - reciprocal_vectors(): standalone function
  - b_matrix(): standalone function

Based on:
  Busing & Levy, Acta Cryst. 22, 457-464 (1967)
  I16 Diffractometer Rotation Matrix document
  Lecture-2-Reciprocal-lattice-notes
"""

from __future__ import annotations

import logging
import math

import numpy as np

from .display import fmt

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Crystal system definitions
# ---------------------------------------------------------------------------

_SYSTEM_FREE_PARAMS: dict[str, tuple[str, ...]] = {
    "cubic": ("a",),
    "tetragonal": ("a", "c"),
    "orthorhombic": ("a", "b", "c"),
    "hexagonal": ("a", "c"),
    "trigonal": ("a", "alpha"),
    "monoclinic": ("a", "b", "c", "beta"),
    "triclinic": ("a", "b", "c", "alpha", "beta", "gamma"),
}
"""Parameters that are independent (not symmetry-constrained) for each system.

These are the parameters reported by ``__str__``.  Ordered from most
constrained (cubic) to most general (triclinic).
"""

CRYSTAL_SYSTEMS = list(_SYSTEM_FREE_PARAMS)
"""The seven crystal systems, derived from :data:`_SYSTEM_FREE_PARAMS`.

Membership and order are always consistent with the free-parameter table.
"""

# Sentinel: marks a parameter as "not supplied by the caller"
_UNSET = object()

# ---------------------------------------------------------------------------
# Validation
# ---------------------------------------------------------------------------

_ANGLE_MIN = 1.0  # degrees — prevents degenerate cells
_ANGLE_MAX = 179.0  # degrees


def _validate_params(a, b, c, alpha, beta, gamma) -> None:
    """
    Validate a complete set of lattice parameters.

    Raises
    ------
    ValueError
        If any length is not strictly positive, any angle is outside
        (_ANGLE_MIN, _ANGLE_MAX), or the parameters would produce a
        degenerate (zero-volume) unit cell.
    """
    for name, val in (("a", a), ("b", b), ("c", c)):
        if val <= 0.0:
            raise ValueError(
                f"Lattice parameter {name!r} must be strictly positive; got {val}."
            )

    for name, val in (("alpha", alpha), ("beta", beta), ("gamma", gamma)):
        if not (_ANGLE_MIN < val < _ANGLE_MAX):
            raise ValueError(
                f"Lattice angle {name!r} must be in ({_ANGLE_MIN}, {_ANGLE_MAX}) degrees; "
                f"got {val}."
            )

    # Check unit cell volume is non-zero (triangle inequality on angles)
    alpha_r = np.deg2rad(alpha)
    beta_r = np.deg2rad(beta)
    gamma_r = np.deg2rad(gamma)
    volume_factor = (
        1
        - np.cos(alpha_r) ** 2
        - np.cos(beta_r) ** 2
        - np.cos(gamma_r) ** 2
        + 2 * np.cos(alpha_r) * np.cos(beta_r) * np.cos(gamma_r)
    )
    if volume_factor <= 0.0:
        raise ValueError(
            f"Lattice angles alpha={alpha}, beta={beta}, gamma={gamma} "
            f"produce a degenerate (zero-volume) unit cell."
        )


# ---------------------------------------------------------------------------
# Crystal system deduction
# ---------------------------------------------------------------------------


def _deduce_system_and_params(
    a, b, c, alpha, beta, gamma
) -> tuple[str, float, float, float, float, float, float]:
    """
    Deduce the crystal system from the supplied parameters and fill in
    symmetry-required defaults for any that were not supplied.

    Caller supplies values or _UNSET for each parameter.  'a' is always
    required.  All other parameters are optional; their defaults depend on
    the deduced crystal system.

    Deduction rules (most constrained first):

    Cubic:         a only
    Trigonal:      a + alpha (alpha != 90)
    Hexagonal:     a + c + gamma=120 (gamma must be explicit)
    Tetragonal:    a + c (no angles)
    Orthorhombic:  a + b + c (no angles)
    Monoclinic:    a + b + c + beta
    Triclinic:     all six

    When all six parameters are supplied, the system is classified by
    examining the symmetry of the values.

    Parameters
    ----------
    a, b, c, alpha, beta, gamma
        Each is either a float (supplied) or _UNSET (not supplied).

    Returns
    -------
    system : str
    a, b, c, alpha, beta, gamma : float (complete, validated)

    Raises
    ------
    ValueError
        If 'a' is not supplied, if the supplied combination is ambiguous
        or inconsistent, or if validated parameters are out of range.
    """
    has_a = a is not _UNSET
    has_b = b is not _UNSET
    has_c = c is not _UNSET
    has_alpha = alpha is not _UNSET
    has_beta = beta is not _UNSET
    has_gamma = gamma is not _UNSET

    if not has_a:
        raise ValueError("Lattice parameter 'a' must always be provided.")

    # --- All six supplied: classify by symmetry of values ---
    if has_a and has_b and has_c and has_alpha and has_beta and has_gamma:
        a, b, c = float(a), float(b), float(c)
        alpha, beta, gamma = float(alpha), float(beta), float(gamma)
        _validate_params(a, b, c, alpha, beta, gamma)

        eps = 1e-9
        alpha_90 = abs(alpha - 90.0) < eps
        beta_90 = abs(beta - 90.0) < eps
        gamma_90 = abs(gamma - 90.0) < eps
        gamma_120 = abs(gamma - 120.0) < eps
        all_angles_90 = alpha_90 and beta_90 and gamma_90
        ab_equal = abs(a - b) < eps
        abc_equal = ab_equal and abs(b - c) < eps
        angles_equal = abs(alpha - beta) < eps and abs(beta - gamma) < eps

        if abc_equal and all_angles_90:
            system = "cubic"
        elif ab_equal and alpha_90 and beta_90 and gamma_120:
            system = "hexagonal"
        elif abc_equal and angles_equal and not all_angles_90:
            system = "trigonal"
        elif ab_equal and all_angles_90:
            system = "tetragonal"
        elif all_angles_90:
            system = "orthorhombic"
        elif alpha_90 and gamma_90 and not beta_90:
            system = "monoclinic"
        else:
            system = "triclinic"

        return system, a, b, c, alpha, beta, gamma

    # --- Trigonal (rhombohedral): a + alpha, no b/c/beta/gamma ---
    if (
        has_a
        and has_alpha
        and not has_b
        and not has_c
        and not has_beta
        and not has_gamma
    ):
        a, alpha = float(a), float(alpha)
        if abs(alpha - 90.0) < 1e-9:
            raise ValueError(
                "For trigonal (rhombohedral setting), alpha must differ from 90 degrees."
            )
        b, c, beta, gamma = a, a, alpha, alpha
        _validate_params(a, b, c, alpha, beta, gamma)
        return "trigonal", a, b, c, alpha, beta, gamma

    # --- Hexagonal: a + c + gamma=120 explicitly ---
    if has_a and has_c and has_gamma and not has_b and not has_alpha and not has_beta:
        a, c, gamma = float(a), float(c), float(gamma)
        if abs(gamma - 120.0) > 1e-9:
            raise ValueError(f"For hexagonal, gamma must be 120 degrees; got {gamma}.")
        b, alpha, beta = a, 90.0, 90.0
        _validate_params(a, b, c, alpha, beta, gamma)
        return "hexagonal", a, b, c, alpha, beta, gamma

    # --- Tetragonal: a + c, no angles ---
    if (
        has_a
        and has_c
        and not has_b
        and not has_alpha
        and not has_beta
        and not has_gamma
    ):
        a, c = float(a), float(c)
        b, alpha, beta, gamma = a, 90.0, 90.0, 90.0
        _validate_params(a, b, c, alpha, beta, gamma)
        return "tetragonal", a, b, c, alpha, beta, gamma

    # --- Cubic: a only ---
    if (
        has_a
        and not has_b
        and not has_c
        and not has_alpha
        and not has_beta
        and not has_gamma
    ):
        a = float(a)
        b, c, alpha, beta, gamma = a, a, 90.0, 90.0, 90.0
        _validate_params(a, b, c, alpha, beta, gamma)
        return "cubic", a, b, c, alpha, beta, gamma

    # --- Orthorhombic: a + b + c, no angles ---
    if has_a and has_b and has_c and not has_alpha and not has_beta and not has_gamma:
        a, b, c = float(a), float(b), float(c)
        alpha, beta, gamma = 90.0, 90.0, 90.0
        _validate_params(a, b, c, alpha, beta, gamma)
        return "orthorhombic", a, b, c, alpha, beta, gamma

    # --- Monoclinic: a + b + c + beta ---
    if has_a and has_b and has_c and has_beta and not has_alpha and not has_gamma:
        a, b, c, beta = float(a), float(b), float(c), float(beta)
        alpha, gamma = 90.0, 90.0
        _validate_params(a, b, c, alpha, beta, gamma)
        return "monoclinic", a, b, c, alpha, beta, gamma

    raise ValueError(
        "Cannot deduce crystal system from the supplied parameters. "
        "Provide at minimum:\n"
        "  a                      -> cubic\n"
        "  a, c                   -> tetragonal\n"
        "  a, c, gamma=120        -> hexagonal\n"
        "  a, alpha (!=90)        -> trigonal\n"
        "  a, b, c                -> orthorhombic\n"
        "  a, b, c, beta          -> monoclinic\n"
        "  a, b, c, alpha, beta, gamma -> triclinic"
    )


# ---------------------------------------------------------------------------
# Lattice class
# ---------------------------------------------------------------------------


class Lattice:
    """
    A crystal lattice described by up to six parameters (a, b, c, α, β, γ).

    The crystal system is deduced automatically from the **minimum set** of
    parameters supplied.  Missing parameters are filled in according to the
    symmetry constraints of the deduced system.  The default (no arguments)
    is a cubic lattice with a = 1 Å.

    .. list-table:: Minimum parameters by crystal system
       :header-rows: 1
       :widths: 20 50

       * - System
         - Minimum parameters required
       * - Cubic
         - ``a``
       * - Tetragonal
         - ``a``, ``c``
       * - Hexagonal
         - ``a``, ``c``, ``gamma=120`` *(must be explicit)*
       * - Trigonal
         - ``a``, ``alpha`` *(alpha ≠ 90)*
       * - Orthorhombic
         - ``a``, ``b``, ``c``
       * - Monoclinic
         - ``a``, ``b``, ``c``, ``beta``
       * - Triclinic
         - ``a``, ``b``, ``c``, ``alpha``, ``beta``, ``gamma``

    All parameters are validated on construction and on every change:
      - Lengths (a, b, c) must be strictly positive.
      - Angles (alpha, beta, gamma) must be in (1, 179) degrees.
      - The combination must produce a non-degenerate unit cell.

    Cartesian lattice vectors, reciprocal lattice vectors, and the B matrix
    are computed lazily: calculated on first access and cached until any
    lattice parameter changes.

    Parameters
    ----------
    a : float, optional
        Lattice parameter a in Angstroms.  Default 1.0.
    b : float, optional
        Lattice parameter b in Angstroms.
    c : float, optional
        Lattice parameter c in Angstroms.
    alpha : float, optional
        Angle between b and c axes in degrees.
    beta : float, optional
        Angle between a and c axes in degrees.
    gamma : float, optional
        Angle between a and b axes in degrees.

    Properties (read/write)
    -----------------------
    a, b, c : float — lattice parameters in Angstroms
    alpha, beta, gamma : float — lattice angles in degrees
    precision : int or None — decimal places for display; None uses the
        package-level default from display.get_precision()

    Properties (read-only, lazy)
    ----------------------------
    system : str — deduced crystal system
    cartesian_lattice_vectors : tuple(a1, a2, a3) — direct lattice vectors (Å)
    reciprocal_lattice_vectors : tuple(b1, b2, b3) — reciprocal vectors (Å⁻¹, with 2π)
    B : ndarray (3,3) — B matrix (Å⁻¹, no 2π)

    Examples
    --------
    >>> Lattice()                                     # cubic, a=1 Å
    >>> Lattice(a=5.0)                                # cubic, a=5 Å
    >>> Lattice(a=3.0, c=6.0)                         # tetragonal
    >>> Lattice(a=4.785, c=12.991, gamma=120.0)       # hexagonal
    >>> Lattice(a=5.0, alpha=60.0)                    # trigonal
    >>> Lattice(a=2.0, b=3.0, c=4.0)                 # orthorhombic
    >>> Lattice(a=5.0, b=6.0, c=7.0, beta=110.0)     # monoclinic
    >>> Lattice(a=3, b=4, c=5, alpha=80, beta=85, gamma=95)  # triclinic
    >>> lat = Lattice(a=4.785, c=12.991, gamma=120.0)
    >>> lat.precision = 3
    >>> print(lat)   # shows 3 decimal places
    """

    def __init__(
        self,
        a: float = 1.0,
        b: float = _UNSET,
        c: float = _UNSET,
        alpha: float = _UNSET,
        beta: float = _UNSET,
        gamma: float = _UNSET,
        precision: int | None = None,
    ):
        system, a_, b_, c_, alpha_, beta_, gamma_ = _deduce_system_and_params(
            a, b, c, alpha, beta, gamma
        )
        self._system = system
        self._a = a_
        self._b = b_
        self._c = c_
        self._alpha = alpha_
        self._beta = beta_
        self._gamma = gamma_

        self._cartesian_lattice_vectors: tuple | None = None
        self._reciprocal_lattice_vectors: tuple | None = None
        self._B: np.ndarray | None = None

        # Display precision: None means use package-level default
        self.precision: int | None = precision

    # ------------------------------------------------------------------
    # Cache management
    # ------------------------------------------------------------------

    def _invalidate(self) -> None:
        """Mark all cached computations as stale."""
        self._cartesian_lattice_vectors = None
        self._reciprocal_lattice_vectors = None
        self._B = None

    def _set_params(self, **kwargs) -> None:
        """
        Update one or more lattice parameters, re-deduce the crystal system,
        and invalidate the cache.  Parameters are always passed to the
        deduction function in the canonical order a, b, c, alpha, beta, gamma.

        All six current values are passed through, with kwargs overriding.
        This allows re-deduction across crystal systems when a parameter
        change alters the symmetry (e.g. setting b != a on a cubic lattice
        causes re-deduction to orthorhombic).

        Individual parameter validation (positive lengths, angle range) is
        performed first, before deduction, so that errors are reported with
        the correct parameter name rather than a generic deduction error.
        """
        # Validate only the parameters being changed, before touching state
        for name in ("a", "b", "c"):
            if name in kwargs and kwargs[name] <= 0.0:
                raise ValueError(
                    f"Lattice parameter {name!r} must be strictly positive; "
                    f"got {kwargs[name]}."
                )
        for name in ("alpha", "beta", "gamma"):
            if name in kwargs:
                val = kwargs[name]
                if not (_ANGLE_MIN < val < _ANGLE_MAX):
                    raise ValueError(
                        f"Lattice angle {name!r} must be in "
                        f"({_ANGLE_MIN}, {_ANGLE_MAX}) degrees; got {val}."
                    )

        # Determine which parameters to pass for re-deduction.
        #
        # If ALL changed parameters are free (independent) for the current
        # system, we pass only the free parameters (using _UNSET for
        # constrained ones) so that symmetry constraints are re-applied
        # cleanly.  This keeps a cubic lattice cubic when only 'a' changes.
        #
        # If ANY changed parameter is NOT free for the current system (e.g.
        # setting 'b' on a cubic lattice), the user is intentionally changing
        # the symmetry.  In that case we pass all six current values plus the
        # override, allowing the all-six classification path to run.
        free = set(_SYSTEM_FREE_PARAMS[self._system])
        changing = set(kwargs)
        all_free = changing.issubset(free)

        if all_free:
            # Pass only free params; constrained ones become _UNSET
            a = kwargs.get("a", self._a) if "a" in free else _UNSET
            b = kwargs.get("b", self._b) if "b" in free else _UNSET
            c = kwargs.get("c", self._c) if "c" in free else _UNSET
            alpha = kwargs.get("alpha", self._alpha) if "alpha" in free else _UNSET
            beta = kwargs.get("beta", self._beta) if "beta" in free else _UNSET
            gamma = kwargs.get("gamma", self._gamma) if "gamma" in free else _UNSET
        else:
            # Pass all six so the all-six classification path can run
            a = kwargs.get("a", self._a)
            b = kwargs.get("b", self._b)
            c = kwargs.get("c", self._c)
            alpha = kwargs.get("alpha", self._alpha)
            beta = kwargs.get("beta", self._beta)
            gamma = kwargs.get("gamma", self._gamma)

        system, a_, b_, c_, alpha_, beta_, gamma_ = _deduce_system_and_params(
            a, b, c, alpha, beta, gamma
        )
        self._system = system
        self._a = a_
        self._b = b_
        self._c = c_
        self._alpha = alpha_
        self._beta = beta_
        self._gamma = gamma_
        self._invalidate()

    # ------------------------------------------------------------------
    # Read/write parameter properties
    # ------------------------------------------------------------------

    @property
    def system(self) -> str:
        """Deduced crystal system name."""
        return self._system

    @property
    def a(self) -> float:
        """Lattice parameter a (Å)."""
        return self._a

    @a.setter
    def a(self, value: float) -> None:
        self._set_params(a=float(value))

    @property
    def b(self) -> float:
        """Lattice parameter b (Å)."""
        return self._b

    @b.setter
    def b(self, value: float) -> None:
        self._set_params(b=float(value))

    @property
    def c(self) -> float:
        """Lattice parameter c (Å)."""
        return self._c

    @c.setter
    def c(self, value: float) -> None:
        self._set_params(c=float(value))

    @property
    def alpha(self) -> float:
        """Lattice angle alpha (degrees)."""
        return self._alpha

    @alpha.setter
    def alpha(self, value: float) -> None:
        self._set_params(alpha=float(value))

    @property
    def beta(self) -> float:
        """Lattice angle beta (degrees)."""
        return self._beta

    @beta.setter
    def beta(self, value: float) -> None:
        self._set_params(beta=float(value))

    @property
    def gamma(self) -> float:
        """Lattice angle gamma (degrees)."""
        return self._gamma

    @gamma.setter
    def gamma(self, value: float) -> None:
        self._set_params(gamma=float(value))

    # ------------------------------------------------------------------
    # Lazy computed properties
    # ------------------------------------------------------------------

    @property
    def cartesian_lattice_vectors(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Direct-lattice vectors (a1, a2, a3) in Angstroms, expressed in
        the BL1967 crystal-Cartesian frame.

        The frame is defined by the BL1967 B-matrix convention: the
        reciprocal vector b1\\* is along crystal-x̂, b2\\* lies in the
        x̂-ŷ plane, b3\\* completes the right-handed triple.  The
        direct vectors are then determined by the duality
        ``aᵢ · bⱼ\\* = 2π δᵢⱼ``.  Their magnitudes equal the lattice
        parameters (``|a1| = a``, ``|a2| = b``, ``|a3| = c``) and the
        angles between them recover (α, β, γ).

        For orthogonal cells (cubic, tetragonal, orthorhombic) this
        reduces to ``a1 = a x̂``, ``a2 = b ŷ``, ``a3 = c ẑ``.  For
        non-orthogonal cells the direct vectors are rotated relative
        to that simple placement so that the reciprocal vectors land
        in the BL1967 frame.  See :func:`lattice_vectors` for the
        construction.

        Cached until a parameter changes.
        """
        if self._cartesian_lattice_vectors is None:
            self._cartesian_lattice_vectors = lattice_vectors(
                self._a,
                self._b,
                self._c,
                self._alpha,
                self._beta,
                self._gamma,
            )
        return self._cartesian_lattice_vectors

    @property
    def reciprocal_lattice_vectors(self) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
        """
        Reciprocal-lattice vectors (b1, b2, b3) in Å⁻¹ (with 2π factor)
        expressed in the BL1967 crystal-Cartesian frame.

        The columns of :attr:`B`.  Satisfies ``bᵢ · aⱼ = 2π δᵢⱼ`` with
        the direct vectors returned by
        :attr:`cartesian_lattice_vectors`.  Cached until a parameter
        changes.
        """
        if self._reciprocal_lattice_vectors is None:
            self._reciprocal_lattice_vectors = reciprocal_vectors(
                self._a,
                self._b,
                self._c,
                self._alpha,
                self._beta,
                self._gamma,
            )
        return self._reciprocal_lattice_vectors

    @property
    def B(self) -> np.ndarray:
        """
        B matrix in Å⁻¹ (with 2π factor, BL1967 / SPEC / hkl_soleil
        convention).

        Transforms Miller indices ``h`` = (h, k, l) to the scattering
        vector expressed in the BL1967 crystal-Cartesian frame::

            Q_c = B @ h

        ``|B @ h| = 2π / d_hkl``.  The orthogonality condition is
        ``bᵢ · aⱼ = 2π δᵢⱼ``.

        Numerical equivalence
        ---------------------
        Returns exactly the same B matrix (to machine precision) as
        SPEC, hkl_soleil (default ``HKL_TAU = 2π``), and diffcalc-core,
        for every crystal system.  See :func:`b_matrix_bl1967` for the
        closed-form derivation and the cross-solver references.

        Other packages (FullProf, CrysFML, parts of CCP4, the hkl
        library compiled with ``HKL_TAU = 1``) use a no-2π variant
        where ``bᵢ · aⱼ = δᵢⱼ`` and ``|B @ h| = 1/d_hkl``; only the
        2π variant is exposed by this package.

        Cached until a lattice parameter changes.
        """
        if self._B is None:
            self._B = b_matrix_bl1967(
                self._a,
                self._b,
                self._c,
                self._alpha,
                self._beta,
                self._gamma,
            )
        return self._B

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------

    def __str__(self) -> str:
        """
        Human-readable representation reporting only the free (independent)
        parameters for the deduced crystal system.

        Display precision is controlled by self.precision (per-instance) or
        the package-level default from display.get_precision() if
        self.precision is None.
        """
        free = _SYSTEM_FREE_PARAMS[self._system]
        param_strs = []
        for p in free:
            val = getattr(self, p)
            unit = "Å" if p in ("a", "b", "c") else "°"
            param_strs.append(f"{p}={fmt(val, self.precision)} {unit}")
        return f"Lattice({self._system}: {', '.join(param_strs)})"

    def __eq__(self, other: object, atol: float | None = None) -> bool:
        """
        True if all six lattice parameters agree within tolerance.

        Parameters
        ----------
        other : object
        atol : float or None
            Absolute tolerance.  If None, derived from the current display
            precision via ``display.precision_atol()``.

        Returns
        -------
        bool
        """
        if not isinstance(other, Lattice):
            return NotImplemented
        from .display import allclose

        return allclose(
            [self._a, self._b, self._c, self._alpha, self._beta, self._gamma],
            [other._a, other._b, other._c, other._alpha, other._beta, other._gamma],
            atol=atol,
        )

    def __repr__(self) -> str:
        return (
            f"Lattice(system={self._system!r}, "
            f"a={self._a}, b={self._b}, c={self._c}, "
            f"alpha={self._alpha}, beta={self._beta}, gamma={self._gamma})"
        )

    # ------------------------------------------------------------------
    # Serialisation
    # ------------------------------------------------------------------

    def to_dict(self) -> dict:
        """
        Return a JSON-serialisable ``dict`` representing this lattice.

        All six parameters are always stored (even constrained ones), so
        the dict is unambiguous regardless of crystal system.

        Returns
        -------
        dict
            Keys: ``"a"``, ``"b"``, ``"c"``, ``"alpha"``, ``"beta"``,
            ``"gamma"`` (all ``float``).

        Examples
        --------
        >>> Lattice(a=4.785, c=12.991, gamma=120).to_dict()
        {'a': 4.785, 'b': 4.785, 'c': 12.991, 'alpha': 90.0, 'beta': 90.0, 'gamma': 120.0}
        """
        return {
            "a": self._a,
            "b": self._b,
            "c": self._c,
            "alpha": self._alpha,
            "beta": self._beta,
            "gamma": self._gamma,
        }

    @classmethod
    def from_dict(cls, d: dict) -> Lattice:
        """
        Reconstruct a ``Lattice`` from a dict produced by :meth:`to_dict`.

        Parameters
        ----------
        d : dict
            Must contain keys ``"a"``, ``"b"``, ``"c"``, ``"alpha"``,
            ``"beta"``, ``"gamma"``.

        Returns
        -------
        Lattice

        Examples
        --------
        >>> Lattice.from_dict({'a': 5.431, 'b': 5.431, 'c': 5.431,
        ...                    'alpha': 90.0, 'beta': 90.0, 'gamma': 90.0})
        Lattice(system='cubic', a=5.431, b=5.431, c=5.431, ...)
        """
        return cls(
            a=d["a"],
            b=d["b"],
            c=d["c"],
            alpha=d["alpha"],
            beta=d["beta"],
            gamma=d["gamma"],
        )


# ---------------------------------------------------------------------------
# Standalone functions (kept for backward compatibility and direct use)
# ---------------------------------------------------------------------------


def b_matrix_bl1967(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> np.ndarray:
    r"""
    Compute the B matrix directly from lattice parameters using the
    Busing & Levy 1967 closed-form expression (eq. 3).

    The columns of B are the reciprocal-lattice vectors b1\*, b2\*, b3\*
    expressed in the BL1967 crystal-Cartesian frame: b1\* along
    crystal-x, b2\* in the crystal-x–crystal-y plane (with positive
    y component), b3\* completing the right-handed orthonormal triple.
    With this choice B is upper triangular::

        B[0, 0] = a*               B[0, 1] = b* cos γ*    B[0, 2] = c* cos β*
        B[1, 0] = 0                B[1, 1] = b* sin γ*    B[1, 2] = -c* sin β* cos α
        B[2, 0] = 0                B[2, 1] = 0            B[2, 2] = 2π / c

    where (a\*, b\*, c\*) are the reciprocal-lattice vector magnitudes
    (including the 2π factor) and (α\*, β\*, γ\*) are the reciprocal-cell
    angles, related to the direct-cell parameters by the standard
    crystallographic identities.

    Numerical equivalence with external solvers
    -------------------------------------------
    Returns exactly the same B matrix (to machine precision) as

    * **SPEC** psic and friends,
    * **hkl_soleil** (the libhkl library) with its default
      ``HKL_TAU = 2π`` compile-time setting,
    * **diffcalc-core** (``Crystal._set_reciprocal_cell``).

    This holds for every crystal system, including non-orthogonal
    cells (hexagonal, trigonal, monoclinic, triclinic).  For
    orthogonal cells (cubic, tetragonal, orthorhombic), the
    direct-lattice vectors `a`, `b`, `c` are parallel to the reciprocal
    vectors a\*, b\*, c\* respectively, and the BL1967 frame coincides
    with the "direct-lattice-a-along-x" frame that earlier versions of
    this package used.  For non-orthogonal cells the two frames
    differ by a rotation in the crystal-Cartesian space; both encode
    the same physical crystal, but only the BL1967 frame is the one
    SPEC and hkl_soleil cross-validate against.

    Parameters
    ----------
    a, b, c : float
        Direct-lattice parameters in Angstroms.
    alpha, beta, gamma : float
        Direct-lattice angles in degrees.
        alpha = angle between b and c axes
        beta  = angle between a and c axes
        gamma = angle between a and b axes

    Returns
    -------
    B : numpy.ndarray, shape (3, 3)
        B matrix in inverse Angstroms (with 2π factor, BL1967 / SPEC /
        hkl_soleil convention).  Upper triangular by construction.

    References
    ----------
    * W. R. Busing & H. A. Levy, *Acta Cryst.* **22**, 457-464 (1967),
      eq. 3.
    * libhkl source: ``hkl/hkl-lattice.c`` (``HKL_TAU = 2π``).
    * diffcalc-core source:
      ``src/diffcalc/ub/crystal.py`` (``_set_reciprocal_cell``).
    """
    al = math.radians(alpha)
    be = math.radians(beta)
    ga = math.radians(gamma)

    # Direct-cell volume scale: D = sqrt(1 - cos²α - cos²β - cos²γ + 2 cos α cos β cos γ).
    D = math.sqrt(
        1.0
        - math.cos(al) ** 2
        - math.cos(be) ** 2
        - math.cos(ga) ** 2
        + 2.0 * math.cos(al) * math.cos(be) * math.cos(ga)
    )

    # Reciprocal-lattice magnitudes (include the 2π factor).
    a_star = 2.0 * math.pi * math.sin(al) / (a * D)
    b_star = 2.0 * math.pi * math.sin(be) / (b * D)
    c_star = 2.0 * math.pi * math.sin(ga) / (c * D)

    # Reciprocal-cell angles via standard crystallographic identities.
    cos_beta_star = (math.cos(ga) * math.cos(al) - math.cos(be)) / (
        math.sin(ga) * math.sin(al)
    )
    cos_gamma_star = (math.cos(al) * math.cos(be) - math.cos(ga)) / (
        math.sin(al) * math.sin(be)
    )
    # sin(β*) and sin(γ*) are non-negative by convention.
    sin_beta_star = math.sqrt(max(0.0, 1.0 - cos_beta_star**2))
    sin_gamma_star = math.sqrt(max(0.0, 1.0 - cos_gamma_star**2))

    return np.array(
        [
            [a_star, b_star * cos_gamma_star, c_star * cos_beta_star],
            [0.0, b_star * sin_gamma_star, -c_star * sin_beta_star * math.cos(al)],
            [0.0, 0.0, 2.0 * math.pi / c],
        ]
    )


def reciprocal_vectors(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return the three reciprocal-lattice vectors (b1\\*, b2\\*, b3\\*)
    expressed in the BL1967 crystal-Cartesian frame.

    The vectors are the columns of :func:`b_matrix_bl1967` — see that
    function's docstring for the full convention and the cross-solver
    equivalence note.  Magnitudes carry the 2π factor (BL1967 / SPEC /
    hkl_soleil default convention).

    Orthogonality with the direct-lattice vectors holds::

        bᵢ\\* · aⱼ = 2π δᵢⱼ

    where ``aⱼ`` is the corresponding direct-lattice vector returned by
    :func:`lattice_vectors` (also expressed in the BL1967
    crystal-Cartesian frame).

    Parameters
    ----------
    a, b, c : float
        Direct-lattice parameters in Angstroms.
    alpha, beta, gamma : float
        Direct-lattice angles in degrees.

    Returns
    -------
    b1, b2, b3 : numpy.ndarray, shape (3,)
        Reciprocal-lattice vectors in inverse Angstroms (with 2π factor),
        expressed as columns of the BL1967 B matrix.
    """
    B = b_matrix_bl1967(a, b, c, alpha, beta, gamma)
    return B[:, 0], B[:, 1], B[:, 2]


def lattice_vectors(
    a: float,
    b: float,
    c: float,
    alpha: float,
    beta: float,
    gamma: float,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """
    Return the three direct-lattice vectors (a1, a2, a3) expressed in the
    BL1967 crystal-Cartesian frame (the same frame as :func:`b_matrix_bl1967`).

    Derived from the reciprocal-of-reciprocal identity
    ``aᵢ · bⱼ\\* = 2π δᵢⱼ``::

        [a1 | a2 | a3] = 2π · (B^T)^{-1} = 2π · (B^{-1})^T

    For orthogonal cells (cubic, tetragonal, orthorhombic) this gives
    ``a1 = a x̂``, ``a2 = b ŷ``, ``a3 = c ẑ`` — equivalent to the
    "direct-lattice-a-along-x" frame used by earlier versions of this
    package.  For non-orthogonal cells the direct vectors are no longer
    placed with ``a1 ∥ x̂``; instead ``b1\\* ∥ x̂`` (by the BL1967
    convention), and the direct vectors land wherever the duality
    relation puts them.  Magnitudes and angles between direct vectors
    are unchanged (``|a1| = a``, ``|a2| = b``, ``|a3| = c``,
    ``angle(a2, a3) = α``, etc.).

    Parameters
    ----------
    a, b, c : float
        Direct-lattice parameters in Angstroms.
    alpha, beta, gamma : float
        Direct-lattice angles in degrees.

    Returns
    -------
    a1, a2, a3 : numpy.ndarray, shape (3,)
        Direct-lattice vectors in Angstroms, expressed in the BL1967
        crystal-Cartesian frame.
    """
    B = b_matrix_bl1967(a, b, c, alpha, beta, gamma)
    A = 2.0 * math.pi * np.linalg.inv(B).T
    return A[:, 0], A[:, 1], A[:, 2]


def b_matrix(
    b1: np.ndarray,
    b2: np.ndarray,
    b3: np.ndarray,
) -> np.ndarray:
    """
    Assemble the B matrix from the reciprocal-lattice vectors.

    With the BL1967 / SPEC / hkl_soleil convention used by this
    package, the B matrix transforms Miller indices ``h`` = (h, k, l)
    to the scattering vector in the BL1967 crystal-Cartesian frame
    (Busing & Levy 1967 eq. 3)::

        Q_c = B @ h
        B = [b1, b2, b3]    (column-stack of the reciprocal vectors)

    Hence ``B @ h = h·b1 + k·b2 + l·b3``.  The orthogonality condition
    is ``bᵢ · aⱼ = 2π δᵢⱼ``, and ``|B @ h| = 2π / d_hkl``.

    For non-orthogonal lattices the columns are not mutually
    orthogonal; the B matrix is upper-triangular by construction (a\\*
    along crystal-x, b\\* in crystal-xy plane, c\\* completing the
    right-handed triple).  See :func:`b_matrix_bl1967` for the
    closed-form construction and the cross-solver equivalence note.

    Parameters
    ----------
    b1, b2, b3 : numpy.ndarray, shape (3,)
        Reciprocal-lattice vectors in inverse Angstroms (with 2π factor),
        as returned by :func:`reciprocal_vectors`.

    Returns
    -------
    B : numpy.ndarray, shape (3, 3)
        B matrix in inverse Angstroms (with 2π factor, BL1967 / SPEC /
        hkl_soleil convention).
    """
    return np.column_stack([b1, b2, b3])
