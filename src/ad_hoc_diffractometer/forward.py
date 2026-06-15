# Copyright (c) 2026-2026 UChicago Argonne, LLC
# SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
r"""
forward.py — Forward diffraction calculation: (h, k, l) → motor angles.

The forward calculation is the inverse of ``inverse()``: given a reciprocal-
lattice point (h, k, l) and a diffraction mode, compute all valid sets of
motor angles that satisfy the Bragg condition.

Entry point
-----------
``AdHocDiffractometer.forward(h, k, l)`` calls ``compute_forward(geometry,
h, k, l)`` which dispatches to the appropriate solver based on the active
mode type.

Supported modes
---------------
Bisecting (ConstraintSet with BisectConstraint)
    Classic Eulerian four-circle bisecting solution (Busing & Levy 1967,
    section "Angle settings").  Valid for psic (eta/delta bisecting),
    fourcv, fourch (omega/ttheta bisecting), and analogous geometries.
    Stage names are read directly from the BisectConstraint — no
    axis-geometry heuristic is used.

    Algorithm:
        1. Q_phi = UB @ (h, k, l)             — target in phi frame
        2. \|Q\| → ttheta via Bragg's law
        3. fixed sample stages set from SampleConstraint values
        4. detector_stage = ttheta            (from BisectConstraint.detector_stage)
        5. sample_stage   = ttheta / 2        (from BisectConstraint.sample_stage)
        6. Remaining free sample stages (chi, phi or kchi, kphi) solved
           from the direction of Q_phi, choosing among the standard
           solution branches (two solutions: chi in [0°,180°] and
           chi in [-180°, 0°]).

True virtual bisecting (ConstraintSet with VirtualBisectConstraint)
    The physically correct bisecting condition for kappa diffractometers
    (kappa4cv, kappa4ch, kappa6c).  The bisecting condition is enforced
    on the *virtual* Eulerian omega pseudoangle (Walko 2016, eq. [16]):
    ``omega_virtual = ttheta/2``.  The literal ``komega = ttheta/2``
    enforced by :class:`~mode.BisectConstraint` is only an approximation
    and is replaced for kappa-bisect modes by
    :class:`~mode.VirtualBisectConstraint`.  See issue #226.

    Algorithm (:func:`_solve_bisecting_kappa_virtual`):
        1. Build a synthetic mode that fixes ``omega_virtual = ttheta/2``
           via a :class:`~mode.SampleConstraint`.
        2. Delegate to :func:`~kappa.solve_kappa_virtual`, which runs
           Newton iteration on virtual ``(chi, phi)`` with the omega
           pseudoangle frozen.
        3. For each converged virtual triple, convert to real motors via
           :func:`~kappa.eulerian_to_kappa` for both kappa branches
           (``±1``).
        4. Polish the result with a finer central-difference Newton on
           ``(chi, phi)`` to push the residual toward machine precision.

    Some reflections that were "accessible" under the previous
    ``komega = ttheta/2`` approximation are **not** accessible under
    true virtual bisecting — those previous solutions had
    ``omega_virtual ≠ ttheta/2`` and were therefore not physically
    bisecting.

Fixed sample angle (ConstraintSet without BisectConstraint)
    One or more sample stages are frozen at declared values.  The
    remaining free stages are solved numerically from the direction
    of Q_phi.

Raises
------
NotImplementedError
    If the active mode is None or is not a supported type for this geometry.
ValueError
    If wavelength or UB matrix is not set.
ValueError
    If (h, k, l) = (0, 0, 0).
EwaldSphereViolation
    If the requested ``|Q|`` exceeds the Ewald sphere (wavelength too long).
ConstraintViolation
    If a returned solution violates a declared constraint beyond tolerance.

References
----------
* Busing & Levy, Acta Cryst. 22, 457-464 (1967) — angle-setting equations
* You, J. Appl. Cryst. 32, 614-623 (1999) — psic geometry
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# ForwardContext — pre-computed intermediates for the inner Newton loop
# ---------------------------------------------------------------------------


class ForwardContext:
    """
    Pre-computed intermediates for the forward solver's inner loops.

    Created once per :func:`~ad_hoc_diffractometer.forward.compute_forward()`
    call, this bundles all constant quantities needed by the Newton-Raphson
    residual evaluations so they are not recomputed on every call to
    :func:`~ad_hoc_diffractometer.orientation.angles_to_phi_vector()`.

    Attributes
    ----------
    sample_stages : list of Stage
    detector_stages : list of Stage
    two_pi_over_lambda : float
    y_eff : numpy.ndarray, shape (3,)
        Effective beam direction (R_inc.T @ y_hat).
    """

    def __init__(self, geometry):
        self.sample_stages = geometry.sample_stages
        self.detector_stages = geometry.detector_stages
        self.two_pi_over_lambda = 2.0 * math.pi / geometry.wavelength

        y_hat = np.asarray(geometry.basis["longitudinal"], dtype=float)
        y_norm = float(np.linalg.norm(y_hat))
        y_hat = y_hat / y_norm
        R_inc = geometry.inclination_matrix
        self.y_eff = R_inc.T @ y_hat

        # Cached rotation matrices (populated by prepare_bisecting)
        self._cached_Z_prefix: np.ndarray | None = None
        self._free_sample_indices: list[int] | None = None
        self._cached_D: np.ndarray | None = None

    def prepare_caching(
        self,
        fixed_angles: dict[str, float],
        free_stage_names: set[str],
    ) -> None:
        """
        Pre-compute rotation matrices for stages whose angles will not
        change during Newton iteration.

        Parameters
        ----------
        fixed_angles : dict
            All angles including fixed ones.
        free_stage_names : set of str
            Names of stages that will vary during iteration.
        """
        from .rotation import _rotation_matrix_normalized

        # Detector: if no detector stage is free, cache D entirely.
        # Composition: BL1967 standard order (Busing & Levy 1967),
        # outermost (floor-most) stage leftmost,
        # so D = R_0 @ R_1 @ ... @ R_{n-1}.
        det_free = any(s.name in free_stage_names for s in self.detector_stages)
        if not det_free:
            D = np.eye(3)
            for s in self.detector_stages:
                angle = fixed_angles.get(s.name, s.angle)
                D = D @ _rotation_matrix_normalized(s._axis_hat, angle)  # noqa: SLF001
            self._cached_D = D
        else:  # pragma: no cover
            # No remaining call site passes a detector stage as free
            # after the kappa virtual-angle solver was rewritten in
            # issue #241; retained as a defensive path for callers
            # that build a ForwardContext directly.
            self._cached_D = None

        # Sample: find the first free stage index and cache the prefix product
        # of the fixed *outer* stages (indices 0..first_free-1).  This prefix
        # is built in the same outermost-leftmost order:
        # Z_prefix = R_0 @ R_1 @ ... @ R_{first_free-1}.
        self._free_sample_indices = []
        for i, s in enumerate(self.sample_stages):
            if s.name in free_stage_names:
                self._free_sample_indices.append(i)

        if self._free_sample_indices:
            first_free = self._free_sample_indices[0]
            if first_free > 0:
                Z_prefix = np.eye(3)
                for i in range(first_free):
                    s = self.sample_stages[i]
                    angle = fixed_angles.get(s.name, s.angle)
                    Z_prefix = Z_prefix @ _rotation_matrix_normalized(
                        s._axis_hat, angle
                    )  # noqa: SLF001
                self._cached_Z_prefix = Z_prefix
            else:
                self._cached_Z_prefix = np.eye(3)
        else:
            self._cached_Z_prefix = None

    def q_phi(self, angles: dict[str, float]) -> np.ndarray:
        """Compute Q_phi using cached rotation matrices where possible."""
        from .orientation import _compute_q_phi_cached

        return _compute_q_phi_cached(
            self.sample_stages,
            self.detector_stages,
            angles,
            self.two_pi_over_lambda,
            self.y_eff,
            self._cached_Z_prefix,
            self._free_sample_indices,
            self._cached_D,
        )

    def q_phi_uncached(self, angles: dict[str, float]) -> np.ndarray:
        """Compute Q_phi without any caching (for validation)."""
        from .orientation import _compute_q_phi

        return _compute_q_phi(
            self.sample_stages,
            self.detector_stages,
            angles,
            self.two_pi_over_lambda,
            self.y_eff,
        )

    def jacobian_analytic(
        self,
        angles: dict[str, float],
        free_names: list[str],
    ) -> np.ndarray:
        """
        Compute the analytic Jacobian of Q_phi with respect to free stage angles.

        Uses the closed-form derivative of the Rodrigues rotation matrix to
        avoid finite-difference evaluations entirely.  The detector rotation
        ``D`` must be cached (all detector stages fixed) — this is the case
        for every call site that passes a ``ForwardContext`` to
        :func:`_solve_two_angles`.

        Parameters
        ----------
        angles : dict[str, float]
            Current motor angles (degrees) for all stages.
        free_names : list of str
            Ordered list of free stage names (length 1 or 2).

        Returns
        -------
        J : numpy.ndarray, shape (3, n_free)
            ``J[:, k]`` is ``dQ_phi / d(stage_k angle in degrees)``.
        """
        from .rotation import _rotation_matrix_and_derivative_normalized
        from .rotation import _rotation_matrix_normalized

        n_free = len(free_names)
        free_set = set(free_names)

        # Q_lab is constant (D and y_eff are fixed).
        D = self._cached_D if self._cached_D is not None else np.eye(3)
        Q_lab = self.two_pi_over_lambda * (D @ self.y_eff - self.y_eff)

        # Build per-stage rotation matrices for all sample stages from the
        # first free index onward (stages before that are in cached_Z_prefix).
        stages = self.sample_stages
        first_free = (
            self._free_sample_indices[0] if self._free_sample_indices else len(stages)
        )

        # Collect R_i and (for free stages) dR_i for indices >= first_free.
        n_tail = len(stages) - first_free
        R_list: list[np.ndarray] = [np.empty((3, 3))] * n_tail
        dR_map: dict[int, np.ndarray] = {}  # absolute index -> dR

        for idx in range(first_free, len(stages)):
            s = stages[idx]
            angle = angles.get(s.name, s.angle)
            if s.name in free_set:
                R_i, dR_i = _rotation_matrix_and_derivative_normalized(
                    s._axis_hat,
                    angle,  # noqa: SLF001
                )
                R_list[idx - first_free] = R_i
                dR_map[idx] = dR_i
            else:
                R_list[idx - first_free] = _rotation_matrix_normalized(
                    s._axis_hat,
                    angle,  # noqa: SLF001
                )

        # Composition (BL1967 standard convention, outermost-leftmost):
        # the full sample rotation is
        #
        #   Z = R_0 @ R_1 @ ... @ R_{N-1}
        #     = Z_prefix @ R_list[0] @ R_list[1] @ ... @ R_list[n_tail-1]
        #
        # where Z_prefix = R_0 @ ... @ R_{first_free-1} is cached and the
        # tail R_list[0..n_tail-1] holds R_{first_free}..R_{N-1}.
        #
        # The derivative of Z w.r.t. the angle of the free stage at
        # relative index k is
        #
        #   dZ/dθ_k = before[k] @ dR_list[k] @ after[k]
        #
        # with
        #
        #   before[k] = Z_prefix @ R_list[0] @ R_list[1] @ ... @ R_list[k-1]
        #   after[k]  = R_list[k+1] @ R_list[k+2] @ ... @ R_list[n_tail-1]
        #
        # Build them by simple forward / backward sweeps.
        Z_prefix = (
            self._cached_Z_prefix if self._cached_Z_prefix is not None else np.eye(3)
        )

        # before[k] for k = 0..n_tail
        before = [np.empty((3, 3))] * (n_tail + 1)
        before[0] = Z_prefix
        for k in range(n_tail):
            before[k + 1] = before[k] @ R_list[k]

        # after[k] for k = 0..n_tail (after[k] is the product of R_list[k+1..])
        after = [np.empty((3, 3))] * (n_tail + 1)
        after[n_tail] = np.eye(3)
        for k in range(n_tail - 1, -1, -1):
            after[k] = R_list[k + 1] @ after[k + 1] if (k + 1) < n_tail else np.eye(3)

        # deg2rad factor: dR is w.r.t. radians, but angles are in degrees.
        deg2rad = np.pi / 180.0

        J = np.zeros((3, n_free))
        for col, name in enumerate(free_names):
            # Find the absolute stage index for this free name
            abs_idx = next(i for i, s in enumerate(stages) if s.name == name)
            k = abs_idx - first_free  # relative index into R_list
            dR_k = dR_map[abs_idx]
            # dZ/dθ = before[k] @ dR_k @ after[k]
            dZ = before[k] @ dR_k @ after[k]
            # dQ_phi/dθ = dZ^T @ Q_lab, with deg2rad conversion
            J[:, col] = deg2rad * (dZ.T @ Q_lab)

        return J


# ---------------------------------------------------------------------------
# Public entry point
# ---------------------------------------------------------------------------


def compute_forward(
    geometry: AdHocDiffractometer,
    h: float,
    k: float,
    l: float,  # noqa: E741
) -> list[dict[str, float]]:
    r"""
    Compute all valid motor-angle solutions for the reciprocal-lattice point
    (h, k, l) in the geometry's active diffraction mode.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer.  Must have ``wavelength`` and ``sample.UB`` set,
        and an active ``mode_name``.
    h, k, l : float
        Miller indices of the target reflection.

    Returns
    -------
    list of dict[str, float]
        Each element is a complete set of motor angles (all stage names as
        keys, values in degrees) that satisfies the Bragg condition under
        the active mode constraints.  May contain zero, one, or multiple
        solutions depending on the geometry and mode.

    Raises
    ------
    ValueError
        If ``geometry.wavelength`` is None.
    ValueError
        If ``geometry.sample.UB`` is None.
    ValueError
        If (h, k, l) == (0, 0, 0).
    EwaldSphereViolation
        If the scattering vector magnitude exceeds the Ewald sphere
        (``|Q|`` > 4π/λ, i.e. Bragg condition cannot be satisfied).
    ConstraintViolation
        If a solution returned by the solver violates a declared constraint
        beyond the display-precision tolerance.
    NotImplementedError
        If no active mode is set or the active mode type is not supported
        for this geometry.
    """
    from .mode import ConstraintSet
    from .mode import EwaldSphereViolation

    # --- Precondition checks ------------------------------------------------

    if geometry.wavelength is None:
        raise ValueError(
            "geometry.wavelength must be set before calling forward(). "
            "Set it with e.g. geometry.wavelength = 1.5406."
        )

    sample = geometry.sample
    if sample.UB is None:
        raise ValueError(
            f"Sample {sample.name!r} has no UB matrix. "
            "Set one with ub_identity(), ub_from_one_reflection(), or "
            "assign sample.UB directly."
        )

    hkl = np.array([float(h), float(k), float(l)], dtype=float)
    if np.allclose(hkl, 0.0):
        raise ValueError("forward(): (h, k, l) = (0, 0, 0) is not a valid reflection.")

    if geometry.mode_name is None:
        raise NotImplementedError(
            f"Geometry {geometry.name!r} has no active diffraction mode. "
            "Set geometry.mode_name to one of: "
            f"{sorted(geometry.modes.keys())}."
        )

    # --- Compute target scattering vector in phi frame ----------------------

    Q_phi = sample.UB @ hkl  # shape (3,)
    Q_mag = float(np.linalg.norm(Q_phi))

    # Bragg's law: |Q| = 4π sin(θ) / λ  →  sin(θ) = |Q|λ / (4π)
    wavelength = geometry.wavelength
    sin_theta = Q_mag * wavelength / (4.0 * math.pi)
    if sin_theta > 1.0:
        raise EwaldSphereViolation(
            q_mag=Q_mag,
            q_max=4.0 * math.pi / wavelength,
            wavelength=wavelength,
        )

    ttheta_rad = 2.0 * math.asin(sin_theta)
    ttheta_deg = math.degrees(ttheta_rad)

    # --- Dispatch to mode-specific solver ------------------------------------

    mode = geometry.mode

    if not isinstance(mode, ConstraintSet):  # pragma: no cover
        raise NotImplementedError(
            f"forward(): mode {geometry.mode_name!r} "
            f"({type(mode).__name__}) is not a ConstraintSet. "
            f"All modes must be ConstraintSet instances."
        )

    if not mode.is_implemented(geometry):
        raise NotImplementedError(
            f"forward(): mode {geometry.mode_name!r} is not yet implemented "
            f"for geometry {geometry.name!r}.  "
            "Check mode.is_implemented(geometry) before calling forward()."
        )

    solutions = _solve_constraint_set(geometry, Q_phi, ttheta_deg, mode)
    _validate_solutions(solutions, mode, geometry)
    _populate_output_extras(geometry, mode, solutions)
    return solutions


# ---------------------------------------------------------------------------
# Output-extras population (issue #292)
# ---------------------------------------------------------------------------


# Map of output-slot key → (callable(geometry, angles) → float, friendly label).
# The callables are imported lazily inside :func:`_populate_output_extras` to
# avoid a circular import at module load time (``reference`` imports from
# ``forward``).
_OUTPUT_EXTRA_KEYS = ("alpha_i", "beta_out", "psi", "omega")


def _populate_output_extras(
    geometry: AdHocDiffractometer,
    mode,
    solutions: list[dict[str, float]],
) -> None:
    """Populate output-slot extras (alpha_i, beta_out, psi, omega) per solution.

    Issue #292.  A subset of declarative modes (psic, sixc, zaxis, s2d2 surface
    modes; the fixed_psi_* family; fixed_omega_*) declare placeholder slots
    in ``mode.extras`` for derived angles that the forward solver does not
    constrain directly.  Before this hook those slots remained at their YAML
    default of ``None`` even after a successful ``forward()`` call.  This
    function fills each declared slot with a list of values aligned with
    ``solutions`` (one float per solution), computed via the corresponding
    helper in :mod:`ad_hoc_diffractometer.reference`.

    Behavior
    --------
    * Only keys actually declared in ``mode.extras`` are touched.
    * A key declared but whose required reference vector is unset on the
      geometry (e.g. ``alpha_i`` without ``surface_normal``) is left as
      ``None``; a debug-level log message records why.
    * Empty ``solutions`` leaves every slot as an empty list.
    * Each successful call **replaces** the prior contents of the slot.
    """
    if mode is None or not getattr(mode, "extras", None):
        return

    relevant = [k for k in _OUTPUT_EXTRA_KEYS if k in mode.extras]
    if not relevant:
        return

    # Lazy imports — ``reference`` imports from ``forward``.
    from .reference import exit_angle as _exit_angle
    from .reference import incidence_angle as _incidence_angle
    from .reference import omega_pseudo as _omega_pseudo
    from .reference import psi_angle as _psi_angle

    computers: dict[str, callable] = {
        "alpha_i": _incidence_angle,
        "beta_out": _exit_angle,
        "psi": _psi_angle,
        "omega": _omega_pseudo,
    }

    for key in relevant:
        compute = computers[key]
        values: list[float] = []
        failure: Exception | None = None
        for angles in solutions:
            try:
                values.append(float(compute(geometry, angles=angles)))
            except Exception as exc:  # noqa: BLE001
                # Underlying call raised (e.g. missing surface_normal /
                # azimuth, or psi undefined when Q ∥ n_ref).
                # Leave the slot unpopulated and record the cause for a
                # single debug log below.  We break immediately so a
                # later good solution does not mask the failure.
                failure = exc
                values = []
                break
        if values:
            mode.extras[key] = values
        else:
            # Reset to None so a stale value from a previous forward() call
            # is not retained, and report the cause once.
            mode.extras[key] = None if solutions else []
            if failure is not None:
                logger.debug(
                    "_populate_output_extras: leaving mode.extras[%r] as None "
                    "for mode %r on geometry %r: %s",
                    key,
                    geometry.mode_name,
                    geometry.name,
                    failure,
                )


# ---------------------------------------------------------------------------
# Constraint-set dispatcher
# ---------------------------------------------------------------------------


def _solve_constraint_set(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Dispatch to the appropriate analytic or numeric solver based on the
    constraint pattern encoded in the ConstraintSet.

    Currently supported patterns:

    - Bisect (any geometry): a :class:`~mode.BisectConstraint` is present,
      with an optional ``DetectorConstraint`` freezing an outer detector
      stage and optional fixed ``SampleConstraint`` values.
    - Fixed sample angle without bisect: one or more fixed sample angles;
      delegates to _solve_fixed_sample.

    Unsupported patterns fall back to the numeric solver.
    """

    # Psi-constant mode (ReferenceConstraint("psi") validation filter).
    # Must be checked FIRST: psic/kappa6c fixed_psi_vertical has a
    # BisectConstraint but needs psi validation before bisecting.
    if _is_psi_mode(geometry, mode):
        return _solve_psi_mode(geometry, Q_phi, ttheta_deg, mode)

    # Double-diffraction mode (simultaneous Bragg for two reflections).
    # Must be checked before has_bisect: double_diffraction modes have a
    # BisectConstraint but use a 4D solver instead of plain bisecting.
    if _is_double_diffraction_mode(geometry, mode):
        return _solve_double_diffraction(geometry, Q_phi, ttheta_deg, mode)

    # Zone mode (You 1999 §6, SPEC `setmode 5`):  Q is confined to the
    # plane spanned by two reciprocal-lattice vectors z0, z1.  Like
    # double_diffraction, this is dispatched by extras keys, not by a
    # dedicated constraint type.
    if _is_zone_mode(geometry, mode):
        return _solve_zone(geometry, Q_phi, ttheta_deg, mode)

    if mode.has_bisect:
        return _solve_bisecting(geometry, Q_phi, ttheta_deg, mode)

    # Virtual kappa angle mode (omega, chi, phi on a kappa geometry).
    if _is_kappa_virtual_mode(geometry, mode):
        return _solve_kappa_virtual(geometry, Q_phi, ttheta_deg, mode)

    # Surface diffraction mode (ReferenceConstraint with surface_normal).
    if _is_surface_mode(geometry, mode):
        return _solve_surface(geometry, Q_phi, ttheta_deg, mode)

    # OMEGA pseudo-angle mode (SPEC psic ``omega-fixed`` family, #264).
    # Must be checked after _is_surface_mode (whose predicate explicitly
    # excludes the "omega" reference name) and before the qaz/fallback
    # branches.
    if _is_omega_mode(geometry, mode):
        return _solve_omega_mode(geometry, Q_phi, ttheta_deg, mode)

    # qaz detector constraint mode (lifting_detector_* family).
    if _is_qaz_mode(geometry, mode):
        return _solve_qaz_mode(geometry, Q_phi, ttheta_deg, mode)

    # Free-detectors mode (issue #264 — both detector stages float to
    # satisfy Bragg jointly with the remaining sample stages, optionally
    # with one ReferenceConstraint such as alpha_i).
    if _is_free_detectors_mode(geometry, mode):
        return _solve_free_detectors(geometry, Q_phi, ttheta_deg, mode)

    # Fixed-angle only (no bisect).
    return _solve_fixed_sample(geometry, Q_phi, ttheta_deg, mode)


# ---------------------------------------------------------------------------
# Mode-specific solvers
# ---------------------------------------------------------------------------


def _is_standard_eulerian_pair(chi_stage, phi_stage) -> bool:
    """
    Return True if the chi and phi stages have orthogonal axis vectors.

    Standard Eulerian geometries (fourcv, fourch, psic, sixc, fivec) have
    chi about the longitudinal axis and phi about a transverse or vertical
    axis — always orthogonal.  Kappa geometries do NOT have orthogonal
    chi/phi (the kappa axis is tilted), so they will return False here and
    use the Newton fallback.

    Parameters
    ----------
    chi_stage, phi_stage : Stage
        The two free sample stages (chi-role and phi-role).

    Returns
    -------
    bool
    """
    dot = abs(float(np.dot(chi_stage._axis_hat, phi_stage._axis_hat)))  # noqa: SLF001
    return dot < 1e-8


def _solve_bisecting_analytic(
    ctx: ForwardContext,
    chi_stage,
    phi_stage,
    angles: dict[str, float],
    Q_phi_target: np.ndarray,
) -> list[tuple[float, float]]:
    r"""
    Analytic solver for bisecting with two free Eulerian angles.

    Derives (chi, phi) from a direct trigonometric decomposition when the
    chi and phi axes are orthogonal (standard Eulerian geometry).

    **Derivation.**

    Under the BL1967 standard convention (Busing & Levy 1967,
    outermost stage leftmost), the sample rotation matrix is:

    .. math::

        Z = Z_{\text{prefix}} \cdot R_\chi \cdot R_\phi

    where :math:`Z_{\text{prefix}}` is the product of all fixed outer
    sample stages.  The scattering vector in the phi frame is:

    .. math::

        Q_\phi = Z^T Q_{\text{lab}}
              = R_\phi^T \cdot R_\chi^T \cdot Z_{\text{prefix}}^T \cdot Q_{\text{lab}}

    Equivalently, multiplying both sides by :math:`R_\chi R_\phi`:

    .. math::

        R(\hat n_\chi, \chi) \cdot R(\hat n_\phi, \phi) \cdot Q_{\phi,\text{target}}
            = Z_{\text{prefix}}^T \cdot Q_{\text{lab}}

    Define the known vectors:

    - :math:`q = Z_{\text{prefix}}^T \cdot Q_{\text{lab}}`
    - :math:`v = Q_{\phi,\text{target}}`

    The equation to solve becomes:

    .. math::

        R(\hat n_\chi, \chi) \cdot R(\hat n_\phi, \phi) \cdot v = q

    Since :math:`\hat n_\chi \perp \hat n_\phi`, define
    :math:`\hat n_3 = \hat n_\chi \times \hat n_\phi`.

    **Step 1 — solve for phi.**  Project the equation onto
    :math:`\hat n_\chi` (which is invariant under
    :math:`R(\hat n_\chi, \chi)`):

    .. math::

        \hat n_\chi \cdot q
            = \hat n_\chi \cdot R(\hat n_\phi, \phi) \cdot v

    By Rodrigues, since :math:`\hat n_\chi \perp \hat n_\phi`:

    .. math::

        \hat n_\chi \cdot R(\hat n_\phi, \phi) \cdot v
            = (\hat n_\chi \cdot v) \cos\phi
            + (\hat n_\chi \cdot (\hat n_\phi \times v)) \sin\phi
            = (\hat n_\chi \cdot v) \cos\phi
            + (\hat n_3 \cdot v) \sin\phi

    so :math:`A \cos\phi + B \sin\phi = C` with
    :math:`A = \hat n_\chi \cdot v`,
    :math:`B = \hat n_3 \cdot v`,
    :math:`C = \hat n_\chi \cdot q`, yielding two :math:`\phi`
    candidates.

    **Step 2 — for each phi, solve for chi.**  Compute
    :math:`w = R(\hat n_\phi, \phi) \cdot v` (known once :math:`\phi`
    is fixed).  The equation :math:`R(\hat n_\chi, \chi) \cdot w = q`
    is a single rotation about :math:`\hat n_\chi` taking :math:`w`
    onto :math:`q`.  Their :math:`\hat n_\chi`-components must agree
    (guaranteed by construction); :math:`\chi` is the angle in the
    plane perpendicular to :math:`\hat n_\chi` between the
    perpendicular projections of :math:`w` and :math:`q`.

    Parameters
    ----------
    ctx : ForwardContext
        Must have ``prepare_caching`` already called with the chi/phi
        stages as free.
    chi_stage, phi_stage : Stage
        The two free sample stages (orthogonal axes).
    angles : dict[str, float]
        Baseline angles with all fixed stages set.
    Q_phi_target : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (Å⁻¹).

    Returns
    -------
    list of (chi_deg, phi_deg)
        Candidate solutions.  May contain 0, 1, or 2 entries.
        All entries have been validated against a single residual check.
    """

    # Axis unit vectors
    n_chi = chi_stage._axis_hat  # noqa: SLF001
    n_phi = phi_stage._axis_hat  # noqa: SLF001

    # Known vectors (BL1967 standard outermost-leftmost convention).
    Z_prefix = ctx._cached_Z_prefix
    D = ctx._cached_D if ctx._cached_D is not None else np.eye(3)

    # Q_lab: scattering vector in the lab frame (all detector & bisect angles fixed)
    Q_lab = ctx.two_pi_over_lambda * (D @ ctx.y_eff - ctx.y_eff)

    # q = Z_prefix^T @ Q_lab  (right-hand side of R_chi R_phi v = q)
    q = Z_prefix.T @ Q_lab
    # v = target Q_phi (left-hand side argument)
    v = np.asarray(Q_phi_target, dtype=float)

    n3 = np.cross(n_chi, n_phi)

    # --- Step 1: solve A cos phi + B sin phi = C for two phi values ---
    A = float(np.dot(n_chi, v))
    B = float(np.dot(n3, v))
    C = float(np.dot(n_chi, q))

    amp = math.sqrt(A * A + B * B)
    if amp < 1e-12:
        # v ∥ n_phi (perpendicular component of v in (n_chi, n_3) plane is
        # zero).  phi is indeterminate; defer to the Newton fallback.
        return []

    cos_arg_phi = C / amp
    if abs(cos_arg_phi) > 1.0 + 1e-8:  # pragma: no cover
        return []  # no real phi solution

    cos_arg_phi = max(-1.0, min(1.0, cos_arg_phi))
    phi0 = math.atan2(B, A)
    delta_phi = math.acos(cos_arg_phi)
    phi_candidates_rad = [phi0 + delta_phi, phi0 - delta_phi]

    # --- Step 2: for each phi, recover chi from the residual rotation ---
    # R(n_chi, chi) @ w = q with w = R(n_phi, phi) @ v.  chi is the
    # signed angle about n_chi taking w_perp onto q_perp in the plane
    # perpendicular to n_chi.
    raw_candidates = []
    for phi_rad in phi_candidates_rad:
        phi_d = math.degrees(phi_rad)
        # Build w = R(n_phi, phi) @ v via Rodrigues:
        #   w = v cos φ + (n_phi × v) sin φ + n_phi (n_phi · v)(1 − cos φ)
        c_ph = math.cos(phi_rad)
        s_ph = math.sin(phi_rad)
        nphi_v = float(np.dot(n_phi, v))
        w = v * c_ph + np.cross(n_phi, v) * s_ph + n_phi * nphi_v * (1.0 - c_ph)

        # Project w and q onto the plane perpendicular to n_chi.
        w_par = float(np.dot(w, n_chi))
        q_par = float(np.dot(q, n_chi))
        w_perp = w - w_par * n_chi
        q_perp = q - q_par * n_chi
        w_perp_norm = float(np.linalg.norm(w_perp))
        q_perp_norm = float(np.linalg.norm(q_perp))
        if w_perp_norm < 1e-12 or q_perp_norm < 1e-12:  # pragma: no cover
            # Degenerate: w or q is parallel to n_chi; chi indeterminate.
            # No shipped geometry / mode / reflection lands here after the
            # issue-#284 kappa equivalent-Eulerian chi axis correction;
            # retained as a defensive fallback for ad-hoc geometries
            # whose target Q_phi happens to project to zero on the
            # plane perpendicular to n_chi.
            chi_d = 0.0
        else:
            # cos chi = (w_perp · q_perp) / (|w_perp| |q_perp|)
            # sin chi = n_chi · (w_perp × q_perp) / (|w_perp| |q_perp|)
            cos_chi = float(np.dot(w_perp, q_perp)) / (w_perp_norm * q_perp_norm)
            sin_chi = float(np.dot(n_chi, np.cross(w_perp, q_perp))) / (
                w_perp_norm * q_perp_norm
            )
            chi_d = math.degrees(math.atan2(sin_chi, cos_chi))

        raw_candidates.append((chi_d, phi_d))

    # Validate each candidate with a single residual evaluation.
    # Normalize angles to (-180, 180] before returning — this matches
    # the range that Newton iteration naturally produces and ensures
    # the angles are within typical stage limits.
    validated = []
    trial = dict(angles)
    chi_name = chi_stage.name
    phi_name = phi_stage.name

    for chi_d, phi_d in raw_candidates:
        trial[chi_name] = chi_d
        trial[phi_name] = phi_d

        # Compute Q_phi with these angles and check residual
        Q_computed = ctx.q_phi(trial)
        residual = float(np.linalg.norm(Q_computed - Q_phi_target))
        if residual < 1e-6:  # pragma: no branch
            # Normalize to (-180, 180]
            chi_d = (chi_d + 180.0) % 360.0 - 180.0
            phi_d = (phi_d + 180.0) % 360.0 - 180.0
            validated.append((chi_d, phi_d))

    # Sort so the "positive chi" branch (chi ≥ 0, or closer to 0°) comes
    # first.  This matches the ordering convention of the Newton solver and
    # ensures consistent branch selection in trajectory scans.
    validated.sort(key=lambda pair: -pair[0])  # descending chi → positive first

    return validated


def _solve_bisecting_kappa_virtual(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
    bisect_c,
    angles: dict[str, float],
) -> list[dict[str, float]]:
    r"""
    True virtual bisecting solver for kappa geometries.

    On a kappa diffractometer the real motor triple
    ``(komega, kappa, kphi)`` is parameterized by *virtual Eulerian
    pseudoangles* ``(omega, chi, phi)``.  *True* virtual bisecting
    means the **virtual** ``omega`` equals ``ttheta/2`` — not the
    previous approximation ``komega = ttheta/2``, which only coincides
    with true bisecting at ``kappa = 0``.

    Implementation
    --------------
    This function is a thin wrapper around
    :func:`~ad_hoc_diffractometer.kappa.solve_kappa_virtual`, which
    handles the entire kappa virtual-angle solve uniformly for both
    bisect (``VirtualBisectConstraint``) and fixed-virtual-angle
    (``SampleConstraint`` on ``omega``/``chi``/``phi``) modes via the
    geometry-aware :func:`~kappa.eulerian_to_kappa_axes` decomposition
    introduced by issue #241.

    The wrapper applies the kappa-side post-processing — cut-points,
    motor-limit checks, and final deduplication on ``(komega, kappa,
    kphi)`` — that the rest of the forward pipeline expects.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        A kappa diffractometer with ``kappa_pseudo_angle_convention``
        set.
    Q_phi : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (Å⁻¹).
    ttheta_deg : float
        Detector angle (2θ) in degrees.
    mode : ConstraintSet
        Active diffraction mode.
    bisect_c : :class:`~mode.VirtualBisectConstraint`
        Active virtual-bisecting constraint (parameter retained for
        signature compatibility with the dispatcher in
        :func:`_solve_bisecting`).
    angles : dict[str, float]
        Baseline angles dict with detector stage already frozen at
        ``ttheta_deg`` (parameter retained for signature
        compatibility).

    Returns
    -------
    list of dict[str, float]
        Valid solutions, possibly empty when the requested reflection
        is not in the bisecting locus or all candidate configurations
        are outside the stage limits.

    References
    ----------
    * D. A. Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016) —
      textbook pseudoangle special case.
    * W. R. Busing & H. A. Levy, *Acta Cryst.* 22, 457-464 (1967) —
      bisecting Eulerian geometry.
    * Issue #241 — geometry-aware decomposition that supersedes the
      Walko-only implementation.
    """
    from .kappa import kappa_to_eulerian_axes
    from .kappa import solve_kappa_virtual

    Q_phi_arr = np.asarray(Q_phi, dtype=float)

    raw = solve_kappa_virtual(geometry, Q_phi_arr, ttheta_deg, mode)

    convention = geometry.kappa_pseudo_angle_convention
    omega_target = ttheta_deg / 2.0

    # ``solve_kappa_virtual`` returns analytic, deduplicated kappa
    # motor triples (the geometry-aware decomposition introduced by
    # issue #241 is exact, and the inner solver deduplicates at 1e-6).
    # The wrapper therefore only needs to:
    #   1. enforce the virtual-bisect omega = ttheta/2 condition (the
    #      equivalent Eulerian solver returns both chi branches, only
    #      one of which is bisecting);
    #   2. apply mode cut-points; and
    #   3. drop solutions outside the hardware stage limits.
    # No residual Q-recheck or wrapper-side dedup is needed; both were
    # removed in issue #245 cleanup as unreachable post-#241 code.
    solutions: list[dict[str, float]] = []
    for sol in raw:
        merged = dict(angles)
        merged.update(sol)
        om_v, _, _ = kappa_to_eulerian_axes(
            merged["komega"], merged["kappa"], merged["kphi"], convention
        )
        residual = (om_v - omega_target + 180.0) % 360.0 - 180.0
        if abs(residual) > 1e-6:
            continue
        _apply_cut_points(merged, mode, geometry)
        if not _check_limits(geometry, merged):
            continue
        solutions.append(merged)
    return solutions


def _solve_bisecting(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Bisecting-mode forward solver for ConstraintSet.

    Implements the Busing & Levy (1967) angle-setting algorithm for an
    Eulerian geometry in the bisecting condition.

    The bisect stage pair (sample_stage, detector_stage) is read directly
    from the :class:`~mode.BisectConstraint` — no geometry heuristics are
    used.  ``bisect_constraint.detector_stage`` receives ``ttheta_deg``;
    ``bisect_constraint.sample_stage`` is set to ``ttheta_deg / 2``.

    Any ``DetectorConstraint`` in the mode (e.g. ``nu=0`` in psic) freezes
    that stage at its declared value.  Any fixed ``SampleConstraint`` values
    are applied before solving for the remaining free stages.

    Two solutions are returned corresponding to the "positive chi"
    branch (chi in [0, 180]) and the "negative chi" branch
    (chi in (-180, 0)).

    Both solutions are validated against stage limits; those that fail are
    dropped.  Cut-points from the mode are applied to each angle.

    **Dispatch**

    1. **Kappa true-virtual bisecting** — if the geometry is a kappa
       diffractometer and the bisect sample stage is ``komega``, the
       solver works in the *virtual* Eulerian frame (omega, chi, phi)
       where true bisecting means ``omega_virtual = ttheta/2``.  The
       analytic solver is invoked on the synthetic Eulerian triple, and
       each solution is converted to real ``(komega, kappa, kphi)`` motor
       angles via :func:`~kappa.eulerian_to_kappa` for both branches
       (kappa ≥ 0 and kappa ≤ 0).  See
       :func:`_solve_bisecting_kappa_virtual`.
    2. **Standard Eulerian fast path** — when the two free sample stages
       have orthogonal axes (fourcv, fourch, psic, sixc, fivec), the
       analytic solver :func:`_solve_bisecting_analytic` yields a
       5-10× speedup over Newton iteration.
    3. **Newton fallback** — for any non-standard axis configuration not
       covered by the fast paths above.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (Å⁻¹).
    ttheta_deg : float
        Detector angle (2θ) in degrees.
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
        Valid solutions (may be empty if all are out of limits).
    """
    from .mode import BisectConstraint

    sample_stages = geometry.sample_stages

    # Build the baseline angle dict: all stages at their current angles.
    angles: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }

    # Apply fixed sample constraints from the ConstraintSet.
    for sc in mode.fixed_sample_constraints:
        if sc.name in geometry._stages:  # noqa: SLF001
            angles[sc.name] = float(sc.value)

    # Freeze any explicitly constrained detector stage (e.g. nu=0 in psic).
    det_constraint = mode.detector_constraint
    if det_constraint is not None and not det_constraint.is_qaz:
        angles[det_constraint.name] = det_constraint.value

    # Read the bisect stage pair from the BisectConstraint — explicit, no heuristic.
    bisect_c = next(c for c in mode.constraints if isinstance(c, BisectConstraint))
    detector_name = bisect_c.detector_stage  # receives ttheta_deg
    sample_bisect_name = bisect_c.sample_stage  # receives ttheta_deg / 2

    angles[detector_name] = ttheta_deg

    # --- Kappa true-virtual bisecting fast path -----------------------------
    #
    # When the bisect constraint is a :class:`VirtualBisectConstraint`,
    # the bisecting condition is on the *virtual* Eulerian omega
    # pseudoangle (``omega_virtual = ttheta/2``), not on the real
    # ``komega`` motor.  Dispatch to a dedicated solver that works in
    # virtual Eulerian space and converts to real kappa motor angles via
    # :func:`~kappa.eulerian_to_kappa`.
    from .mode import VirtualBisectConstraint as _VBC

    if isinstance(bisect_c, _VBC):
        return _solve_bisecting_kappa_virtual(
            geometry, Q_phi, ttheta_deg, mode, bisect_c, angles
        )

    angles[sample_bisect_name] = ttheta_deg / 2.0

    # The remaining free sample stages are those not fixed by any constraint.
    constrained_names = set(mode.constrained_stages(geometry))
    constrained_names.add(detector_name)
    if det_constraint is not None and not det_constraint.is_qaz:
        constrained_names.add(det_constraint.name)
    free_sample = [s for s in sample_stages if s.name not in constrained_names]

    if len(free_sample) == 0:
        # All sample stages constrained — return the fixed solution if valid.
        _apply_cut_points(angles, mode, geometry)
        if _check_limits(geometry, angles):
            return [angles]
        return []

    if len(free_sample) == 1:
        # One free sample stage (e.g. phi in fixed_chi mode).
        # Solve 1D: find the single angle that satisfies Q_phi.
        return _solve_one_free_angle(geometry, angles, free_sample[0], Q_phi, mode)

    # The two remaining free stages are chi_stage and phi_stage.
    # By convention (BL1967, You1999) the stacking order is:
    #   [outermost ... chi_stage, phi_stage] (phi closest to sample)
    # The last two free stages in stacking order play the chi/phi roles.
    chi_stage = free_sample[-2]
    phi_stage = free_sample[-1]

    # Create ForwardContext and cache fixed-stage rotation matrices
    ctx = ForwardContext(geometry)
    free_names = {chi_stage.name, phi_stage.name}
    ctx.prepare_caching(angles, free_names)

    chi_name = chi_stage.name
    phi_name = phi_stage.name

    # --- Analytic fast path for standard Eulerian geometries ---
    if _is_standard_eulerian_pair(chi_stage, phi_stage):
        analytic_results = _solve_bisecting_analytic(
            ctx,
            chi_stage,
            phi_stage,
            angles,
            Q_phi,
        )
        if analytic_results:
            solutions = []
            for chi_d, phi_d in analytic_results:
                sol = dict(angles)
                sol[chi_name] = chi_d
                sol[phi_name] = phi_d
                _apply_cut_points(sol, mode, geometry)

                # De-duplicate
                duplicate = False
                for existing in solutions:
                    if all(abs(existing.get(k, 0) - sol.get(k, 0)) < 1e-4 for k in sol):
                        duplicate = True
                        break
                if not duplicate and _check_limits(geometry, sol):
                    solutions.append(sol)
            if solutions:
                return solutions
        # If analytic solver returned nothing (degenerate case), fall through
        # to the Newton solver below.

    # --- Newton fallback for non-standard axes or degenerate cases ---

    from .orientation import angles_to_phi_vector as _a2phi

    Q_norm = float(np.linalg.norm(Q_phi))
    q_hat = Q_phi / Q_norm

    # Normalised chi axis vector (use pre-normalized from stage)
    chi_ax = chi_stage._axis_hat  # noqa: SLF001

    # Component of q_hat along the chi axis direction
    q_along_chi = float(np.dot(q_hat, chi_ax))
    q_along_chi = max(-1.0, min(1.0, q_along_chi))

    # Analytic seed: chi_abs is the angle between Q_phi and the chi-perp plane
    chi_abs_deg = math.degrees(math.asin(abs(q_along_chi)))

    # In-plane direction for phi seeding
    q_in_plane = q_hat - q_along_chi * chi_ax
    q_in_plane_norm = float(np.linalg.norm(q_in_plane))
    if q_in_plane_norm > 1e-8:
        phi_seed_deg = math.degrees(
            math.atan2(float(q_in_plane[1]), float(q_in_plane[0]))
        )
    else:  # pragma: no cover
        # Degenerate: q parallel to chi axis (q_in_plane ≈ 0).
        phi_seed_deg = 0.0

    # Seed (chi, phi) pairs.  Analytic seeds are tried first (most likely
    # to converge quickly).  Grid seeds are appended as fallback.
    # Early termination: stop after finding enough unique solutions.
    analytic_seeds = []
    for chi_seed in [
        chi_abs_deg,
        -chi_abs_deg,
        180.0 - chi_abs_deg,
        -(180.0 - chi_abs_deg),
    ]:
        for phi_seed in [phi_seed_deg, phi_seed_deg + 180.0]:
            analytic_seeds.append((chi_seed, phi_seed))

    # Coarse 90°-grid seeds to catch solutions far from the analytic guess
    grid_seeds = [
        (chi_s, phi_s)
        for chi_s in (0.0, 90.0, -90.0, 180.0)
        for phi_s in (0.0, 90.0, 180.0, 270.0)
    ]

    seed_pairs = analytic_seeds + grid_seeds

    solutions = []
    _MAX_SOLUTIONS = 4
    # Early termination: after we've found at least 2 solutions, stop
    # after _MAX_STALE consecutive stale seeds.  Before 2 solutions, we
    # try all seeds to ensure both branches are found.
    _MIN_SOLUTIONS = 2
    _MAX_STALE = 6
    stale_count = 0

    for chi_start, phi_start in seed_pairs:
        sol_angles = _solve_two_angles(
            geometry=geometry,
            fixed_angles=angles,
            free_stage_1=chi_name,
            free_stage_2=phi_name,
            start_1=chi_start,
            start_2=phi_start,
            Q_phi_target=Q_phi,
            a2phi_fn=_a2phi,
            _ctx=ctx,
        )
        if sol_angles is None:  # pragma: no branch
            stale_count += 1  # pragma: no cover
            if (
                len(solutions) >= _MIN_SOLUTIONS and stale_count >= _MAX_STALE
            ):  # pragma: no cover
                break  # pragma: no cover
            continue  # pragma: no cover

        sol = dict(angles)
        sol[chi_name] = sol_angles[0]
        sol[phi_name] = sol_angles[1]

        _apply_cut_points(sol, mode, geometry)

        # De-duplicate: skip if this solution is the same as an existing one
        duplicate = False
        for existing in solutions:
            if all(abs(existing.get(k, 0) - sol.get(k, 0)) < 1e-4 for k in sol):
                duplicate = True
                break

        if duplicate or not _check_limits(geometry, sol):
            stale_count += 1
            if len(solutions) >= _MIN_SOLUTIONS and stale_count >= _MAX_STALE:
                break
            continue

        solutions.append(sol)
        stale_count = 0  # reset: we found a new solution
        if len(solutions) >= _MAX_SOLUTIONS:  # pragma: no cover
            break  # pragma: no cover

    return solutions


_CARDINAL_AXES = (
    np.array([1.0, 0.0, 0.0]),
    np.array([-1.0, 0.0, 0.0]),
    np.array([0.0, 1.0, 0.0]),
    np.array([0.0, -1.0, 0.0]),
    np.array([0.0, 0.0, 1.0]),
    np.array([0.0, 0.0, -1.0]),
)


def _is_cardinal_axis(n: np.ndarray, atol: float = 1e-12) -> bool:
    """
    Return True if ``n`` is exactly (within ``atol``) ``±XHAT``, ``±YHAT``,
    or ``±ZHAT``.

    The strict cardinal-axis check gates the analytic 1-free-angle fast
    path (issue #227).  Geometries with axes that are slightly off-axis
    (e.g. mis-aligned mounts, kappa-style tilted axes) will fall back to
    the Newton solver.

    Parameters
    ----------
    n : numpy.ndarray, shape (3,)
        Unit-norm axis vector to test.
    atol : float, optional
        Absolute tolerance for the equality check.  Default ``1e-12``,
        matching the cleanliness of basis vectors that come from
        :data:`~constants.XHAT` / ``YHAT`` / ``ZHAT`` after stage
        construction.

    Returns
    -------
    bool
    """
    return any(np.allclose(n, axis, atol=atol, rtol=0.0) for axis in _CARDINAL_AXES)


def _solve_one_free_angle_analytic(
    ctx: ForwardContext,
    free_stage,
    angles: dict[str, float],
    Q_phi_target: np.ndarray,
) -> float | None:
    r"""
    Analytic single-``atan2`` solver for one free sample-stage rotation.

    Under the BL1967 standard convention (Busing & Levy 1967,
    outermost stage leftmost), when all detector stages and all sample
    stages **except one** have fixed angles, the forward equation
    reduces to a single rotation:

    .. math::

        Z = R_{\text{before}} \cdot R(\hat n, \theta) \cdot R_{\text{after}}

    where :math:`R_{\text{before}}` is the product of fixed sample stages
    *outer* to the free stage (already cached as ``ctx._cached_Z_prefix``)
    and :math:`R_{\text{after}}` is the product of fixed sample stages
    *inner* to the free stage.  The target satisfies
    :math:`Q_\phi = Z^T Q_{\text{lab}}`, so

    .. math::

        R(\hat n, \theta) \cdot q = u

    with :math:`q = R_{\text{after}} \cdot Q_{\phi,\text{target}}` and
    :math:`u = R_{\text{before}}^T \cdot Q_{\text{lab}}`.  The unknown
    rotation about a known axis :math:`\hat n` is recovered by projecting
    both vectors onto the plane perpendicular to :math:`\hat n` and
    reading off the angle with a single ``atan2``.

    Parameters
    ----------
    ctx : ForwardContext
        Must have :meth:`~ForwardContext.prepare_caching` already called
        with ``free_stage.name`` as the only free name.
    free_stage : Stage
        The single free sample stage.  ``free_stage._axis_hat`` must be
        a unit vector.
    angles : dict[str, float]
        Baseline angles with all fixed stages set.
    Q_phi_target : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (Å⁻¹).

    Returns
    -------
    float or None
        Angle in degrees that satisfies the rotation equation, or
        ``None`` when:

        * the free stage's axis is not the *first* free index (i.e. there
          are non-trivial fixed sample stages after it that
          :class:`ForwardContext` does not currently cache as a suffix),
        * the projections onto the plane perpendicular to ``n`` are
          degenerate (``q`` or ``u`` is parallel to ``n``), or
        * the parallel-component consistency check fails (``Q_phi_target``
          is not reachable by varying the free angle alone).

        In all such cases the caller falls back to Newton.
    """
    from .rotation import _rotation_matrix_normalized

    # ForwardContext currently caches only the prefix (stages strictly
    # before the first free stage).  Stages after the free one are not
    # cached, so build R_after on the fly here.
    sample_stages = ctx.sample_stages
    free_idx = next(
        (i for i, s in enumerate(sample_stages) if s.name == free_stage.name),
        None,
    )
    if free_idx is None:  # pragma: no cover
        return None

    # Build R_after in BL1967 standard outermost-leftmost order:
    # R_after = R_{free_idx+1} @ R_{free_idx+2} @ ... @ R_{N-1}.
    R_after = np.eye(3)
    for i in range(free_idx + 1, len(sample_stages)):
        s = sample_stages[i]
        a = angles.get(s.name, s.angle)
        R_after = R_after @ _rotation_matrix_normalized(s._axis_hat, a)  # noqa: SLF001

    Z_prefix = ctx._cached_Z_prefix if ctx._cached_Z_prefix is not None else np.eye(3)
    D = ctx._cached_D if ctx._cached_D is not None else np.eye(3)
    Q_lab = ctx.two_pi_over_lambda * (D @ ctx.y_eff - ctx.y_eff)

    # Solve R(n, θ) q = u with
    #   q = R_after @ Q_phi_target
    #   u = Z_prefix^T @ Q_lab
    # (Derivation: Z = Z_prefix @ R(n, θ) @ R_after, so
    #  Q_phi = Z^T Q_lab = R_after^T R(n, θ)^T Z_prefix^T Q_lab; multiply
    #  through by R(n, θ) R_after on the left to get the equation above.)
    n = free_stage._axis_hat  # noqa: SLF001
    q = R_after @ np.asarray(Q_phi_target, dtype=float)
    u = Z_prefix.T @ Q_lab

    # Parallel components must agree: rotation about n preserves the
    # n-component.  If they disagree beyond tolerance the target is not
    # reachable by varying only the free angle.
    q_par = float(np.dot(q, n))
    u_par = float(np.dot(u, n))
    if abs(q_par - u_par) > 1e-7:
        return None

    q_perp = q - q_par * n
    u_perp = u - u_par * n
    q_perp_norm = float(np.linalg.norm(q_perp))
    if q_perp_norm < 1e-12:
        return None  # degenerate: q ∥ n  (any θ rotates q to itself)

    # Build orthonormal basis (e1, e2) in the plane ⟂ n.
    e1 = q_perp / q_perp_norm
    e2 = np.cross(n, e1)

    # Match magnitudes: rotation preserves |q_perp| = |u_perp|.
    u_perp_norm = float(np.linalg.norm(u_perp))
    if abs(q_perp_norm - u_perp_norm) > 1e-7:  # pragma: no cover
        # When the parallel-component check above passes (q_par == u_par)
        # and the targets are derived from a consistent geometry/UB/wavelength
        # configuration (the only call path in the package), rotation-
        # preservation of perpendicular magnitudes is automatic; this branch
        # is a defensive guard for caller-constructed contexts only.
        return None  # not reachable by rotation alone

    cos_t = float(np.dot(u_perp, e1))
    sin_t = float(np.dot(u_perp, e2))
    return math.degrees(math.atan2(sin_t, cos_t))


def _solve_one_free_angle(
    geometry,
    fixed_angles: dict,
    free_stage,
    Q_phi_target: np.ndarray,
    mode,
) -> list[dict[str, float]]:
    """
    Numerically solve for one free stage angle such that
    ``angles_to_phi_vector(geometry, **all_angles) == Q_phi_target``.

    Used when only one sample stage is free (e.g. phi in fixed_chi mode).

    **Fast path (issue #227).**  When the free stage's axis vector is
    exactly one of the cardinal directions (``±XHAT``, ``±YHAT``,
    ``±ZHAT``), :func:`_solve_one_free_angle_analytic` derives the angle
    via a single ``atan2`` call.  This eliminates Newton iteration
    entirely (no Jacobian, no seed sweep) and is 10-20× faster than the
    Newton fallback.

    For non-cardinal axes the function falls back to the previous Newton
    solver, which seeds from four quadrants and deduplicates.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    fixed_angles : dict
        Baseline angles (including all constrained stages).
    free_stage : Stage
        The single free stage.
    Q_phi_target : numpy.ndarray, shape (3,)
    mode : DiffractionMode

    Returns
    -------
    list of dict[str, float]
    """
    from .orientation import angles_to_phi_vector as _a2phi

    # Create ForwardContext and cache fixed-stage rotation matrices
    ctx = ForwardContext(geometry)
    ctx.prepare_caching(fixed_angles, {free_stage.name})

    solutions: list[dict[str, float]] = []

    # --- Analytic fast path (issue #227) ---
    if _is_cardinal_axis(free_stage._axis_hat):  # noqa: SLF001
        theta = _solve_one_free_angle_analytic(
            ctx, free_stage, fixed_angles, Q_phi_target
        )
        if theta is not None:
            # Validate the analytic solution against the actual residual.
            # The atan2 result is exact in exact arithmetic, so this check
            # only fails for pathological numerical edge cases.
            trial = dict(fixed_angles)
            trial[free_stage.name] = theta
            Q_check = ctx.q_phi(trial)
            if (
                float(np.linalg.norm(Q_check - Q_phi_target)) < 1e-7
            ):  # pragma: no branch
                # Normalise to (-180, 180]
                theta_n = (theta + 180.0) % 360.0 - 180.0
                sol = dict(fixed_angles)
                sol[free_stage.name] = theta_n
                _apply_cut_points(sol, mode, geometry)
                if _check_limits(geometry, sol):  # pragma: no branch
                    solutions.append(sol)
                return solutions
        # Fall through to Newton if the fast path could not validate.

    # --- Newton fallback (non-cardinal axes or degenerate fast path) ---
    for start in [0.0, 90.0, -90.0, 180.0]:
        result = _solve_two_angles(
            geometry=geometry,
            fixed_angles=fixed_angles,
            free_stage_1=free_stage.name,
            free_stage_2=free_stage.name,  # dummy — same stage
            start_1=start,
            start_2=start,
            Q_phi_target=Q_phi_target,
            a2phi_fn=_a2phi,
            _one_dimensional=True,
            _ctx=ctx,
        )
        if result is None:  # pragma: no branch
            continue  # pragma: no cover
        sol = dict(fixed_angles)
        sol[free_stage.name] = result[0]
        _apply_cut_points(sol, mode, geometry)
        duplicate = any(
            abs(existing.get(free_stage.name, 0) - sol[free_stage.name]) < 1e-4
            for existing in solutions
        )
        if not duplicate and _check_limits(geometry, sol):
            solutions.append(sol)
    return solutions


def _solve_two_angles(
    geometry,
    fixed_angles: dict,
    free_stage_1: str,
    free_stage_2: str,
    start_1: float,
    start_2: float,
    Q_phi_target: np.ndarray,
    a2phi_fn,
    max_iter: int = 200,
    tol: float = 1e-10,
    _one_dimensional: bool = False,
    _ctx: ForwardContext | None = None,
) -> tuple[float, float] | None:
    """
    Numerically solve for one or two free stage angles such that
    ``angles_to_phi_vector(geometry, **all_angles) == Q_phi_target``.

    Uses Gauss-Newton iteration on the 3-component residual.  For the
    2D case (default) the system is 3×2 (overdetermined but consistent).
    For the 1D case (``_one_dimensional=True``), only ``free_stage_1`` is
    varied; ``free_stage_2`` is ignored.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    fixed_angles : dict
        All stage angles including the constrained ones.  This dict is
        used as the base; free stages are overwritten.
    free_stage_1, free_stage_2 : str
        Names of the free stages.  When ``_one_dimensional=True``,
        only ``free_stage_1`` is varied.
    start_1, start_2 : float
        Starting guess in degrees.
    Q_phi_target : numpy.ndarray, shape (3,)
    a2phi_fn : callable
        ``angles_to_phi_vector`` function (used as fallback if no context).
    max_iter, tol : int, float
        Convergence parameters.
    _one_dimensional : bool
        If True, solve for only ``free_stage_1`` (1D problem).
    _ctx : ForwardContext or None
        Pre-computed context for fast Q_phi evaluation.  When provided,
        bypasses the stateful ``a2phi_fn`` and uses cached rotation
        matrices instead.

    Returns
    -------
    (angle_1, angle_2) : tuple of float, or None if no convergence.
        When ``_one_dimensional=True``, ``angle_2`` equals ``start_2``.
    """
    n_free = 1 if _one_dimensional else 2
    x = np.array([start_1, start_2][:n_free], dtype=float)

    # Use a mutable trial dict that is reused across iterations to avoid
    # allocating a new dict on every residual call.
    trial = dict(fixed_angles)

    if _ctx is not None:

        def residual(x):
            trial[free_stage_1] = float(x[0])
            if not _one_dimensional:
                trial[free_stage_2] = float(x[1])
            return _ctx.q_phi(trial) - Q_phi_target
    else:

        def residual(x):
            trial[free_stage_1] = float(x[0])
            if not _one_dimensional:
                trial[free_stage_2] = float(x[1])
            Q = a2phi_fn(geometry, **trial)
            return Q - Q_phi_target

    def jacobian_fd(x, r0, h=1e-4):  # pragma: no cover
        """Finite-difference Jacobian (3 × n_free) — fallback when no ctx."""
        J = np.zeros((3, n_free))  # pragma: no cover
        for i in range(n_free):  # pragma: no cover
            xp = x.copy()  # pragma: no cover
            xp[i] += h  # pragma: no cover
            J[:, i] = (residual(xp) - r0) / h  # pragma: no cover
        return J  # pragma: no cover

    # Build the ordered list of free stage names for the analytic Jacobian.
    if _ctx is not None:
        if _one_dimensional:
            _free_names = [free_stage_1]
        else:
            _free_names = [free_stage_1, free_stage_2]

    for _ in range(max_iter):
        r = residual(x)
        if np.linalg.norm(r) < tol:
            a1 = float(x[0])
            a2 = float(x[1]) if not _one_dimensional else start_2
            return (a1, a2)
        if _ctx is not None:
            J = _ctx.jacobian_analytic(trial, _free_names)
        else:
            J = jacobian_fd(x, r)  # pragma: no cover
        try:
            dx, _, _, _ = np.linalg.lstsq(J, -r, rcond=None)
        except np.linalg.LinAlgError:  # pragma: no cover
            return None  # pragma: no cover
        step = np.linalg.norm(dx)
        if step > 30.0:
            dx = dx * (30.0 / step)
        x = x + dx

    r = residual(x)  # pragma: no cover
    if np.linalg.norm(r) < 1e-6:  # pragma: no cover
        a1 = float(x[0])  # pragma: no cover
        a2 = float(x[1]) if not _one_dimensional else start_2  # pragma: no cover
        return (a1, a2)  # pragma: no cover
    return None  # pragma: no cover


def _is_kappa_virtual_mode(
    geometry: AdHocDiffractometer,
    mode,
) -> bool:
    """
    Return True when the mode contains a virtual Eulerian angle constraint
    (omega, chi, or phi) on a kappa geometry.

    Delegates to :func:`~kappa.is_kappa_virtual_mode`.
    """
    from .kappa import is_kappa_virtual_mode

    return is_kappa_virtual_mode(geometry, mode)


def _solve_kappa_virtual(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Kappa virtual-angle dispatcher.

    Wraps :func:`~kappa.solve_kappa_virtual` (which performs the
    geometry-aware Eulerian-equivalent solve introduced by issue #241)
    with the kappa-side post-processing: virtual-angle validation
    against the per-geometry
    :class:`~kappa.KappaPseudoAngleConvention`, cut-points, motor
    limits, and deduplication.
    """
    from .kappa import KAPPA_VIRTUAL_ANGLES
    from .kappa import kappa_to_eulerian_axes
    from .kappa import solve_kappa_virtual
    from .mode import SampleConstraint

    raw = solve_kappa_virtual(geometry, Q_phi, ttheta_deg, mode)
    if not raw:  # pragma: no cover
        return []

    convention = geometry.kappa_pseudo_angle_convention

    sample_stages = geometry.sample_stages
    kappa_idx = next(
        (i for i, s in enumerate(sample_stages) if s.name == "kappa"), None
    )
    komega_name = sample_stages[kappa_idx - 1].name if kappa_idx else None
    kphi_name = sample_stages[kappa_idx + 1].name if kappa_idx else None

    virtual_constraints = [
        c
        for c in mode.constraints
        if isinstance(c, SampleConstraint) and c.name in KAPPA_VIRTUAL_ANGLES
    ]

    solutions: list[dict[str, float]] = []
    for angle_dict in raw:
        angles = dict(angle_dict)

        # Wrap kappa angles into (-180, 180] (standard motor cut).
        if komega_name and kphi_name:  # pragma: no branch
            for kname in (komega_name, "kappa", kphi_name):
                if kname in angles:  # pragma: no branch
                    a = angles[kname]
                    a = (a + 180.0) % 360.0 - 180.0
                    if a == -180.0:  # pragma: no cover
                        a = 180.0
                    angles[kname] = a

        # Validate any fixed virtual-angle constraints by inverting
        # the (komega, kappa, kphi) triple via the geometry-aware
        # decomposition.
        if virtual_constraints and convention is not None:  # pragma: no branch
            try:
                om_v, chi_v, phi_v = kappa_to_eulerian_axes(
                    angles[komega_name],
                    angles["kappa"],
                    angles[kphi_name],
                    convention,
                )
                virtual_vals = {"omega": om_v, "chi": chi_v, "phi": phi_v}
                if not all(
                    abs(virtual_vals[c.name] - float(c.value)) < 1e-4
                    for c in virtual_constraints
                ):
                    continue  # pragma: no cover
            except (ValueError, KeyError):  # pragma: no cover
                continue

        # Verify geometric consistency: Q_computed must match Q_target.
        kv_ctx = ForwardContext(geometry)
        Q_computed = kv_ctx.q_phi_uncached(angles)
        if not np.allclose(Q_computed, Q_phi, atol=1e-6):  # pragma: no cover
            continue

        _apply_cut_points(angles, mode, geometry)
        if not _check_limits(geometry, angles):
            continue  # pragma: no cover

        duplicate = False
        for existing in solutions:
            if all(
                abs(existing.get(kk, 0) - angles.get(kk, 0)) < 1e-4 for kk in angles
            ):  # pragma: no cover
                duplicate = True
                break
        if not duplicate:  # pragma: no branch
            solutions.append(angles)

    return solutions


# ---------------------------------------------------------------------------
# Psi-constant solver (ReferenceConstraint("psi") validation filter)
# ---------------------------------------------------------------------------


def _compute_natural_psi(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
) -> float | None:
    """
    Compute the natural azimuthal angle ψ from the phi frame only.

    ψ is a pure phi-frame quantity: it depends on ``Q_phi = UB @ hkl``,
    ``n_phi = UB @ n_hkl`` (the azimuthal reference in the phi frame),
    and ``y_hat`` (the incident-beam direction from the geometry's basis).
    **No motor angles are involved.**

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``sample.UB`` and ``azimuth`` set.
    Q_phi : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (``UB @ hkl``).

    Returns
    -------
    float or None
        ψ in degrees (−180°, +180°], or ``None`` when ψ is undefined
        (Q ∥ incident beam, or reference vector ∥ Q).
    """
    n_hkl = geometry.azimuth
    if n_hkl is None:
        return None  # pragma: no cover

    n_phi = geometry.sample.UB @ np.asarray(n_hkl, dtype=float)
    y_raw = np.asarray(
        geometry.basis.get("longitudinal", np.array([0.0, 1.0, 0.0])),
        dtype=float,
    )
    y_hat = y_raw / np.linalg.norm(y_raw)

    Q_mag = float(np.linalg.norm(Q_phi))
    if Q_mag < 1e-14:
        return None  # pragma: no cover
    Q_hat = Q_phi / Q_mag

    # Project n and y onto the plane perpendicular to Q
    n_perp = n_phi - np.dot(n_phi, Q_hat) * Q_hat
    y_perp = y_hat - np.dot(y_hat, Q_hat) * Q_hat

    n_perp_mag = float(np.linalg.norm(n_perp))
    y_perp_mag = float(np.linalg.norm(y_perp))

    if n_perp_mag < 1e-10 or y_perp_mag < 1e-10:
        return None  # ψ undefined (reference ∥ Q, or Q ∥ incident beam)

    n_perp_hat = n_perp / n_perp_mag
    y_perp_hat = y_perp / y_perp_mag

    cos_psi = float(np.clip(np.dot(y_perp_hat, n_perp_hat), -1.0, 1.0))
    sin_psi = float(np.dot(Q_hat, np.cross(y_perp_hat, n_perp_hat)))
    return math.degrees(math.atan2(sin_psi, cos_psi))


def _is_psi_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when the mode contains a psi ReferenceConstraint and the reference is set."""
    from .mode import ReferenceConstraint

    if geometry.azimuth is None:
        return False
    return any(
        isinstance(c, ReferenceConstraint) and c.name == "psi" for c in mode.constraints
    )


def _solve_psi_mode(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Forward solver for psi-constant modes (validation filter).

    For a given (h,k,l) and UB, the azimuthal angle ψ is a pure phi-frame
    quantity — the same for every Bragg solution.  This solver:

    1. Computes the natural ψ from ``Q_phi`` (no motor angles).
    2. Compares with ``psi_target`` from the mode's
       :class:`~mode.ReferenceConstraint`.
    3. If ψ is undefined for the reflection (Q ∥ azimuthal reference, or
       Q ∥ incident beam): emits a :class:`UserWarning` and returns ``[]``.
    4. If natural ψ and target disagree beyond 0.1°: emits a
       :class:`UserWarning` naming the natural ψ value and returns ``[]``.
    5. If they agree: delegates to the appropriate existing solver
       (bisecting, kappa-virtual, or synthetic bisecting) and returns
       all solutions.

    The warnings (added in issue #278) make the "ψ is fixed by UB and
    hkl, not by motors" semantics visible to callers who previously saw
    only a silent empty list.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
    """
    import warnings

    from .mode import BisectConstraint
    from .mode import ConstraintSet
    from .mode import ReferenceConstraint

    # Extract the psi target value
    rc = next(
        c
        for c in mode.constraints
        if isinstance(c, ReferenceConstraint) and c.name == "psi"
    )
    psi_target = float(rc.value)

    # Compute natural psi from the phi frame (motor-angle independent)
    natural_psi = _compute_natural_psi(geometry, Q_phi)
    if natural_psi is None:
        warnings.warn(
            f"forward(): ψ is undefined for this reflection in geometry "
            f"{geometry.name!r} — Q is parallel to "
            f"azimuth={geometry.azimuth} (or to the "
            f"incident beam).  Choose a different reflection or change "
            f"geometry.azimuth.  Returning [].",
            UserWarning,
            stacklevel=5,
        )
        return []  # ψ undefined for this reflection

    # Compare natural psi with target (tolerance 0.1° — generous enough
    # to handle float rounding, tight enough to be physically meaningful)
    diff = abs(natural_psi - psi_target)
    # Handle wraparound: e.g. -179.9 vs 180.1
    if diff > 180.0:
        diff = 360.0 - diff
    if diff > 0.1:
        warnings.warn(
            f"forward(): mode {geometry.mode_name!r} targets ψ = "
            f"{psi_target:.4f}° but the natural ψ for this reflection is "
            f"{natural_psi:.4f}°.  ψ is fixed by UB and (h, k, l); no motor "
            f"configuration can change it for a given reflection.  Either "
            f"set the constraint value to {natural_psi:.4f}° (or use "
            f"ad_hoc_diffractometer.reference.natural_psi(g, h, k, l) to "
            f"discover the natural value), or pick a different reflection.  "
            f"Returning [].",
            UserWarning,
            stacklevel=5,
        )
        return []  # this (h,k,l) is not accessible at the stored ψ

    # ψ is satisfied — delegate to the appropriate existing solver.
    # The psi constraint is automatically satisfied by ALL Bragg solutions.

    if mode.has_bisect:
        # kappa6c: mode already has a BisectConstraint
        return _solve_bisecting(geometry, Q_phi, ttheta_deg, mode)

    if _is_kappa_virtual_mode(geometry, mode):  # pragma: no cover
        # kappa4cv, kappa4ch: kappa virtual angle mode
        return _solve_kappa_virtual(geometry, Q_phi, ttheta_deg, mode)

    # psic-family fixed_psi_* modes (issue #264 C1/C2 revision): the
    # bisect was dropped in favor of a SampleConstraint + a
    # DetectorConstraint pinning the scattering plane (nu = 0 for
    # vertical, delta = 0 for horizontal) plus the psi reference.  Once
    # ψ is validated, delegate to ``_solve_fixed_sample`` which handles
    # the remaining free sample stages plus the active detector at
    # ttheta.
    if (
        any(s.name == "chi" for s in geometry.sample_stages)
        and mode.detector_constraint is not None
        and len(mode.fixed_sample_constraints) >= 1
    ):
        stripped = ConstraintSet(
            [
                c
                for c in mode.constraints
                if not (isinstance(c, ReferenceConstraint) and c.name == "psi")
            ],
            computed=mode.computed,
            extras=mode.extras,
            cut_points=mode.cut_points,
        )
        return _solve_fixed_sample(geometry, Q_phi, ttheta_deg, stripped)

    # fourcv, fourch, kappa4cv, kappa4ch: no BisectConstraint in the mode.
    # Build a synthetic bisecting mode using the geometry's natural
    # bisect pair: first sample stage + first detector stage.
    sample_stages = geometry.sample_stages
    detector_stages = geometry.detector_stages

    if not sample_stages or not detector_stages:
        return []  # pragma: no cover

    bisect_sample = sample_stages[0].name  # omega in fourcv/fourch
    bisect_detector = detector_stages[-1].name  # ttheta in fourcv/fourch

    # Collect any non-psi constraints from the original mode
    other_constraints = [
        c
        for c in mode.constraints
        if not (isinstance(c, ReferenceConstraint) and c.name == "psi")
    ]

    synthetic = ConstraintSet(
        [BisectConstraint(bisect_sample, bisect_detector)] + other_constraints,
        computed=mode.computed,
    )

    return _solve_bisecting(geometry, Q_phi, ttheta_deg, synthetic)


# ---------------------------------------------------------------------------
# Double-diffraction solver (simultaneous Bragg for two reflections)
# ---------------------------------------------------------------------------


def _is_double_diffraction_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when mode.extras contains h2, k2, l2 keys."""
    extras = mode.extras
    return extras is not None and all(k in extras for k in ("h2", "k2", "l2"))


def _solve_double_diffraction(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Decomposed solver for double-diffraction modes.

    Finds motor angles where both the primary reflection (h1,k1,l1) and a
    secondary reflection (h2,k2,l2) simultaneously satisfy the Ewald sphere
    condition.

    Instead of a 4D Newton-Raphson system, the problem is decomposed into
    sequential subproblems:

    1. **Detector angle** -- known from Bragg's law (``ttheta_deg``).
    2. **Outer sample angle scan** -- 1D sweep over the outermost free
       sample stage (omega, eta, mu, or komega).
    3. **(chi, phi) from Q direction** -- for each outer angle, solve the
       primary Bragg condition for the two inner sample stages.  Uses the
       analytic decomposition (``_solve_bisecting_analytic``) for standard
       Eulerian geometries, or the 2D Gauss-Newton solver for kappa axes.
    4. **Ewald sphere filter** -- accept only candidates where the secondary
       reflection hkl2 lies on the Ewald sphere.  A 1D Newton refinement
       on the outer sample angle drives the residual to machine precision.

    This replaces the former 4D Newton system (~50 seeds x 300 iterations x
    5 residual evaluations per step) with a 1D scan + analytic 2D solve +
    scalar filter, yielding a 100-1000x speedup.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]

    Raises
    ------
    ValueError
        If h2, k2, l2 are not set in ``mode.extras`` (still ``REQUIRED``).
    """
    from .mode import REQUIRED
    from .rotation import rotation_matrix

    extras = mode.extras

    # Validate that h2, k2, l2 are set to numeric values
    hkl2_values = [extras.get("h2"), extras.get("k2"), extras.get("l2")]
    if any(v is REQUIRED for v in hkl2_values):
        raise ValueError(
            "double_diffraction mode requires h2, k2, l2 to be set in "
            "mode.extras before calling forward(). "
            "Set them with e.g. "
            "g.modes['double_diffraction'].extras['h2'] = 1.0"
        )

    hkl2 = np.array([float(v) for v in hkl2_values], dtype=float)
    Q2_phi = geometry.sample.UB @ hkl2

    # Incident beam: ki = (2pi/lambda) * longitudinal_hat
    wavelength = geometry.wavelength
    k_mag = 2.0 * math.pi / wavelength
    y_raw = np.asarray(
        geometry.basis.get("longitudinal", np.array([0.0, 1.0, 0.0])),
        dtype=float,
    )
    ki = k_mag * (y_raw / np.linalg.norm(y_raw))
    ki_sq = float(np.dot(ki, ki))

    # Build baseline angles: all stages at their current positions
    all_stages = list(geometry._stages.values())  # noqa: SLF001
    angles_base: dict[str, float] = {s.name: s.angle for s in all_stages}

    # Apply fixed sample constraints (mu=0, eta=0, etc.)
    for sc in mode.fixed_sample_constraints:
        if sc.name in geometry._stages:  # noqa: SLF001  # pragma: no branch
            angles_base[sc.name] = float(sc.value)

    # Apply fixed detector constraints (nu=0, delta=0, etc.)
    det_c = mode.detector_constraint
    if det_c is not None and not det_c.is_qaz:
        angles_base[det_c.name] = det_c.value

    # The 4 free stages: always [3 sample, 1 detector]
    free_names = list(mode.computed)

    if len(free_names) != 4:  # pragma: no cover
        logger.warning(
            "double_diffraction expects 4 free stages, got %d: %s",
            len(free_names),
            free_names,
        )
        return []

    # Classify the 4 free stages into sample and detector
    sample_names = {s.name for s in geometry.sample_stages}
    det_names = {s.name for s in geometry.detector_stages}
    free_sample = [n for n in free_names if n in sample_names]
    free_det = [n for n in free_names if n in det_names]

    if len(free_sample) != 3 or len(free_det) != 1:  # pragma: no cover
        logger.warning(
            "double_diffraction expects 3 sample + 1 detector free stages, "
            "got %d sample + %d detector: %s",
            len(free_sample),
            len(free_det),
            free_names,
        )
        return []

    # Step 1: Detector angle is known from Bragg's law
    det_stage_name = free_det[0]
    angles_base[det_stage_name] = ttheta_deg

    # The 3 free sample stages in stacking order:
    #   outer_stage (omega/eta/mu/komega) -- scanned in 1D
    #   chi_stage (chi/kappa) -- solved from Q direction
    #   phi_stage (phi/kphi) -- solved from Q direction
    sample_stages = geometry.sample_stages
    outer_stage_name = free_sample[0]
    chi_stage_name = free_sample[1]
    phi_stage_name = free_sample[2]

    # Look up the Stage objects for the two inner stages
    chi_stage = geometry.stage(chi_stage_name)
    phi_stage = geometry.stage(phi_stage_name)

    # Check if the inner pair is standard Eulerian (orthogonal axes)
    is_eulerian = _is_standard_eulerian_pair(chi_stage, phi_stage)

    def _build_Z(angles: dict[str, float]) -> np.ndarray:
        """Compute sample rotation matrix from angle values.

        Composition follows the BL1967 standard convention (Busing &
        Levy 1967): outermost (floor-most) stage leftmost,
        so Z = R_0 @ R_1 @ ... @ R_{N-1}.
        """
        Z = np.eye(3)
        for s in sample_stages:
            Z = Z @ rotation_matrix(s.axis, angles[s.name])
        return Z

    def _ewald_residual(angles: dict[str, float]) -> float:
        """Scalar Ewald sphere residual for secondary reflection."""
        Z = _build_Z(angles)
        Q2_lab = Z @ Q2_phi
        kf2 = ki + Q2_lab
        return float(np.dot(kf2, kf2)) - ki_sq

    # Pre-compute quantities for the fast Eulerian inner solve.
    # For Eulerian geometries, we inline the analytic decomposition to avoid
    # creating a new ForwardContext per call.
    if is_eulerian:  # pragma: no branch — kappa tested by slow_benchmark
        from .rotation import _rotation_matrix_normalized

        n_chi = chi_stage._axis_hat  # noqa: SLF001
        n_phi = phi_stage._axis_hat  # noqa: SLF001
        n3 = np.cross(n_chi, n_phi)

        # Pre-compute the outer-stage-independent pieces:
        # The detector rotation D is fixed (ttheta is set).
        # Z_prefix depends on the outer angle and any fixed stages before it.
        # We need: v = Q_lab = (2pi/lambda) * (D @ y_eff - y_eff)
        # and q = Z_prefix(outer) @ Q_phi_target
        #
        # Build D and y_eff once.
        _dd_ctx_base = ForwardContext(geometry)
        _dd_ctx_base.prepare_caching(
            angles_base, {outer_stage_name, chi_stage_name, phi_stage_name}
        )
        D = _dd_ctx_base._cached_D if _dd_ctx_base._cached_D is not None else np.eye(3)
        y_eff = _dd_ctx_base.y_eff
        two_pi_over_lambda = _dd_ctx_base.two_pi_over_lambda
        v = two_pi_over_lambda * (D @ y_eff - y_eff)  # Q_lab (constant)

        # Find which sample stages are in the Z_prefix (before the outer stage)
        outer_stage_idx = next(
            i for i, s in enumerate(sample_stages) if s.name == outer_stage_name
        )
        outer_stage_obj = sample_stages[outer_stage_idx]

        # Pre-compute the Z_prefix for stages before the outer stage.
        # Composition: BL1967 standard outermost-leftmost, so
        # Z_pre_outer = R_0 @ R_1 @ ... @ R_{outer_stage_idx-1}.
        Z_pre_outer = np.eye(3)
        for i in range(outer_stage_idx):
            s = sample_stages[i]
            angle = angles_base.get(s.name, s.angle)
            Z_pre_outer = Z_pre_outer @ _rotation_matrix_normalized(s._axis_hat, angle)  # noqa: SLF001

    def _solve_inner_eulerian_fast(
        outer_deg: float,
    ) -> list[tuple[float, float]]:
        """
        Fast analytic (chi, phi) solve for Eulerian geometries.

        Inlines the algebra of :func:`_solve_bisecting_analytic` for the
        double-diffraction outer-scan loop.  Under the BL1967 standard
        convention (Busing & Levy 1967, outermost-leftmost):

            Z_prefix = Z_pre_outer @ R_outer
            R_chi R_phi v = q
              with  v = Q_phi (target),  q = Z_prefix^T @ Q_lab

        Step 1 — solve for phi (two candidates) by projecting onto
        n_chi (invariant under R_chi).  Step 2 — for each phi, recover
        chi as the residual rotation about n_chi taking
        w = R(n_phi, phi) v onto q.

        Returns list of (chi_deg, phi_deg) pairs.  Returns an empty list
        when v ∥ n_phi (phi indeterminate); the caller falls back to a
        chi-scan via ``_find_degenerate_outers`` /
        ``_solve_degenerate_outer``.
        """
        R_outer = _rotation_matrix_normalized(
            outer_stage_obj._axis_hat,
            outer_deg,  # noqa: SLF001
        )
        # Textbook order: Z_prefix = Z_pre_outer @ R_outer.
        Z_prefix = Z_pre_outer @ R_outer

        # q_loc = Z_prefix^T @ Q_lab  (constant within an outer-angle step)
        # v_target = Q_phi_target     (constant)
        v_target = Q_phi  # alias for clarity
        q_loc = Z_prefix.T @ v  # v is the Q_lab vector defined at line 2064

        # --- Step 1: A cos phi + B sin phi = C ---
        A = float(np.dot(n_chi, v_target))
        B = float(np.dot(n3, v_target))
        C = float(np.dot(n_chi, q_loc))

        amp = math.sqrt(A * A + B * B)
        if amp < 1e-12:
            return []  # v ∥ n_phi; defer to degenerate path

        cos_arg_phi = C / amp
        if abs(cos_arg_phi) > 1.0 + 1e-8:
            return []
        cos_arg_phi = max(-1.0, min(1.0, cos_arg_phi))

        phi0 = math.atan2(B, A)
        delta_phi = math.acos(cos_arg_phi)
        phi_candidates_rad = [phi0 + delta_phi, phi0 - delta_phi]

        # --- Step 2: for each phi, recover chi from the residual ---
        results = []
        for phi_rad in phi_candidates_rad:
            phi_d = math.degrees(phi_rad)
            c_ph = math.cos(phi_rad)
            s_ph = math.sin(phi_rad)
            nphi_v_target = float(np.dot(n_phi, v_target))
            w = (
                v_target * c_ph
                + np.cross(n_phi, v_target) * s_ph
                + n_phi * nphi_v_target * (1.0 - c_ph)
            )

            w_par = float(np.dot(w, n_chi))
            q_par = float(np.dot(q_loc, n_chi))
            w_perp = w - w_par * n_chi
            q_perp = q_loc - q_par * n_chi
            w_perp_norm = float(np.linalg.norm(w_perp))
            q_perp_norm = float(np.linalg.norm(q_perp))
            if w_perp_norm < 1e-12 or q_perp_norm < 1e-12:  # pragma: no cover
                chi_d = 0.0
            else:
                cos_chi = float(np.dot(w_perp, q_perp)) / (w_perp_norm * q_perp_norm)
                sin_chi = float(np.dot(n_chi, np.cross(w_perp, q_perp))) / (
                    w_perp_norm * q_perp_norm
                )
                chi_d = math.degrees(math.atan2(sin_chi, cos_chi))

            # Normalize to (-180, 180]
            chi_d = (chi_d + 180.0) % 360.0 - 180.0
            phi_d = (phi_d + 180.0) % 360.0 - 180.0
            results.append((chi_d, phi_d))

        return results

    def _find_degenerate_outers() -> list[float]:
        """
        Find outer angles where the analytic decomposition is degenerate.

        Under the BL1967 standard convention (Busing & Levy 1967) the
        degeneracy condition is ``Q_phi ∥ n_phi`` (the *target*
        scattering vector lies along the innermost stage axis), which
        is **independent of the outer angle**.  Either every outer
        angle is degenerate or none is.  Return a sweep of candidate
        outer angles in the degenerate case; otherwise return an empty
        list.

        The sweep includes both a dense scan and the bisecting outer
        angles ``± ttheta_deg / 2``.  The bisecting values are
        appended because they are exactly the outer angles at which
        the primary Bragg condition is solvable in the degenerate
        case (every other outer angle gives ``q_par ≠ w_par`` in the
        inner ``_chi_to_trial`` and is silently rejected).
        """
        A = float(np.dot(n_chi, Q_phi))
        B = float(np.dot(n3, Q_phi))
        if math.sqrt(A * A + B * B) >= 0.1:
            return []
        # Degenerate: phi is indeterminate at every outer angle.  Include
        # the bisecting outer values so the caller can find at least the
        # bisecting-locus solutions, plus a dense sweep elsewhere.
        outers = [-180.0 + i * 0.5 for i in range(720)]
        outers.extend([ttheta_deg / 2.0, -ttheta_deg / 2.0])
        return outers

    def _solve_degenerate_outer(
        outer_deg: float,
    ) -> list[dict[str, float]]:
        """
        At degenerate outer angles (Q_phi ∥ n_phi), scan phi to find
        solutions that satisfy both the Bragg and Ewald conditions.

        Under the corrected outermost-leftmost composition, the
        degeneracy is phi-indeterminate (R(n_phi, phi) @ Q_phi = Q_phi
        for every phi).  We scan phi, compute chi from the residual
        rotation R(n_chi, chi) @ w = q where w = R(n_phi, phi) @ v, and
        detect Ewald sign changes between adjacent phi values.  A
        bisection then refines the phi to the exact root.
        """
        angles = dict(angles_base)
        angles[outer_stage_name] = outer_deg

        ctx = ForwardContext(geometry)
        ctx.prepare_caching(angles, {chi_stage_name, phi_stage_name})

        # Textbook outermost-leftmost composition:
        #   Z_prefix = Z_pre_outer @ R_outer
        R_outer = _rotation_matrix_normalized(
            outer_stage_obj._axis_hat,
            outer_deg,  # noqa: SLF001
        )
        Z_prefix_local = Z_pre_outer @ R_outer
        q_loc = Z_prefix_local.T @ v  # v here is the Q_lab vector
        v_target = Q_phi  # for clarity
        nphi_v = float(np.dot(n_phi, v_target))

        def _chi_to_trial(phi_deg_f: float) -> dict[str, float] | None:
            """Given a phi value, compute chi and return a trial dict."""
            phi_rad = math.radians(phi_deg_f)
            c_ph = math.cos(phi_rad)
            s_ph = math.sin(phi_rad)
            # w = R(n_phi, phi) @ v_target (Rodrigues)
            w = (
                v_target * c_ph
                + np.cross(n_phi, v_target) * s_ph
                + n_phi * nphi_v * (1.0 - c_ph)
            )
            # chi = angle about n_chi from w_perp to q_perp
            w_par = float(np.dot(w, n_chi))
            q_par = float(np.dot(q_loc, n_chi))
            w_perp = w - w_par * n_chi
            q_perp = q_loc - q_par * n_chi
            w_perp_norm = float(np.linalg.norm(w_perp))
            q_perp_norm = float(np.linalg.norm(q_perp))
            if w_perp_norm < 1e-12 or q_perp_norm < 1e-12:  # pragma: no cover
                return None
            cos_chi = float(np.dot(w_perp, q_perp)) / (w_perp_norm * q_perp_norm)
            sin_chi = float(np.dot(n_chi, np.cross(w_perp, q_perp))) / (
                w_perp_norm * q_perp_norm
            )
            phi_d = phi_deg_f  # the scanned phi
            chi_d = math.degrees(math.atan2(sin_chi, cos_chi))

            trial = dict(angles)
            trial[chi_stage_name] = chi_d
            trial[phi_stage_name] = phi_d

            # Verify Bragg condition
            Q_computed = ctx.q_phi(trial)
            if float(np.linalg.norm(Q_computed - Q_phi)) > 1e-3:
                return None
            return trial

        # Scan the free angle (phi when v ∥ n_phi) at 2-degree
        # intervals, detect sign changes in the Ewald residual.
        candidates = []
        prev_scan: float | None = None
        prev_ew: float = 0.0

        for scan_int in range(-180, 180, 2):
            scan_f = float(scan_int)
            trial = _chi_to_trial(scan_f)
            if trial is None:
                prev_scan = None
                continue
            ew = _ewald_residual(trial)

            if abs(ew) < 1e-3:  # pragma: no cover
                candidates.append(trial)  # pragma: no cover
            elif prev_scan is not None and prev_ew * ew < 0:
                # Sign change: bisect to find root
                lo, hi = prev_scan, scan_f
                ew_lo = prev_ew
                for _ in range(40):
                    mid = (lo + hi) / 2.0
                    trial_mid = _chi_to_trial(mid)
                    if trial_mid is None:  # pragma: no cover
                        break  # pragma: no cover
                    ew_mid = _ewald_residual(trial_mid)
                    if abs(ew_mid) < 1e-9:
                        candidates.append(trial_mid)
                        break
                    if ew_lo * ew_mid < 0:
                        hi = mid
                    else:
                        lo = mid
                        ew_lo = ew_mid
                else:  # pragma: no cover
                    # Didn't converge to 1e-9 but may be close enough
                    trial_final = _chi_to_trial((lo + hi) / 2.0)
                    if trial_final is not None:
                        if abs(_ewald_residual(trial_final)) < 1e-5:
                            candidates.append(trial_final)

            prev_scan = scan_f
            prev_ew = ew

        return candidates

    def _solve_inner_for_outer(  # pragma: no cover — kappa only, tested by slow_benchmark
        outer_deg: float,
    ) -> list[dict[str, float]]:
        """
        Given a fixed outer sample angle, solve (chi, phi) from Q direction.

        Only used in the kappa seed-based path.

        Returns a list of candidate angle dicts (typically 0 or 2).
        """
        angles = dict(angles_base)
        angles[outer_stage_name] = outer_deg

        candidates = []

        if is_eulerian:
            # Fast analytic path (no ForwardContext creation)
            pairs = _solve_inner_eulerian_fast(outer_deg)
            for chi_d, phi_d in pairs:
                sol = dict(angles)
                sol[chi_stage_name] = chi_d
                sol[phi_stage_name] = phi_d
                candidates.append(sol)
        else:
            # Kappa/non-standard: 2D Gauss-Newton for (chi, phi)
            from .orientation import angles_to_phi_vector as _a2phi

            ctx = ForwardContext(geometry)
            ctx.prepare_caching(angles, {chi_stage_name, phi_stage_name})

            for chi_start, phi_start in [
                (0.0, 0.0),
                (90.0, 0.0),
                (-90.0, 0.0),
                (0.0, 90.0),
            ]:
                result = _solve_two_angles(
                    geometry=geometry,
                    fixed_angles=angles,
                    free_stage_1=chi_stage_name,
                    free_stage_2=phi_stage_name,
                    start_1=chi_start,
                    start_2=phi_start,
                    Q_phi_target=Q_phi,
                    a2phi_fn=_a2phi,
                    _ctx=ctx,
                )
                if result is not None:
                    sol = dict(angles)
                    sol[chi_stage_name] = result[0]
                    sol[phi_stage_name] = result[1]
                    # De-duplicate within this outer angle
                    dup = any(
                        abs(c[chi_stage_name] - sol[chi_stage_name]) < 1e-3
                        and abs(c[phi_stage_name] - sol[phi_stage_name]) < 1e-3
                        for c in candidates
                    )
                    if not dup:
                        candidates.append(sol)

        return candidates

    def _solve_inner_seeded(
        outer_deg: float, chi_hint: float, phi_hint: float
    ) -> dict[str, float] | None:
        """
        Solve (chi, phi) for a given outer angle, seeded from a previous
        nearby solution.  Used during 1D Newton refinement to track a branch
        cheaply (single seed instead of multi-seed scan).
        """
        angles = dict(angles_base)
        angles[outer_stage_name] = outer_deg

        if is_eulerian:
            pairs = _solve_inner_eulerian_fast(outer_deg)
            if not pairs:
                return None  # pragma: no cover
            # Pick the branch closest to the hint
            best = min(
                pairs,
                key=lambda p: abs(p[0] - chi_hint) + abs(p[1] - phi_hint),
            )
            sol = dict(angles)
            sol[chi_stage_name] = best[0]
            sol[phi_stage_name] = best[1]
            return sol
        else:  # pragma: no cover — kappa path tested by slow_benchmark
            from .orientation import angles_to_phi_vector as _a2phi

            ctx = ForwardContext(geometry)
            ctx.prepare_caching(angles, {chi_stage_name, phi_stage_name})
            result = _solve_two_angles(
                geometry=geometry,
                fixed_angles=angles,
                free_stage_1=chi_stage_name,
                free_stage_2=phi_stage_name,
                start_1=chi_hint,
                start_2=phi_hint,
                Q_phi_target=Q_phi,
                a2phi_fn=_a2phi,
                _ctx=ctx,
            )
            if result is None:
                return None
            sol = dict(angles)
            sol[chi_stage_name] = result[0]
            sol[phi_stage_name] = result[1]
            return sol

    found_solutions: list[dict[str, float]] = []
    seen_keys: list[np.ndarray] = []

    if is_eulerian:
        # ---------------------------------------------------------------
        # Eulerian fast path: dense 1D scan + sign-change detection
        # ---------------------------------------------------------------
        # The inlined analytic inner solve is very cheap (~0.1ms per point),
        # so we can afford a dense scan.  For each grid point we compute
        # the Ewald residual for both (chi, phi) branches and detect sign
        # changes between adjacent points.  A 1D Newton refinement then
        # drives the outer angle to the exact root.

        _SCAN_STEP = 0.5  # degrees
        _SCAN_LO = -180.0
        _SCAN_HI = 180.0
        n_scan = int((_SCAN_HI - _SCAN_LO) / _SCAN_STEP)

        # prev_branches stores (outer_deg, [(chi, phi, ewald_res), ...])
        prev_branches: list[tuple[float, float, float]] | None = None
        prev_outer: float = 0.0

        for i_scan in range(n_scan):
            outer_deg = _SCAN_LO + i_scan * _SCAN_STEP
            pairs = _solve_inner_eulerian_fast(outer_deg)
            if not pairs:
                prev_branches = None
                continue

            # Compute Ewald residual for each branch
            cur_branches = []
            for chi_d, phi_d in pairs:
                trial = dict(angles_base)
                trial[outer_stage_name] = outer_deg
                trial[chi_stage_name] = chi_d
                trial[phi_stage_name] = phi_d
                ew = _ewald_residual(trial)
                cur_branches.append((chi_d, phi_d, ew))

                # Direct hit at grid point
                if abs(ew) < 1e-3:  # pragma: no cover
                    refined = _refine_dd_outer_seeded(
                        outer_deg,
                        trial,
                        chi_d,
                        phi_d,
                        chi_stage_name,
                        phi_stage_name,
                        _solve_inner_seeded,
                        _ewald_residual,
                    )
                    if refined is not None:
                        _collect_dd_solution(
                            refined,
                            free_names,
                            found_solutions,
                            seen_keys,
                            mode,
                            geometry,
                        )

            # Detect sign changes between adjacent grid points
            if prev_branches is not None:
                for pc, pp, pr in prev_branches:
                    for cc, cp, cr in cur_branches:
                        # Same branch: (chi, phi) should be close
                        if abs(pc - cc) < 15.0 and abs(pp - cp) < 15.0:
                            if pr * cr < 0:
                                # Sign change — refine from the midpoint
                                mid = (prev_outer + outer_deg) / 2.0
                                mid_pairs = _solve_inner_eulerian_fast(mid)
                                if mid_pairs:  # pragma: no branch
                                    # Pick the branch closest to prev
                                    best = min(
                                        mid_pairs,
                                        key=lambda p: abs(p[0] - pc) + abs(p[1] - pp),
                                    )
                                    mid_sol = dict(angles_base)
                                    mid_sol[outer_stage_name] = mid
                                    mid_sol[chi_stage_name] = best[0]
                                    mid_sol[phi_stage_name] = best[1]
                                    refined = _refine_dd_outer_seeded(
                                        mid,
                                        mid_sol,
                                        best[0],
                                        best[1],
                                        chi_stage_name,
                                        phi_stage_name,
                                        _solve_inner_seeded,
                                        _ewald_residual,
                                    )
                                    if refined is not None:  # pragma: no branch
                                        _collect_dd_solution(
                                            refined,
                                            free_names,
                                            found_solutions,
                                            seen_keys,
                                            mode,
                                            geometry,
                                        )

            prev_branches = cur_branches
            prev_outer = outer_deg

        # Handle degenerate outer angles where the analytic solver fails
        # (q parallel to n_chi).  At these points chi is indeterminate
        # and solutions may exist at specific chi values that satisfy the
        # Ewald constraint.
        degen_outers = _find_degenerate_outers()
        for degen_outer in degen_outers:
            degen_cands = _solve_degenerate_outer(degen_outer)
            for cand in degen_cands:
                # Refine the outer angle to drive the Ewald residual to zero
                refined = _refine_dd_outer_seeded(
                    degen_outer,
                    cand,
                    cand[chi_stage_name],
                    cand[phi_stage_name],
                    chi_stage_name,
                    phi_stage_name,
                    _solve_inner_seeded,
                    _ewald_residual,
                )
                if refined is not None:  # pragma: no branch
                    _collect_dd_solution(
                        refined,
                        free_names,
                        found_solutions,
                        seen_keys,
                        mode,
                        geometry,
                    )

    else:  # pragma: no cover
        # ---------------------------------------------------------------
        # Kappa/non-standard: seed-based 1D Newton
        # ---------------------------------------------------------------
        # The 2D Newton inner solve is expensive, so we use targeted seeds
        # (bisecting solutions + coarse grid) and track branches cheaply.

        outer_seeds: list[tuple[float, float, float]] = []

        # Bisecting seeds
        from .mode import BisectConstraint as _BC
        from .mode import ConstraintSet as _CS

        bisect_sample = free_sample[0]
        bisect_det = free_det[-1]
        other_constraints = list(mode.constraints)
        synth = _CS(
            [_BC(bisect_sample, bisect_det)] + other_constraints,
            computed=mode.computed,
        )
        try:
            bisect_sols = _solve_bisecting(geometry, Q_phi, ttheta_deg, synth)
            for sol in bisect_sols:
                outer_seeds.append(
                    (
                        sol[outer_stage_name],
                        sol[chi_stage_name],
                        sol[phi_stage_name],
                    )
                )
        except Exception:  # pragma: no cover
            pass  # pragma: no cover

        # Coarse grid seeds (every 30 degrees — reduced for speed)
        half_tth = ttheta_deg / 2.0
        outer_seeds.append((half_tth, 0.0, 0.0))
        outer_seeds.append((-half_tth, 0.0, 0.0))
        for angle in range(-180, 180, 30):
            outer_seeds.append((float(angle), 0.0, 0.0))

        for seed_outer, _chi_hint, _phi_hint in outer_seeds:
            cands = _solve_inner_for_outer(seed_outer)
            for cand in cands:
                refined = _refine_dd_outer_seeded(
                    seed_outer,
                    cand,
                    cand[chi_stage_name],
                    cand[phi_stage_name],
                    chi_stage_name,
                    phi_stage_name,
                    _solve_inner_seeded,
                    _ewald_residual,
                )
                if refined is not None:
                    _collect_dd_solution(
                        refined,
                        free_names,
                        found_solutions,
                        seen_keys,
                        mode,
                        geometry,
                    )

    return found_solutions


def _refine_dd_outer_seeded(
    outer_deg: float,
    candidate: dict[str, float],
    chi_hint: float,
    phi_hint: float,
    chi_name: str,
    phi_name: str,
    solve_inner_seeded_fn,
    ewald_fn,
    max_iter: int = 20,
    h: float = 1e-4,
) -> dict[str, float] | None:
    """
    1D Newton refinement on the outer sample angle to satisfy the Ewald
    sphere constraint for the secondary reflection.

    Tracks a single (chi, phi) branch using the seeded inner solver, which
    uses the previous solution as the starting guess.  This avoids the
    multi-seed overhead of the full inner solver.

    Parameters
    ----------
    outer_deg : float
        Starting value for the outer sample angle.
    candidate : dict[str, float]
        Pre-computed candidate at ``outer_deg``.
    chi_hint, phi_hint : float
        Initial (chi, phi) hint for branch tracking.
    chi_name, phi_name : str
        Stage names for chi and phi.
    solve_inner_seeded_fn : callable
        ``(outer_deg, chi_hint, phi_hint) -> dict | None``
    ewald_fn : callable
        ``(angles) -> float``
    max_iter : int
    h : float
        Finite-difference step for the derivative.

    Returns
    -------
    dict[str, float] or None
        Refined solution dict, or None if refinement failed.
    """
    x = outer_deg
    sol = candidate
    c_hint = chi_hint
    p_hint = phi_hint

    for _ in range(max_iter):
        r = ewald_fn(sol)
        if abs(r) < 1e-9:
            return sol
        # Finite-difference derivative dr/d(outer)
        sol_p = solve_inner_seeded_fn(x + h, c_hint, p_hint)
        if sol_p is None:  # pragma: no cover
            return None  # pragma: no cover
        r_p = ewald_fn(sol_p)
        dr = (r_p - r) / h
        if abs(dr) < 1e-15:  # pragma: no cover
            return None  # flat — no root here  # pragma: no cover
        dx = -r / dr
        if abs(dx) > 5.0:  # pragma: no cover
            dx = 5.0 * (1.0 if dx > 0 else -1.0)  # pragma: no cover
        x += dx
        sol = solve_inner_seeded_fn(x, c_hint, p_hint)
        if sol is None:  # pragma: no cover
            return None  # pragma: no cover
        c_hint = sol[chi_name]
        p_hint = sol[phi_name]

    # Final check  # pragma: no cover
    r = ewald_fn(sol)  # pragma: no cover
    if abs(r) < 1e-5:  # pragma: no cover
        return sol  # pragma: no cover
    return None  # pragma: no cover


def _collect_dd_solution(
    sol: dict[str, float],
    free_names: list[str],
    found_solutions: list[dict[str, float]],
    seen_keys: list[np.ndarray],
    mode,
    geometry,
) -> None:
    """De-duplicate, apply cut-points, check limits, and append if valid."""
    key = np.array([sol[n] for n in free_names], dtype=float)
    if any(np.allclose(key, sk, atol=1e-3) for sk in seen_keys):  # pragma: no cover
        return  # pragma: no cover
    seen_keys.append(key)

    _apply_cut_points(sol, mode, geometry)
    if _check_limits(geometry, sol):  # pragma: no branch
        found_solutions.append(sol)


# ---------------------------------------------------------------------------
# Zone-mode solver (You 1999 §6, SPEC `setmode 5`)
# ---------------------------------------------------------------------------


def _is_zone_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when ``mode.extras`` contains both ``z0`` and ``z1`` keys.

    Zone modes confine the scattering vector Q to the plane spanned by two
    reciprocal-lattice vectors z0 and z1.  Like double-diffraction modes,
    zone modes are dispatched by their ``extras`` schema rather than by a
    dedicated constraint class.
    """
    extras = mode.extras
    return extras is not None and ("z0" in extras and "z1" in extras)


def _solve_zone(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """Forward solver for zone modes.

    The zone plane is fixed by two reciprocal-lattice vectors ``z0`` and
    ``z1`` supplied via ``mode.extras``.  For a single ``forward(h, k, l)``
    call the in-plane condition is a property of the requested (h, k, l)
    and the orientation matrix only — it does not constrain the motor
    angles directly.  The solver therefore:

    1. Validates ``z0``, ``z1`` are concrete and non-degenerate.
    2. Computes the zone-plane normal in the φ frame:
       ``n_zone_phi = (UB @ z0) × (UB @ z1)``.
    3. Verifies that the requested ``Q_phi`` lies in the plane.  If not,
       returns ``[]`` and emits a warning (matching the soft-failure
       pattern adopted for ``fixed_psi`` in #176).
    4. Records the in-plane residual in ``mode.extras['in_plane_residual']``.
    5. Builds a synthetic :class:`~mode.ConstraintSet` that adds a
       :class:`~mode.BisectConstraint` (or
       :class:`~mode.VirtualBisectConstraint` for kappa6c) so the system
       is fully constrained, then delegates to the existing bisecting
       solver.  Bisecting is the canonical SPEC-zone "br" position; the
       remaining azimuthal DOF can be exercised by a separate scan
       function (``cz``/``mz`` macros in SPEC), which is out of scope for
       a single ``forward()`` call.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
        Target scattering vector in the φ frame: ``UB @ (h, k, l)``.
    ttheta_deg : float
        Two-theta angle from Bragg's law.
    mode : ConstraintSet
        Active mode with ``extras['z0']``, ``extras['z1']`` set.

    Returns
    -------
    list of dict[str, float]
        Bisecting solutions filtered for plane membership.  Empty list when
        the requested (h, k, l) does not lie in the zone plane.

    Raises
    ------
    ValueError
        If ``z0`` or ``z1`` are still :data:`~mode.REQUIRED` placeholders,
        or if they are zero vectors, or if they are parallel (so the plane
        normal is undefined).
    """
    from .mode import REQUIRED
    from .mode import BisectConstraint
    from .mode import ConstraintSet
    from .mode import VirtualBisectConstraint

    extras = mode.extras
    z0_raw = extras.get("z0")
    z1_raw = extras.get("z1")

    # --- Validate z0/z1 -------------------------------------------------
    if z0_raw is REQUIRED or z1_raw is REQUIRED:
        raise ValueError(
            "zone mode requires z0 and z1 to be set in mode.extras "
            "before calling forward(). "
            "Set them with e.g. "
            "g.modes['zone_vertical'].extras['z0'] = (1, 0, 0); "
            "g.modes['zone_vertical'].extras['z1'] = (0, 1, 0)"
        )
    try:
        z0 = np.asarray(z0_raw, dtype=float).reshape(3)
        z1 = np.asarray(z1_raw, dtype=float).reshape(3)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"zone mode: z0 and z1 must each be a 3-element sequence of "
            f"Miller indices; got z0={z0_raw!r}, z1={z1_raw!r}."
        ) from exc

    if float(np.linalg.norm(z0)) < 1e-12 or float(np.linalg.norm(z1)) < 1e-12:
        raise ValueError(
            "zone mode: z0 and z1 must be non-zero reciprocal-lattice vectors; "
            f"got z0={z0_raw!r}, z1={z1_raw!r}."
        )

    # --- Compute the zone-plane normal in the φ frame -------------------
    UB = geometry.sample.UB
    z0_phi = UB @ z0
    z1_phi = UB @ z1
    n_zone_phi = np.cross(z0_phi, z1_phi)
    n_zone_norm = float(np.linalg.norm(n_zone_phi))
    if n_zone_norm < 1e-12:
        raise ValueError(
            "zone mode: z0 and z1 are parallel (or one is zero) so the "
            "zone-plane normal is undefined; "
            f"got z0={z0_raw!r}, z1={z1_raw!r}."
        )
    n_zone_hat = n_zone_phi / n_zone_norm

    # --- Plane-membership prefilter -------------------------------------
    in_plane_residual = float(abs(np.dot(Q_phi, n_zone_hat)))
    extras["in_plane_residual"] = in_plane_residual

    q_norm = float(np.linalg.norm(Q_phi))
    # Tolerance scaled by |Q| so the test is meaningful for both
    # short and long reciprocal vectors.
    tol = max(1e-8, 1e-6 * q_norm)
    if in_plane_residual > tol:
        logger.warning(
            "zone mode: requested Q (|Q|=%.6g) is not in the zone plane "
            "defined by z0=%s, z1=%s "
            "(|Q · n_zone| = %.3e > tolerance %.3e); returning no solutions.",
            q_norm,
            tuple(z0_raw) if hasattr(z0_raw, "__iter__") else z0_raw,
            tuple(z1_raw) if hasattr(z1_raw, "__iter__") else z1_raw,
            in_plane_residual,
            tol,
        )
        return []

    # --- Build the synthetic bisecting ConstraintSet --------------------
    # Vertical zone modes: bisect = (sample-stage rotating about transverse,
    #   detector-stage rotating about transverse) i.e. (eta or komega, delta).
    # Horizontal zone modes: bisect = (mu, nu).  ``mu`` is a real motor on
    #   both psic and kappa6c, so a literal BisectConstraint suffices.
    sample_names = {s.name for s in geometry.sample_stages}
    is_vertical = "nu" in {c.name for c in mode.constraints if hasattr(c, "name")}
    # ``mode.constraints`` for zone_vertical contains DetectorConstraint("nu",0)
    # and SampleConstraint("mu",0); for zone_horizontal it has
    # DetectorConstraint("delta",0) and SampleConstraint("eta"/"komega",0).
    # Distinguish on the detector-side fixed name:
    det_c = mode.detector_constraint
    is_vertical = det_c is not None and det_c.name == "nu"

    if is_vertical:
        # Vertical: bisect detector "delta" with the sample stage that
        # rotates about the same axis (transverse).  On psic that is
        # "eta"; on kappa6c true bisecting is on the *virtual* omega.
        if "komega" in sample_names and geometry.kappa_alpha_deg is not None:
            bisect_cs = VirtualBisectConstraint("omega", "delta")
        else:
            bisect_cs = BisectConstraint("eta", "delta")
        synth_constraints = [bisect_cs] + list(mode.constraints)
    else:
        # Horizontal: bisect detector "nu" with sample "mu" (a real motor
        # on both psic and kappa6c, rotating about the vertical axis).
        bisect_cs = BisectConstraint("mu", "nu")
        synth_constraints = [bisect_cs] + list(mode.constraints)

    synth_mode = ConstraintSet(
        synth_constraints,
        computed=mode.computed,
        cut_points=mode.cut_points,
    )

    # --- Delegate to the existing bisecting solver ----------------------
    # ``_solve_bisecting`` itself dispatches to the kappa-virtual solver
    # when the BisectConstraint is a VirtualBisectConstraint, so the
    # single call covers both psic (literal bisect) and kappa6c (virtual
    # bisect) zone modes uniformly.
    return _solve_bisecting(geometry, Q_phi, ttheta_deg, synth_mode)


# ---------------------------------------------------------------------------
# Surface diffraction solvers (ReferenceConstraint modes)
# ---------------------------------------------------------------------------


def _is_surface_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when mode has a surface ReferenceConstraint and surface_normal is set.

    The ``"psi"`` and ``"omega"`` ReferenceConstraints are NOT surface modes —
    they are handled by :func:`_is_psi_mode` / :func:`_solve_psi_mode` and
    :func:`_is_omega_mode` / :func:`_solve_omega_mode` respectively.

    psic-family modes that leave **two free sample stages** along with both
    detector stages free are routed to :func:`_solve_free_detectors`
    instead (the new general solver added for issue #264) — the legacy
    ``_solve_surface`` only rocks a single sample stage and assumes that
    fixing the active detector to ``ttheta`` plus the other constraints
    is enough to satisfy Bragg.
    """
    from .mode import ReferenceConstraint

    if geometry.surface_normal is None:
        return False
    has_surface_ref = any(
        isinstance(c, ReferenceConstraint) and c.name not in {"psi", "omega"}
        for c in mode.constraints
    )
    if not has_surface_ref:  # pragma: no cover
        return False
    # Defer to :func:`_is_free_detectors_mode` for psic-family modes
    # with 2 free sample stages and 2 free detector stages (issue
    # #264 — the B3 mode ``fixed_alpha_i_fixed_chi_fixed_phi``).  The
    # ``chi`` stage check restricts this exclusion to psic; existing
    # zaxis/s2d2/sixc surface modes (which also have free detectors)
    # stay on :func:`_solve_surface` where the legacy 1-D Newton works
    # well thanks to the constrained sample stack.
    has_chi = any(s.name == "chi" for s in geometry.sample_stages)
    if has_chi:
        fixed_sample_names = {c.name for c in mode.fixed_sample_constraints}
        n_free_sample = sum(
            1 for s in geometry.sample_stages if s.name not in fixed_sample_names
        )
        if (
            n_free_sample >= 2
            and mode.detector_constraint is None
            and len(geometry.detector_stages) >= 2
        ):
            return False
    return True


def _solve_surface(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Surface diffraction forward solver.

    Supports ReferenceConstraint modes where the constraint is:

    - ``"alpha_i"`` — incidence angle fixed at target value
    - ``"beta_out"`` — exit angle fixed at target value
    - ``"a_eq_b"`` — symmetric reflection: alpha_i = beta_out

    The solver builds a baseline angles dict (applying all fixed sample/detector
    constraints and setting the detector stage to ttheta_deg), then performs a
    1D Newton-Raphson root-find over the remaining free stage to satisfy the
    surface reference constraint.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
    """

    from .mode import ReferenceConstraint

    # Extract the reference constraint
    rc = next(c for c in mode.constraints if isinstance(c, ReferenceConstraint))
    target_name = rc.name  # "alpha_i", "beta_out", or "a_eq_b"
    target_value = rc.value  # float or True

    # Build baseline angles dict with all fixed constraints applied
    angles: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }

    # Apply fixed sample constraints
    for c in mode.fixed_sample_constraints:
        if c.name in geometry._stages:  # noqa: SLF001  # pragma: no branch
            angles[c.name] = float(c.value)

    # Apply detector constraint if present
    det_constraint = mode.detector_constraint
    if det_constraint is not None and not det_constraint.is_qaz:
        angles[det_constraint.name] = det_constraint.value

    # Select the *active* detector stage — the one that should carry ttheta.
    # When the mode's DetectorConstraint pins a specific detector stage
    # (e.g. ``delta=0`` in psic horizontal modes), the active stage is the
    # OTHER detector stage (e.g. ``nu``).  Picking the last stage
    # unconditionally would overwrite the pinned value (issue #279).
    pinned_det: str | None = None
    if det_constraint is not None and not det_constraint.is_qaz:
        pinned_det = det_constraint.name
    free_det_stages = [s for s in geometry.detector_stages if s.name != pinned_det]
    if not free_det_stages:
        # Every detector stage is pinned — ``_solve_surface`` cannot place
        # ttheta anywhere.  Caller should have routed to a different solver.
        return []
    # Prefer the last free detector stage to preserve historical behavior
    # for the common single-detector-stage geometries.
    det_stage = free_det_stages[-1]
    angles[det_stage.name] = ttheta_deg

    # Identify the one free sample stage (the one not fixed by any constraint)
    constrained_names = set(mode.constrained_stages(geometry))
    constrained_names.add(det_stage.name)
    if det_constraint is not None and not det_constraint.is_qaz:
        constrained_names.add(det_constraint.name)
    free_sample = [s for s in geometry.sample_stages if s.name not in constrained_names]

    if not free_sample:  # pragma: no cover
        # All stages fixed — evaluate the constraint and return if satisfied
        if _surface_residual(angles, geometry, target_name, target_value) < 1e-6:
            _apply_cut_points(angles, mode, geometry)
            if _check_limits(geometry, angles):
                return [angles]
        return []

    if len(free_sample) > 1:  # pragma: no branch
        # More than one free sample stage — use the first one (rocking stage)
        # and leave others at current values.  This handles cases where the
        # reference constraint only restricts one angle.
        pass

    free_stage = free_sample[0]

    # 1D Newton-Raphson over the free stage angle
    def residual(angle_val: float) -> float:
        trial = dict(angles)
        trial[free_stage.name] = angle_val
        return _surface_residual(trial, geometry, target_name, target_value)

    solutions = []
    # Seed from a grid of starting values
    lim_lo, lim_hi = free_stage.limits
    step = max(5.0, (lim_hi - lim_lo) / 40.0)
    seeds = list(np.arange(lim_lo + step / 2, lim_hi, step))

    seen_angles: list[float] = []
    for seed in seeds:
        x = float(seed)
        r0 = residual(x)
        for _ in range(60):
            h = 1e-4
            dr = (residual(x + h) - residual(x - h)) / (2 * h)
            if abs(dr) < 1e-12:
                break
            dx = -r0 / dr
            dx = float(np.clip(dx, -10.0, 10.0))
            x += dx
            r0 = residual(x)
            if abs(r0) < 1e-8:
                break
        if abs(r0) > 1e-5:
            continue
        # Check it's within limits
        if not (lim_lo <= x <= lim_hi):
            continue
        # Deduplicate
        if any(abs(x - xs) < 1e-3 for xs in seen_angles):
            continue
        seen_angles.append(x)
        sol = dict(angles)
        sol[free_stage.name] = x
        _apply_cut_points(sol, mode, geometry)
        if _check_limits(geometry, sol):  # pragma: no branch
            solutions.append(sol)

    return solutions


def _surface_residual(
    angles: dict[str, float],
    geometry: AdHocDiffractometer,
    target_name: str,
    target_value,
) -> float:
    """
    Compute the surface constraint residual for the given angles.

    Returns a float residual in degrees (zero = constraint satisfied).
    """
    from .reference import exit_angle as _beta_out
    from .reference import incidence_angle as _alpha_i

    if target_name == "alpha_i":
        ai = _alpha_i(geometry, angles=angles)
        return ai - float(target_value)
    if target_name == "beta_out":
        bo = _beta_out(geometry, angles=angles)
        return bo - float(target_value)
    # a_eq_b: alpha_i = beta_out
    ai = _alpha_i(geometry, angles=angles)
    bo = _beta_out(geometry, angles=angles)
    return ai - bo


# ---------------------------------------------------------------------------
# OMEGA pseudo-angle solver (SPEC psic omega-fixed family)
# ---------------------------------------------------------------------------


def _is_omega_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """
    Return True when the mode has an ``"omega"`` ReferenceConstraint and the
    geometry implements the OMEGA pseudo-angle (i.e. has a ``chi`` stage).

    OMEGA is the SPEC ``Q[6]`` pseudo-angle: angle between Q and the plane
    of the chi circle.  See
    :func:`~ad_hoc_diffractometer.reference.omega_pseudo`.
    """
    from .mode import ReferenceConstraint

    has_chi = any(s.name == "chi" for s in geometry.sample_stages)
    if not has_chi:
        return False
    return any(
        isinstance(c, ReferenceConstraint) and c.name == "omega"
        for c in mode.constraints
    )


def _solve_omega_mode(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Forward solver for ``ReferenceConstraint("omega", value)`` modes.

    OMEGA is the SPEC ``Q[6]`` pseudo-angle — the angle between the
    scattering vector Q and the plane of the chi circle.  When the
    target value is zero, OMEGA = 0 means Q lies in the chi-circle
    plane and the diffractometer is in the bisecting condition; non-zero
    OMEGA tilts Q out of that plane.

    Strategy
    --------
    The OMEGA pseudo-angle depends only on the **outer** sample stages
    (those below chi in the stack — ``mu`` and ``eta`` in psic) and the
    detector stages.  In each of the supported psic mode topologies
    exactly one outer stage is fixed and the other is free; we drive
    that free outer stage to satisfy ``OMEGA = target`` using a 1-D
    Newton root-find.  The inner sample stages ``(chi, phi)`` and the
    active detector angle are determined at each Newton step by the
    Bragg condition through the ordinary fixed-sample solver.

    For the special case ``target = 0`` the bisecting condition (``OMEGA
    = 0 ⇔ Q in the chi-circle plane``) is exact and the solver
    short-circuits to ``_solve_bisecting`` via a synthetic
    ``BisectConstraint``.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
    """
    from .mode import ConstraintSet
    from .mode import ReferenceConstraint
    from .mode import SampleConstraint
    from .reference import omega_pseudo

    # Extract the omega target value
    rc = next(
        c
        for c in mode.constraints
        if isinstance(c, ReferenceConstraint) and c.name == "omega"
    )
    omega_target = float(rc.value)

    # Identify the outer (pre-chi) sample stages and which of them is
    # currently free under the mode's SampleConstraints.
    fixed_names = {c.name for c in mode.fixed_sample_constraints}
    sample_stages = list(geometry.sample_stages)
    chi_index = next((i for i, s in enumerate(sample_stages) if s.name == "chi"), None)
    if chi_index is None:  # pragma: no cover
        return []  # _is_omega_mode would have rejected this geometry
    outer_stages = sample_stages[:chi_index]
    outer_free = [s for s in outer_stages if s.name not in fixed_names]

    # ---- Special case: target = 0 ⇒ exact bisecting --------------------
    if abs(omega_target) < 1e-9 and outer_free:
        synth = _synthetic_bisecting_for_omega(mode, geometry)
        if synth is not None:  # pragma: no branch
            return _solve_bisecting(geometry, Q_phi, ttheta_deg, synth)

    # ---- All outer stages fixed: OMEGA is determined; verify only ------
    if not outer_free:  # pragma: no cover
        # The mode supplies enough fixed sample constraints to pin OMEGA
        # entirely; reduce to a fixed-sample solve and accept solutions
        # that match the requested target.  Not used by the current
        # YAML modes (B1/B2 always leave one outer stage free).
        synth = _omega_to_fixed_sample_mode(mode)
        candidates = _solve_constraint_set_no_omega(geometry, Q_phi, ttheta_deg, synth)
        return [
            cand
            for cand in candidates
            if abs(omega_pseudo(geometry, angles=cand) - omega_target) < 1e-3
        ]

    # ---- General case: 1-D Newton on the free outer stage --------------
    rocking = outer_free[0]
    base_constraints = [
        c
        for c in mode.constraints
        if not (isinstance(c, ReferenceConstraint) and c.name == "omega")
    ]

    def trial_solutions(xi: float) -> list[dict[str, float]]:
        """Solve the Bragg condition with ``rocking`` fixed at ``xi``."""
        trial_constraints = list(base_constraints)
        trial_constraints.append(SampleConstraint(rocking.name, xi))
        trial_mode = ConstraintSet(trial_constraints, computed=mode.computed)
        return _solve_constraint_set_no_omega(geometry, Q_phi, ttheta_deg, trial_mode)

    # Seed sweep across the rocking stage's range
    lim_lo, lim_hi = rocking.limits
    span = lim_hi - lim_lo
    n_seeds = 24
    seed_step = span / n_seeds
    seeds = [lim_lo + (i + 0.5) * seed_step for i in range(n_seeds)]

    solutions: list[dict[str, float]] = []
    seen_x: list[float] = []

    for seed in seeds:
        x = float(seed)
        sols = trial_solutions(x)
        if not sols:  # pragma: no cover
            continue
        # Pick the solution that minimizes the chi-axis sign-flip
        sol = sols[0]
        r0 = omega_pseudo(geometry, angles=sol) - omega_target
        # Newton iteration with a numerical derivative
        for _ in range(50):  # pragma: no branch
            if abs(r0) < 1e-7:
                break
            h = 1e-3
            sols_p = trial_solutions(x + h)
            sols_m = trial_solutions(x - h)
            if not sols_p or not sols_m:  # pragma: no cover
                break
            r_p = omega_pseudo(geometry, angles=sols_p[0]) - omega_target
            r_m = omega_pseudo(geometry, angles=sols_m[0]) - omega_target
            dr = (r_p - r_m) / (2 * h)
            if abs(dr) < 1e-12:  # pragma: no cover
                break
            dx = -r0 / dr
            dx = float(np.clip(dx, -15.0, 15.0))
            x += dx
            sols = trial_solutions(x)
            if not sols:  # pragma: no cover
                break
            sol = sols[0]
            r0 = omega_pseudo(geometry, angles=sol) - omega_target
        if abs(r0) > 1e-4:  # pragma: no cover
            continue
        if not (lim_lo <= x <= lim_hi):  # pragma: no cover
            continue
        # Deduplicate by rocking-stage value
        if any(abs(x - xs) < 1e-3 for xs in seen_x):  # pragma: no cover
            continue
        seen_x.append(x)
        if not _check_limits(geometry, sol):  # pragma: no cover
            continue
        _apply_cut_points(sol, mode, geometry)
        solutions.append(sol)

    return solutions


def _solve_constraint_set_no_omega(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """Run the constraint-set dispatcher on ``mode`` after stripping any
    ``ReferenceConstraint("omega", ...)`` so the inner solve cannot recurse.

    The omega solver re-enters the dispatcher to handle the inner Bragg
    solve; without this guard the dispatcher would route the recursive
    call back into :func:`_solve_omega_mode` and loop forever.
    """
    from .mode import ConstraintSet
    from .mode import ReferenceConstraint

    if not any(
        isinstance(c, ReferenceConstraint) and c.name == "omega"
        for c in mode.constraints
    ):
        # Already free of omega — call the dispatcher directly
        return _solve_constraint_set(geometry, Q_phi, ttheta_deg, mode)

    stripped = ConstraintSet(  # pragma: no cover
        [
            c
            for c in mode.constraints
            if not (isinstance(c, ReferenceConstraint) and c.name == "omega")
        ],
        computed=mode.computed,
        extras=mode.extras,
        cut_points=mode.cut_points,
    )
    return _solve_constraint_set(
        geometry, Q_phi, ttheta_deg, stripped
    )  # pragma: no cover


def _omega_to_fixed_sample_mode(mode):  # pragma: no cover
    """Return ``mode`` with the omega ReferenceConstraint stripped.

    Used only by the "all outer stages fixed" branch of
    :func:`_solve_omega_mode`, which the current YAML modes never trigger.
    """
    from .mode import ConstraintSet
    from .mode import ReferenceConstraint

    return ConstraintSet(
        [
            c
            for c in mode.constraints
            if not (isinstance(c, ReferenceConstraint) and c.name == "omega")
        ],
        computed=mode.computed,
        extras=mode.extras,
        cut_points=mode.cut_points,
    )


def _synthetic_bisecting_for_omega(mode, geometry):
    """
    Build a synthetic ConstraintSet that replaces the omega ReferenceConstraint
    with a BisectConstraint on the geometry's chi-circle bisect pair.

    Returns ``None`` if no bisect pair can be identified for this mode.
    """
    from .mode import ConstraintSet
    from .mode import DetectorConstraint
    from .mode import ReferenceConstraint

    bisect = _bisect_pair_for(geometry, mode)
    if bisect is None:  # pragma: no cover
        return None
    others = [
        c
        for c in mode.constraints
        if not (isinstance(c, ReferenceConstraint) and c.name == "omega")
    ]
    # Drop any DetectorConstraint on the same detector stage that the
    # bisect would constrain — bisect drives that stage analytically.
    others = [
        c
        for c in others
        if not (
            isinstance(c, DetectorConstraint)
            and getattr(c, "name", None) == bisect.detector_stage
        )
    ]
    return ConstraintSet(
        [bisect] + others,
        computed=mode.computed,
        extras=mode.extras,
        cut_points=mode.cut_points,
    )


def _bisect_pair_for(geometry, mode):
    """
    Return a ``BisectConstraint`` on the geometry's chi-circle bisect pair.

    For psic the bisect pair is determined by which sample stage is fixed:

    - If ``mu`` is fixed at 0 (vertical scattering plane): ``eta ↔ delta``.
    - If ``eta`` is fixed at 0 (horizontal scattering plane): ``mu ↔ nu``.

    Returns ``None`` if neither condition is satisfied (e.g. an unsupported
    psic-like geometry).
    """
    from .mode import BisectConstraint

    fixed = {c.name: c.value for c in mode.fixed_sample_constraints}
    if fixed.get("mu") == 0.0:
        return BisectConstraint("eta", "delta")
    if fixed.get("eta") == 0.0:
        return BisectConstraint("mu", "nu")
    return None  # pragma: no cover


def _is_qaz_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when the mode contains a qaz DetectorConstraint."""
    from .mode import DetectorConstraint

    return any(isinstance(c, DetectorConstraint) and c.is_qaz for c in mode.constraints)


def _solve_qaz_mode(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Forward solver for modes with a ``DetectorConstraint("qaz", value)``.

    These are the ``lifting_detector_*`` family of modes (You 1999).  All
    sample stages are frozen by explicit :class:`~mode.SampleConstraint`
    entries; the solver finds the two detector stage angles ``(nu, delta)``
    that simultaneously satisfy the Bragg condition and the qaz constraint.

    Two equations, two unknowns ``(nu, delta)``:

    1. Bragg: ``cos(nu) * cos(delta) = cos(ttheta)``
    2. qaz:   ``atan2(tan(delta), sin(nu)) = qaz_target``

    The detector angles are solved analytically (no Newton iteration).
    For ``qaz = 90°`` (vertical scattering plane, the common case):
    ``nu = 0``, ``delta = ±ttheta``.  For the general case,
    ``sin²(nu) = (1 - cos²(ttheta)) / (1 + cos²(ttheta)·tan²(qaz))``
    and ``delta = atan(sin(nu)·tan(qaz))``.

    After finding the detector angles, any remaining free sample stages
    (those not frozen by the fixed sample constraints) are solved
    numerically via :func:`_solve_one_free_angle` or :func:`_solve_two_angles`.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
    """
    import math

    from .mode import DetectorConstraint

    # Identify the qaz constraint
    qaz_c = next(
        c for c in mode.constraints if isinstance(c, DetectorConstraint) and c.is_qaz
    )
    target_qaz_deg = float(qaz_c.value)

    # Identify detector stage names (outer = nu-like, inner = delta-like)
    det_stages = geometry.detector_stages
    nu_stage = det_stages[0]
    delta_stage = det_stages[-1]

    # Build baseline angles from current stage positions
    angles_base: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }

    # Apply fixed sample constraints
    for sc in mode.fixed_sample_constraints:
        if sc.name in geometry._stages:  # noqa: SLF001  # pragma: no branch
            angles_base[sc.name] = float(sc.value)

    # ---- Analytic detector-angle solver (issue #224) ----
    #
    # Two equations, two unknowns (nu, delta):
    #   1. Bragg:  cos(nu) * cos(delta) = cos(ttheta)
    #   2. qaz:    atan2(tan(delta), sin(nu)) = qaz_target
    #
    # From (2): tan(delta) = sin(nu) * tan(qaz).
    # Substituting into (1) via cos(delta) = 1/sqrt(1 + tan^2(delta)):
    #   cos(nu) / sqrt(1 + sin^2(nu)*tan^2(qaz)) = cos(ttheta)
    # which is directly solvable for nu, then delta follows from (2).
    #
    # Special case qaz = 90°:  nu = 0, delta = ±ttheta  (trivial).

    ttheta_rad = math.radians(ttheta_deg)
    cos_ttheta = math.cos(ttheta_rad)
    qaz_rad = math.radians(target_qaz_deg)

    # Collect all physically valid (nu, delta) pairs analytically.
    detector_pairs: list[tuple[float, float]] = []

    abs_qaz_mod = abs(target_qaz_deg % 360)
    _is_qaz_90 = abs(abs_qaz_mod - 90.0) < 1e-9 or abs(abs_qaz_mod - 270.0) < 1e-9

    if _is_qaz_90:
        # Fast path: qaz = 90° (or 270°) => nu = 0, delta has sign from qaz.
        # When nu = 0, qaz = atan2(tan(delta), sin(0)) = atan2(tan(delta), 0),
        # which equals +90° when delta > 0 and -90° when delta < 0.
        # Only the sign matching the target qaz is physically valid.
        qaz_sign = 1.0 if abs_qaz_mod < 180.0 else -1.0
        detector_pairs.append((0.0, qaz_sign * abs(ttheta_deg)))
    else:  # pragma: no cover — no demo geometry currently uses non-90° qaz
        # General analytic solution.
        # From Bragg + qaz combined:
        #   cos^2(nu) = cos^2(ttheta) * (1 + sin^2(nu)*tan^2(qaz))
        # Let s = sin(nu), c = cos(nu), T = tan^2(qaz), C = cos^2(ttheta):
        #   1 - s^2 = C*(1 + s^2*T)  =>  s^2*(1 + C*T) = 1 - C
        #   s^2 = (1 - C) / (1 + C*T)
        cos2_ttheta = cos_ttheta * cos_ttheta  # pragma: no cover
        tan_qaz = math.tan(qaz_rad)  # pragma: no cover
        tan2_qaz = tan_qaz * tan_qaz  # pragma: no cover
        denom = 1.0 + cos2_ttheta * tan2_qaz  # pragma: no cover
        if abs(denom) < 1e-15:  # pragma: no cover
            sin2_nu = 1.0  # pragma: no cover
        else:  # pragma: no cover
            sin2_nu = (1.0 - cos2_ttheta) / denom  # pragma: no cover

        if sin2_nu < -1e-12:  # pragma: no cover
            sin2_nu = -1.0  # pragma: no cover
        sin2_nu = min(1.0, sin2_nu)  # pragma: no cover

        if sin2_nu < 0.0:  # pragma: no cover
            pass  # pragma: no cover
        else:  # pragma: no cover
            sin_nu = math.sqrt(sin2_nu)  # pragma: no cover
            nu_candidates = []  # pragma: no cover
            if sin_nu < 1e-12:  # pragma: no cover
                nu_candidates.append(0.0)  # pragma: no cover
            else:  # pragma: no cover
                nu_candidates.append(
                    math.degrees(math.asin(sin_nu))
                )  # pragma: no cover
                nu_candidates.append(
                    math.degrees(math.asin(-sin_nu))
                )  # pragma: no cover

            for nu_deg in nu_candidates:  # pragma: no cover
                nu_r = math.radians(nu_deg)  # pragma: no cover
                cos_nu = math.cos(nu_r)  # pragma: no cover
                if abs(cos_nu) < 1e-15:  # pragma: no cover
                    continue  # pragma: no cover
                cos_delta = cos_ttheta / cos_nu  # pragma: no cover
                cos_delta = max(-1.0, min(1.0, cos_delta))  # pragma: no cover
                tan_delta = math.sin(nu_r) * tan_qaz  # pragma: no cover
                delta_deg = math.degrees(math.atan2(tan_delta, 1.0))  # pragma: no cover
                actual_cos_delta = math.cos(math.radians(delta_deg))  # pragma: no cover
                if actual_cos_delta * cos_delta < -1e-6:  # pragma: no cover
                    delta_abs = math.degrees(math.acos(cos_delta))  # pragma: no cover
                    delta_deg = math.copysign(delta_abs, tan_delta)  # pragma: no cover
                detector_pairs.append((nu_deg, delta_deg))  # pragma: no cover
                if abs(delta_deg) > 1e-10:  # pragma: no cover
                    neg_cos_delta = math.cos(
                        math.radians(-delta_deg)
                    )  # pragma: no cover
                    if abs(neg_cos_delta - cos_delta) < 1e-6:  # pragma: no cover
                        detector_pairs.append((nu_deg, -delta_deg))  # pragma: no cover

    found_solutions = []
    seen_keys: list[tuple[float, float]] = []

    for nu_sol, delta_sol in detector_pairs:
        # De-duplicate detector angle pairs
        duplicate = any(
            abs(nu_sol - ks[0]) < 1e-4 and abs(delta_sol - ks[1]) < 1e-4
            for ks in seen_keys
        )
        if duplicate:  # pragma: no cover — analytic solver produces no duplicates
            continue  # pragma: no cover
        seen_keys.append((nu_sol, delta_sol))

        # Build the full angle dict for this detector solution
        angles = dict(angles_base)
        angles[nu_stage.name] = nu_sol
        angles[delta_stage.name] = delta_sol

        # Check detector stage limits
        if not nu_stage.in_limits(nu_sol) or not delta_stage.in_limits(delta_sol):
            continue

        # Identify free sample stages (not frozen by any constraint)
        constrained_names = set(mode.constrained_stages(geometry))
        constrained_names.add(nu_stage.name)
        constrained_names.add(delta_stage.name)
        free_sample = [
            s for s in geometry.sample_stages if s.name not in constrained_names
        ]

        if len(free_sample) == 0:  # pragma: no cover
            # All sample stages fixed — validate Q direction
            from .orientation import angles_to_phi_vector as _a2phi  # pragma: no cover

            Q_computed = _a2phi(geometry, **angles)  # pragma: no cover
            if not np.allclose(Q_computed, Q_phi, atol=1e-3):  # pragma: no cover
                continue  # pragma: no cover
            _apply_cut_points(angles, mode, geometry)  # pragma: no cover
            if _check_limits(geometry, angles):  # pragma: no cover
                found_solutions.append(angles)  # pragma: no cover

        elif len(free_sample) == 1:  # pragma: no cover
            sub = _solve_one_free_angle(  # pragma: no cover
                geometry, angles, free_sample[0], Q_phi, mode
            )
            found_solutions.extend(sub)  # pragma: no cover

        else:
            # Two+ free sample stages: solve for chi and phi
            from .orientation import angles_to_phi_vector as _a2phi

            chi_stage = free_sample[-2]
            phi_stage = free_sample[-1]

            # Create ForwardContext for this detector angle pair
            qaz_ctx = ForwardContext(geometry)
            qaz_ctx.prepare_caching(angles, {chi_stage.name, phi_stage.name})

            Q_norm = float(np.linalg.norm(Q_phi))
            q_hat = Q_phi / Q_norm
            chi_ax = chi_stage._axis_hat  # noqa: SLF001
            q_along_chi = float(np.dot(q_hat, chi_ax))
            q_along_chi = max(-1.0, min(1.0, q_along_chi))
            chi_abs_deg = math.degrees(math.asin(abs(q_along_chi)))
            q_in_plane = q_hat - q_along_chi * chi_ax
            if np.linalg.norm(q_in_plane) > 1e-8:
                phi_seed_deg = math.degrees(
                    math.atan2(float(q_in_plane[1]), float(q_in_plane[0]))
                )
            else:  # pragma: no cover
                phi_seed_deg = 0.0  # pragma: no cover
            seed_pairs_inner = [
                (chi_s, phi_s)
                for chi_s in [chi_abs_deg, -chi_abs_deg, 180.0 - chi_abs_deg]
                for phi_s in [phi_seed_deg, phi_seed_deg + 180.0]
            ] + [
                (cs, ps)
                for cs in (0.0, 90.0, -90.0)
                for ps in (0.0, 90.0, 180.0, 270.0)
            ]
            for chi_start_i, phi_start_i in seed_pairs_inner:
                sol_angles = _solve_two_angles(
                    geometry=geometry,
                    fixed_angles=angles,
                    free_stage_1=chi_stage.name,
                    free_stage_2=phi_stage.name,
                    start_1=chi_start_i,
                    start_2=phi_start_i,
                    Q_phi_target=Q_phi,
                    a2phi_fn=_a2phi,
                    _ctx=qaz_ctx,
                )
                if sol_angles is None:
                    continue
                sol = dict(angles)
                sol[chi_stage.name] = sol_angles[0]
                sol[phi_stage.name] = sol_angles[1]
                _apply_cut_points(sol, mode, geometry)
                dup = any(
                    all(abs(ex.get(k, 0) - sol.get(k, 0)) < 1e-4 for k in sol)
                    for ex in found_solutions
                )
                if not dup and _check_limits(geometry, sol):
                    found_solutions.append(sol)

    return found_solutions


# ---------------------------------------------------------------------------
# Free-detectors solver (issue #264 — both detector stages free)
# ---------------------------------------------------------------------------


def _is_free_detectors_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when both detector stages are free under ``mode``.

    Used by :func:`_solve_free_detectors` to dispatch psic-family modes
    that fix multiple sample stages and let the detector arm orient
    itself entirely from the Bragg condition (and any reference
    constraint such as ``alpha_i``).  Examples (issue #264):

    - ``lifting_detector_eta`` (3 sample fixed, 1 sample + 2 detector free)
    - revised ``lifting_detector_phi`` / ``lifting_detector_mu``
      (after step C of #264 — same shape, different fixed sample stage)
    - ``fixed_alpha_i_fixed_chi_fixed_phi`` (2 sample fixed + alpha_i,
      2 sample + 2 detector free)

    The predicate is intentionally conservative: it requires the
    psic-family geometry signature (a sample stage named ``"chi"``),
    two detector stages with neither pinned by ``DetectorConstraint``,
    no qaz/bisect constraint, and a free-sample count compatible with
    the equation count (Bragg gives 3 equations; each
    :class:`~mode.ReferenceConstraint` other than ``"omega"`` adds 1).
    Existing surface geometries (zaxis, s2d2, sixc) lack ``chi`` and
    continue to use :func:`_solve_surface`.
    """
    from .mode import ReferenceConstraint

    detectors = geometry.detector_stages
    if len(detectors) < 2:
        return False

    # Restrict to psic-family geometries to avoid disturbing the
    # zaxis/s2d2/sixc surface mode dispatch path.
    if not any(s.name == "chi" for s in geometry.sample_stages):
        return False

    det_constraint = mode.detector_constraint
    if det_constraint is not None:
        # qaz is handled by _is_qaz_mode; any other detector constraint
        # pins one detector stage and is handled by _solve_fixed_sample.
        return False

    if mode.has_bisect:  # pragma: no cover
        return False

    fixed_sample_names = {c.name for c in mode.fixed_sample_constraints}
    free_sample = [
        s for s in geometry.sample_stages if s.name not in fixed_sample_names
    ]
    n_free_sample = len(free_sample)

    # Count reference equations (excluding "omega", which has its own
    # solver at _solve_omega_mode).
    n_ref_eqs = sum(
        1
        for c in mode.constraints
        if isinstance(c, ReferenceConstraint) and c.name != "omega"
    )

    # Equations available: 3 (Bragg) + n_ref_eqs.
    # Unknowns: n_free_sample + 2 (both detectors free).
    n_unknowns = n_free_sample + 2
    n_equations = 3 + n_ref_eqs

    return n_unknowns == n_equations and n_free_sample >= 1


def _solve_free_detectors(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Numerical Newton solver for psic-family modes with both detector
    stages free.

    Adds support for the issue #264 modes that drop the qaz constraint
    and let nu and delta float jointly to satisfy the Bragg condition:

    - ``lifting_detector_eta`` (B4)
    - revised ``lifting_detector_phi`` / ``lifting_detector_mu`` (C3, C4)
    - ``fixed_alpha_i_fixed_chi_fixed_phi`` (B3, with one alpha_i row)

    Variables: every free sample stage plus both detector stages.
    Equations: 3 from Bragg (``Q_phi(angles) == Q_phi_target``) plus 1 per
    non-``omega`` :class:`~mode.ReferenceConstraint` in the mode.

    The solver runs a finite-difference Levenberg-Marquardt-flavoured
    Newton iteration from a small grid of seed points and de-duplicates
    by the rounded free-stage values.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
        Magnitude of 2θ implied by ``Q_phi`` (used as a default seed for
        the active detector stage; not enforced as a constraint).
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
    """
    from .mode import ReferenceConstraint

    sample_stages = list(geometry.sample_stages)
    detector_stages = list(geometry.detector_stages)

    # Baseline angles (apply fixed sample constraints; detectors free)
    angles_base: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }
    fixed_sample_names = set()
    for c in mode.fixed_sample_constraints:  # pragma: no branch
        if c.name in geometry._stages:  # noqa: SLF001  # pragma: no branch
            angles_base[c.name] = float(c.value)
            fixed_sample_names.add(c.name)

    free_sample = [s for s in sample_stages if s.name not in fixed_sample_names]
    free_det = list(detector_stages)
    free_stages = free_sample + free_det
    free_names = [s.name for s in free_stages]
    n_free = len(free_names)

    # Reference constraints contributing extra equations
    ref_constraints = [
        c
        for c in mode.constraints
        if isinstance(c, ReferenceConstraint) and c.name != "omega"
    ]
    n_ref = len(ref_constraints)
    n_eqs = 3 + n_ref

    if n_eqs != n_free:  # pragma: no cover
        # _is_free_detectors_mode already enforces this, but be defensive
        return []

    def residual(x: np.ndarray) -> np.ndarray:
        """Combined Bragg + reference residual."""
        trial = dict(angles_base)
        for name, val in zip(free_names, x, strict=False):
            trial[name] = float(val)
        from .orientation import _compute_q_phi

        two_pi_over_lambda = 2.0 * np.pi / geometry.wavelength
        y_hat = np.asarray(geometry.basis["longitudinal"], dtype=float)
        y_hat = y_hat / np.linalg.norm(y_hat)
        Q_trial = _compute_q_phi(
            sample_stages, detector_stages, trial, two_pi_over_lambda, y_hat
        )
        r = np.zeros(n_eqs)
        r[:3] = Q_trial - Q_phi
        for i, rc in enumerate(ref_constraints):
            r[3 + i] = _surface_residual(trial, geometry, rc.name, rc.value)
        return r

    def jacobian_fd(x: np.ndarray, r0: np.ndarray, h: float = 1e-4) -> np.ndarray:
        J = np.zeros((n_eqs, n_free))
        for j in range(n_free):
            xp = x.copy()
            xp[j] += h
            J[:, j] = (residual(xp) - r0) / h
        return J

    # Seed grid: combine bisecting-style seeds for sample stages with
    # ttheta-anchored seeds for detector stages.
    sample_seeds: list[list[float]] = []
    for s in free_sample:
        lo, hi = s.limits
        center = max(lo, min(hi, 0.0))
        candidates = [center, ttheta_deg / 2.0, -ttheta_deg / 2.0, 30.0, -30.0]
        sample_seeds.append(
            sorted({float(c) for c in candidates if lo <= c <= hi}) or [center]
        )
    det_seed_pairs: list[tuple[float, float]] = []
    nu_lo, nu_hi = free_det[0].limits
    delta_lo, delta_hi = free_det[1].limits
    # Common starting points for (nu, delta): in-plane and lifted variants
    candidate_pairs = [
        (0.0, ttheta_deg),
        (0.0, -ttheta_deg),
        (ttheta_deg, 0.0),
        (-ttheta_deg, 0.0),
        (ttheta_deg / 2.0, ttheta_deg / 2.0),
        (-ttheta_deg / 2.0, ttheta_deg / 2.0),
    ]
    for nu0, d0 in candidate_pairs:  # pragma: no branch
        if nu_lo <= nu0 <= nu_hi and delta_lo <= d0 <= delta_hi:  # pragma: no branch
            det_seed_pairs.append((float(nu0), float(d0)))
    if not det_seed_pairs:  # pragma: no cover
        det_seed_pairs.append((0.0, float(np.clip(ttheta_deg, delta_lo, delta_hi))))

    # Build cartesian product of seed combinations (capped to keep the
    # work bounded for high-dimensional cases).
    from itertools import product

    sample_combos = list(product(*sample_seeds)) if sample_seeds else [()]
    seed_iter = list(product(sample_combos, det_seed_pairs))

    solutions: list[dict[str, float]] = []
    seen: list[tuple[float, ...]] = []

    for sample_combo, (nu0, d0) in seed_iter:
        x0 = list(sample_combo) + [nu0, d0]
        x = np.array(x0, dtype=float)

        # Levenberg-Marquardt with a damping parameter for stability.
        lam = 1e-3
        for _ in range(80):
            r = residual(x)
            err = float(np.linalg.norm(r))
            if err < 1e-9:
                break
            J = jacobian_fd(x, r)
            JtJ = J.T @ J
            try:
                dx = np.linalg.solve(
                    JtJ + lam * np.eye(n_free),
                    -J.T @ r,
                )
            except np.linalg.LinAlgError:  # pragma: no cover
                break
            step = float(np.linalg.norm(dx))
            if step > 30.0:
                dx = dx * (30.0 / step)
            x_new = x + dx
            r_new = residual(x_new)
            err_new = float(np.linalg.norm(r_new))
            if err_new < err:
                lam = max(lam * 0.5, 1e-9)
                x = x_new
            else:
                lam = min(lam * 2.0, 1e3)

        r_final = residual(x)
        if float(np.linalg.norm(r_final)) > 1e-5:
            continue

        # Build the solution dict
        sol = dict(angles_base)
        for name, val in zip(free_names, x, strict=False):
            # Normalize to (-180, 180]
            sol[name] = ((float(val) + 180.0) % 360.0) - 180.0

        # Limits check
        if not _check_limits(geometry, sol):  # pragma: no cover
            continue

        # Apply cut points first so equivalent wrap representatives
        # (e.g. -180 vs +180) collapse before dedup.
        _apply_cut_points(sol, mode, geometry)

        # Deduplicate using a modular key so that representative angles
        # one full turn apart (e.g. -180 vs +180) are recognized as the
        # same physical setting independent of cut-point configuration.
        # Round to 1 decimal place (0.1°) — finer than the typical
        # numerical noise from the Newton iteration but well below any
        # physically meaningful resolution.  The two-step round-then-mod
        # collapses near-zero negatives (e.g. -1e-15) that Python's ``%``
        # would otherwise wrap to ~360.
        key = tuple(round(round(sol[n], 1) % 360.0, 1) for n in free_names)
        if key in seen:
            continue
        seen.append(key)

        solutions.append(sol)

    return solutions


def _solve_fixed_sample(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Fixed-sample-angle solver (no bisect condition).

    Freezes all named sample stages at their constraint values, sets the
    detector stage to ``ttheta_deg``, and solves for any remaining free
    sample stages numerically.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : ConstraintSet

    Returns
    -------
    list of dict[str, float]
    """
    sample_stages = geometry.sample_stages

    # Build the baseline angle dict.
    angles: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }

    # Apply fixed sample constraints.
    for sc in mode.fixed_sample_constraints:
        if sc.name in geometry._stages:  # noqa: SLF001
            angles[sc.name] = float(sc.value)

    # Apply detector constraint if present.
    det_constraint = mode.detector_constraint
    if det_constraint is not None and not det_constraint.is_qaz:
        angles[det_constraint.name] = det_constraint.value
        active_det_name = None
        for ds in geometry.detector_stages:
            if ds.name != det_constraint.name:
                active_det_name = ds.name
                break
        if active_det_name is None:
            active_det_name = geometry.detector_stages[-1].name
    else:
        active_det_name = geometry.detector_stages[-1].name

    angles[active_det_name] = ttheta_deg

    # Free sample stages: those not constrained.
    constrained_names = set(mode.constrained_stages(geometry))
    constrained_names.add(active_det_name)
    if det_constraint is not None and not det_constraint.is_qaz:
        constrained_names.add(det_constraint.name)
    free_sample = [s for s in sample_stages if s.name not in constrained_names]

    if len(free_sample) == 0:
        _apply_cut_points(angles, mode, geometry)
        if _check_limits(geometry, angles):
            return [angles]
        return []

    if len(free_sample) == 1:
        return _solve_one_free_angle(geometry, angles, free_sample[0], Q_phi, mode)

    chi_stage = free_sample[-2]
    phi_stage = free_sample[-1]

    from .orientation import angles_to_phi_vector as _a2phi

    # Create ForwardContext and cache fixed-stage rotation matrices
    ctx = ForwardContext(geometry)
    free_names = {chi_stage.name, phi_stage.name}
    ctx.prepare_caching(angles, free_names)

    grid_seeds = [
        (chi_s, phi_s)
        for chi_s in (0.0, 45.0, 90.0, -90.0, 135.0, 180.0)
        for phi_s in (0.0, 90.0, 180.0, 270.0)
    ]
    solutions = []
    _MAX_SOLUTIONS = 4
    for chi_start, phi_start in grid_seeds:
        sol_angles = _solve_two_angles(
            geometry=geometry,
            fixed_angles=angles,
            free_stage_1=chi_stage.name,
            free_stage_2=phi_stage.name,
            start_1=chi_start,
            start_2=phi_start,
            Q_phi_target=Q_phi,
            a2phi_fn=_a2phi,
            _ctx=ctx,
        )
        if sol_angles is None:  # pragma: no branch
            continue  # pragma: no cover
        sol = dict(angles)
        sol[chi_stage.name] = sol_angles[0]
        sol[phi_stage.name] = sol_angles[1]
        _apply_cut_points(sol, mode, geometry)
        duplicate = False
        for existing in solutions:
            if all(abs(existing.get(k, 0) - sol.get(k, 0)) < 1e-4 for k in sol):
                duplicate = True
                break
        if not duplicate and _check_limits(geometry, sol):
            solutions.append(sol)
            if len(solutions) >= _MAX_SOLUTIONS:
                break
    return solutions


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _apply_cut_points(
    angles: dict[str, float],
    mode,
    geometry: AdHocDiffractometer,
) -> None:
    """
    Apply cut-points (in-place) to all angles in the solution dict.

    Mode-level cut-points take priority over geometry-level cut-points.
    """
    for stage_name in list(angles.keys()):
        if stage_name in mode.cut_points:
            angles[stage_name] = mode.apply_cut_point(stage_name, angles[stage_name])
        elif stage_name in geometry.cut_points:
            cut = geometry.cut_points[stage_name]
            remainder = (angles[stage_name] - cut) % 360.0
            angles[stage_name] = cut + remainder


def _validate_solutions(
    solutions: list[dict[str, float]],
    mode,
    geometry: AdHocDiffractometer,
) -> None:
    """
    Validate that every returned solution actually satisfies all constraints
    in the mode, raising ``ValueError`` if any constraint is violated beyond
    the display precision tolerance.

    This catches solver bugs and virtual-angle mismatches (e.g. a kappa
    geometry where the declared ``SampleConstraint`` value cannot be achieved
    by the approximate solver).

    Uses :func:`~.display.precision_atol` as the tolerance — half a unit in
    the last displayed decimal place, consistent with the package's display
    precision setting.

    Parameters
    ----------
    solutions : list of dict[str, float]
        Motor-angle dicts returned by the solver.
    mode : ConstraintSet
        The active constraint set.
    geometry : AdHocDiffractometer
        The diffractometer (needed for :meth:`~.mode.ConstraintSet.is_implemented`
        checks and stage-name resolution).

    Raises
    ------
    ConstraintViolation
        If any solution violates a constraint beyond tolerance.
    """
    from .display import precision_atol
    from .mode import BisectConstraint
    from .mode import ConstraintViolation
    from .mode import SampleConstraint

    if not solutions:
        return

    atol = precision_atol()

    for i, sol in enumerate(solutions):
        for c in mode.constraints:
            # Only validate constraints that have a numeric residual we can check.
            # BisectConstraint and SampleConstraint are the checkable ones;
            # DetectorConstraint is applied by the solver directly (ttheta_deg);
            # ReferenceConstraint is not yet implemented.
            if isinstance(c, BisectConstraint | SampleConstraint):
                try:
                    residual = c.evaluate(sol, geometry)
                except KeyError:
                    # Stage not in solution dict (virtual angle not yet computed).
                    continue
                if abs(residual) > atol:
                    if isinstance(c, BisectConstraint):
                        constraint_repr = (
                            f"BisectConstraint({c.sample_stage!r}, {c.detector_stage!r}): "
                            f"expected {c.sample_stage}="
                            f"{sol.get(c.detector_stage, '?')!s}/2, "
                            f"got {c.sample_stage}={sol.get(c.sample_stage, '?')!s}"
                        )
                    else:
                        constraint_repr = (
                            f"SampleConstraint({c.name!r}, {c.value!r}): "
                            f"got {c.name}={sol.get(c.name, '?')!s}"
                        )
                    raise ConstraintViolation(
                        solution_index=i,
                        constraint_repr=constraint_repr,
                        residual=residual,
                        tolerance=atol,
                    )


def _check_limits(
    geometry: AdHocDiffractometer,
    angles: dict[str, float],
) -> bool:
    """
    Return True if all angles are within the stage limits of their stage.

    Stages not present in the angles dict are not checked.
    """
    for name, angle in angles.items():
        try:
            stage = geometry.stage(name)
        except KeyError:
            continue
        if not stage.in_limits(angle):
            return False
    return True
