# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
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

    Created once per ``compute_forward()`` call, this bundles all constant
    quantities needed by the Newton-Raphson residual evaluations so they
    are not recomputed on every call to ``angles_to_phi_vector``.

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

        # Detector: if no detector stage is free, cache D entirely
        det_free = any(s.name in free_stage_names for s in self.detector_stages)
        if not det_free:
            D = np.eye(3)
            for s in self.detector_stages:
                angle = fixed_angles.get(s.name, s.angle)
                D = _rotation_matrix_normalized(s._axis_hat, angle) @ D  # noqa: SLF001
            self._cached_D = D
        else:
            self._cached_D = None

        # Sample: find the first free stage index and cache the prefix product
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
                    Z_prefix = (
                        _rotation_matrix_normalized(s._axis_hat, angle) @ Z_prefix
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

        # Build suffix products:  suffix[k] = R_{k-1} @ ... @ R_0
        # where indices are relative to first_free.
        # suffix[0] = cached_Z_prefix (product of stages 0..first_free-1).
        # suffix[k] = R_{first_free+k-1} @ suffix[k-1]
        suffix = [np.empty((3, 3))] * (n_tail + 1)
        suffix[0] = (
            self._cached_Z_prefix if self._cached_Z_prefix is not None else np.eye(3)
        )
        for k in range(n_tail):
            suffix[k + 1] = R_list[k] @ suffix[k]

        # Build prefix products:  prefix[k] = R_{N-1} @ ... @ R_{k+1}
        # where indices are relative to first_free.
        # prefix[n_tail] = I  (nothing after the last stage)
        # prefix[k] = prefix[k+1] @ R_{k+1}
        #
        # Actually we need prefix in absolute terms:
        # prefix_abs[j] = R_{N-1} @ ... @ R_{j+1}
        # For j = first_free + k:
        #   prefix[k] = R_{n_tail-1} @ ... @ R_{k+1}  (relative indices)
        prefix = [np.empty((3, 3))] * (n_tail + 1)
        prefix[n_tail] = np.eye(3)
        for k in range(n_tail - 1, -1, -1):
            prefix[k] = prefix[k + 1] @ R_list[k + 1] if k + 1 < n_tail else np.eye(3)
        # Recalculate more carefully:
        # prefix[k] should be the product of R_list[k+1] ... R_list[n_tail-1]
        # applied left to right in stacking order (outermost first).
        # Actually the Z chain is: Z = R_{N-1} @ R_{N-2} @ ... @ R_0
        # So for relative index k, R_list[k] = R_{first_free + k}.
        # Z_tail = R_list[n_tail-1] @ R_list[n_tail-2] @ ... @ R_list[0]
        #        = suffix[n_tail] (already computed above, without prefix)
        #
        # dZ/dθ_j for j = first_free + k:
        #   = (R_list[n_tail-1] @ ... @ R_list[k+1]) @ dR_list[k] @ (R_list[k-1] @ ... @ R_list[0] @ Z_prefix)
        #   = prefix[k] @ dR_list[k] @ suffix[k]
        #
        # Rebuild prefix properly:
        prefix[n_tail] = np.eye(3)
        for k in range(n_tail - 1, -1, -1):
            prefix[k] = prefix[k + 1] @ R_list[k + 1] if (k + 1) < n_tail else np.eye(3)

        # deg2rad factor: dR is w.r.t. radians, but angles are in degrees.
        deg2rad = np.pi / 180.0

        J = np.zeros((3, n_free))
        for col, name in enumerate(free_names):
            # Find the absolute stage index for this free name
            abs_idx = next(i for i, s in enumerate(stages) if s.name == name)
            k = abs_idx - first_free  # relative index into R_list
            dR_k = dR_map[abs_idx]
            # dZ/dθ = prefix[k] @ dR_k @ suffix[k]
            dZ = prefix[k] @ dR_k @ suffix[k]
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
    return solutions


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

    if mode.has_bisect:
        return _solve_bisecting(geometry, Q_phi, ttheta_deg, mode)

    # Virtual kappa angle mode (omega, chi, phi on a kappa geometry).
    if _is_kappa_virtual_mode(geometry, mode):
        return _solve_kappa_virtual(geometry, Q_phi, ttheta_deg, mode)

    # Surface diffraction mode (ReferenceConstraint with surface_normal).
    if _is_surface_mode(geometry, mode):
        return _solve_surface(geometry, Q_phi, ttheta_deg, mode)

    # qaz detector constraint mode (lifting_detector_* family).
    if _is_qaz_mode(geometry, mode):
        return _solve_qaz_mode(geometry, Q_phi, ttheta_deg, mode)

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

    The sample rotation matrix is:

    .. math::

        Z = R_\phi \cdot R_\chi \cdot Z_{\text{prefix}}

    where :math:`Z_{\text{prefix}}` is the product of all fixed outer
    sample stages.  The scattering vector in the phi frame is:

    .. math::

        Q_\phi = Z^T Q_{\text{lab}}
              = Z_{\text{prefix}}^T \cdot R_\chi^T \cdot R_\phi^T \cdot Q_{\text{lab}}

    Define the known vectors:

    - :math:`q = Z_{\text{prefix}} \cdot Q_{\phi,\text{target}}`
    - :math:`v = Q_{\text{lab}}`

    The equation to solve is:

    .. math::

        R(\hat n_\chi, -\chi) \cdot R(\hat n_\phi, -\phi) \cdot v = q

    Since :math:`\hat n_\chi \perp \hat n_\phi`, define
    :math:`\hat n_3 = \hat n_\chi \times \hat n_\phi`.  Projecting onto
    the :math:`(\hat n_\phi, \hat n_3)` plane gives a 2×2 linear system
    for :math:`(\cos\chi, \sin\chi)`:

    .. math::

        \begin{pmatrix} A & B \\ -B & A \end{pmatrix}
        \begin{pmatrix} \cos\chi \\ \sin\chi \end{pmatrix}
        = \begin{pmatrix} C \\ E \end{pmatrix}

    where
    :math:`A = \hat n_\phi \cdot q`,
    :math:`B = \hat n_\phi \cdot (\hat n_\chi \times q)`,
    :math:`C = \hat n_\phi \cdot v`,
    :math:`E = \hat n_3 \cdot v`.
    This yields a unique :math:`\chi` via ``atan2``.

    Phi is then solved from the :math:`\hat n_\chi` component:

    .. math::

        \hat n_\chi \cdot R(\hat n_\phi, -\phi) \cdot v = \hat n_\chi \cdot q

    which gives :math:`\alpha\cos\phi + \beta\sin\phi = \gamma`, solved
    with the standard ``atan2`` / ``acos`` decomposition.

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

    # Known vectors
    Z_prefix = ctx._cached_Z_prefix
    D = ctx._cached_D if ctx._cached_D is not None else np.eye(3)

    # Q_lab: scattering vector in the lab frame (all detector & bisect angles fixed)
    v = ctx.two_pi_over_lambda * (D @ ctx.y_eff - ctx.y_eff)  # Q_lab

    # q = Z_prefix @ Q_phi_target
    q = Z_prefix @ Q_phi_target

    # --- Step 1: Solve for chi (two solutions) ---
    #
    # Define w = R(n_phi, -phi) @ v.  Then R(n_chi, -chi) @ w = q.
    # Since R(n_phi, -phi) preserves the n_phi component:
    #   n_phi · w = n_phi · v  (constant L)
    #
    # From w = R(n_chi, chi) @ q (inverse of R(n_chi, -chi) @ w = q):
    #   n_phi · [R(n_chi, chi) @ q] = L
    #
    # Expanding via Rodrigues (n_chi ⊥ n_phi):
    #   cos(chi)*(n_phi·q) + sin(chi)*(n_phi·(n_chi×q)) = n_phi·v
    #
    # This is A*cos(chi) + B*sin(chi) = C, yielding two chi values.

    A = float(np.dot(n_phi, q))
    B = float(np.dot(n_phi, np.cross(n_chi, q)))
    C = float(np.dot(n_phi, v))

    R_chi_amp = math.sqrt(A * A + B * B)
    if R_chi_amp < 1e-12:
        return []  # degenerate: q ∥ n_chi

    cos_arg_chi = C / R_chi_amp
    if abs(cos_arg_chi) > 1.0 + 1e-8:  # pragma: no cover
        return []  # no solution
    cos_arg_chi = max(-1.0, min(1.0, cos_arg_chi))

    chi0 = math.atan2(B, A)
    delta_chi = math.acos(cos_arg_chi)

    chi_candidates_rad = [chi0 + delta_chi, chi0 - delta_chi]

    # --- Step 2: For each chi, solve for phi (unique) ---
    #
    # From w = R(n_phi, -phi) @ v and n_chi · w = n_chi · q:
    #   n_chi · [R(n_phi, -phi) @ v] = n_chi · q
    #
    # Expanding via Rodrigues (n_chi ⊥ n_phi):
    #   cos(phi)*(n_chi·v) - sin(phi)*(n_chi·(n_phi×v)) = n_chi · q
    #
    # But this gives only one equation — not enough for a unique phi.
    # We also have the n3 component (n3 = n_chi × n_phi):
    #   n3 · [R(n_phi, -phi) @ v] = n3 · [R(n_chi, chi) @ q]
    #
    # The RHS depends on chi (already known), the LHS depends on phi.
    # Together with the n_chi equation, this gives a 2×2 system for
    # (cos(phi), sin(phi)), yielding a unique phi per chi.
    #
    # n_chi equation:  (n_chi·v)*cos(phi) - (n_chi·(n_phi×v))*sin(phi) = n_chi·q
    # n3 equation:     (n3·v)*cos(phi) - (n3·(n_phi×v))*sin(phi) = n3·w_chi
    #
    # Using triple product identities (n_chi ⊥ n_phi, n3 = n_chi × n_phi):
    #   n_chi·(n_phi×v) = (n_chi×n_phi)·v = n3·v
    #   n3·(n_phi×v)    = v·(n3×n_phi) = v·n_chi   [since n3×n_phi = n_chi]
    #
    # So the system is:
    #   (n_chi·v)*cos(phi) - (n3·v)*sin(phi) = n_chi·q       ...(i)
    #   (n3·v)*cos(phi) + (n_chi·v)*sin(phi) = n3·w_chi      ...(ii)
    #
    # Wait, sign: n3·(n_phi×v) = v·(n3×n_phi). And n3×n_phi:
    #   n3 = n_chi×n_phi, so n3×n_phi = (n_chi×n_phi)×n_phi
    #   = n_chi*(n_phi·n_phi) - n_phi*(n_chi·n_phi) = n_chi
    # So n3·(n_phi×v) = v·n_chi = n_chi·v.
    #
    # n3 LHS: R(n_phi, -phi)@v dotted with n3:
    #   cos(phi)*(n3·v) + (1-cos(phi))*(n_phi·v)*(n_phi·n3) - sin(phi)*(n3·(n_phi×v))
    #   = cos(phi)*(n3·v) - sin(phi)*(n_chi·v)     [since n_phi·n3=0]
    #
    # n3 RHS: R(n_chi, chi)@q dotted with n3:
    #   cos(chi)*(n3·q) + sin(chi)*(n3·(n_chi×q))
    #   n3·(n_chi×q) = q·(n3×n_chi) = q·(-n_phi) = -(n_phi·q)
    #   [since n3×n_chi = (n_chi×n_phi)×n_chi = n_phi]
    #   Wait: BAC-CAB: (n_chi×n_phi)×n_chi = n_phi(n_chi·n_chi) - n_chi(n_phi·n_chi) = n_phi
    #   So n3×n_chi = n_phi, and n3·(n_chi×q) = q·(n3×n_chi) = q·n_phi = n_phi·q
    #   Hmm, let me redo: a·(b×c) = det(a,b,c) = c·(a×b).
    #   n3·(n_chi×q) = q·(n3×n_chi).
    #   n3×n_chi = (n_chi×n_phi)×n_chi = n_phi (BAC-CAB, orthonormal)
    #   So n3·(n_chi×q) = q·n_phi = A (= n_phi·q)
    #
    # So n3 RHS = cos(chi)*(n3·q) + sin(chi)*A
    # And n3 LHS = cos(phi)*(n3·v) - sin(phi)*(n_chi·v)
    #
    # The 2×2 system:
    #   [n_chi·v, -(n3·v)] [cos(phi)]   [n_chi·q        ]
    #   [n3·v,    n_chi·v ] [sin(phi)] = [n3·R_chi(chi)@q]
    #
    # Determinant = (n_chi·v)² + (n3·v)² = |v_perp|² (component of v ⊥ n_phi)
    # This is non-zero unless v ∥ n_phi (degenerate).

    n3 = np.cross(n_chi, n_phi)
    nchi_v = float(np.dot(n_chi, v))
    n3_v = float(np.dot(n3, v))
    nchi_q = float(np.dot(n_chi, q))

    det_phi = nchi_v * nchi_v + n3_v * n3_v
    if det_phi < 1e-20:  # pragma: no cover
        return []  # v ∥ n_phi — phi indeterminate

    # Pre-compute n3·q and n_phi·q for the chi-dependent RHS
    n3_q = float(np.dot(n3, q))
    # A = n_phi · q (already computed above)

    raw_candidates = []
    for chi_rad in chi_candidates_rad:
        chi_d = math.degrees(chi_rad)
        c_chi = math.cos(chi_rad)
        s_chi = math.sin(chi_rad)

        # n3 · R(n_chi, chi) @ q = cos(chi)*(n3·q) + sin(chi)*(n_phi·q)
        rhs_n3 = c_chi * n3_q + s_chi * A

        # Solve 2×2 system for (cos_phi, sin_phi):
        #   nchi_v * cos_phi - n3_v * sin_phi = nchi_q
        #   n3_v * cos_phi + nchi_v * sin_phi = rhs_n3
        cos_phi = (nchi_v * nchi_q + n3_v * rhs_n3) / det_phi
        sin_phi = (nchi_v * rhs_n3 - n3_v * nchi_q) / det_phi
        phi_d = math.degrees(math.atan2(sin_phi, cos_phi))

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
        if residual < 1e-6:
            # Normalize to (-180, 180]
            chi_d = (chi_d + 180.0) % 360.0 - 180.0
            phi_d = (phi_d + 180.0) % 360.0 - 180.0
            validated.append((chi_d, phi_d))

    # Sort so the "positive chi" branch (chi ≥ 0, or closer to 0°) comes
    # first.  This matches the ordering convention of the Newton solver and
    # ensures consistent branch selection in trajectory scans.
    validated.sort(key=lambda pair: -pair[0])  # descending chi → positive first

    return validated


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

    When the two free stages have orthogonal axes (standard Eulerian
    geometry), an analytic fast path (``_solve_bisecting_analytic``) is
    used instead of Newton iteration, yielding a 5-10× speedup.  The
    Newton solver is retained as fallback for non-standard (e.g. kappa)
    axis configurations.

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
    else:
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
    Seeds from 4 quadrants and deduplicates.

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

    solutions = []
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
    Kappa virtual-angle solver.

    Delegates the geometric computation to :func:`~kappa.solve_kappa_virtual`,
    then applies cut-points, limits, and deduplication.
    """
    from .kappa import solve_kappa_virtual
    from .mode import DetectorConstraint
    from .mode import SampleConstraint

    # Get (komega, kappa, kphi, det_name, ttheta, other_constraints) raw tuples
    raw = solve_kappa_virtual(geometry, Q_phi, ttheta_deg, mode)
    if not raw:  # pragma: no cover
        return []

    from .kappa import KAPPA_VIRTUAL_ANGLES
    from .kappa import kappa_to_eulerian

    # Extract non-virtual other constraints for applying to angles dict
    other_constraints = [
        c
        for c in mode.constraints
        if not (isinstance(c, SampleConstraint) and c.name in KAPPA_VIRTUAL_ANGLES)
    ]

    # Identify kappa stage names for virtual angle validation
    sample_stages = geometry.sample_stages
    kappa_idx = next(
        (i for i, s in enumerate(sample_stages) if s.name == "kappa"), None
    )
    komega_name = sample_stages[kappa_idx - 1].name if kappa_idx else None
    kphi_name = sample_stages[kappa_idx + 1].name if kappa_idx else None
    alpha_deg = geometry.kappa_alpha_deg

    # Virtual angle constraints that must be verified
    virtual_constraints = [
        c
        for c in mode.constraints
        if isinstance(c, SampleConstraint) and c.name in KAPPA_VIRTUAL_ANGLES
    ]

    solutions = []
    for angle_dict in raw:
        angles = dict(angle_dict)

        # Apply non-virtual constraints (e.g. DetectorConstraint on kappa6c modes)
        for c in other_constraints:
            if isinstance(c, SampleConstraint):  # pragma: no cover
                if c.name in geometry._stages:  # noqa: SLF001
                    angles[c.name] = float(c.value)
            elif isinstance(c, DetectorConstraint) and not c.is_qaz:  # pragma: no cover
                angles[c.name] = c.value

        # Normalise kappa angles to (-180, 180] to satisfy typical stage limits
        if komega_name and kphi_name:  # pragma: no branch
            for kname in (komega_name, "kappa", kphi_name):
                if kname in angles:  # pragma: no branch
                    a = angles[kname]
                    a = (a + 180.0) % 360.0 - 180.0
                    angles[kname] = a

        # Validate virtual angle constraints using kappa_to_eulerian
        if (
            virtual_constraints and komega_name and kphi_name and alpha_deg is not None
        ):  # pragma: no branch
            try:
                om_v, chi_v, phi_v = kappa_to_eulerian(
                    angles[komega_name],
                    angles["kappa"],
                    angles[kphi_name],
                    alpha_deg=alpha_deg,
                )
                virtual_vals = {"omega": om_v, "chi": chi_v, "phi": phi_v}
                valid = all(
                    abs(virtual_vals[c.name] - float(c.value)) < 1e-2
                    for c in virtual_constraints
                )
                if not valid:
                    continue
            except (ValueError, KeyError):  # pragma: no cover
                continue

        # Verify geometric consistency: Q_computed must match Q_target
        kv_ctx = ForwardContext(geometry)
        Q_computed = kv_ctx.q_phi_uncached(angles)
        if not np.allclose(Q_computed, Q_phi, atol=1e-3):
            continue

        _apply_cut_points(angles, mode, geometry)
        if not _check_limits(geometry, angles):
            continue  # pragma: no cover

        duplicate = False
        for existing in solutions:
            if all(
                abs(existing.get(kk, 0) - angles.get(kk, 0)) < 1e-4 for kk in angles
            ):
                duplicate = True
                break
        if not duplicate:
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
        Must have ``sample.UB`` and ``azimuthal_reference`` set.
    Q_phi : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (``UB @ hkl``).

    Returns
    -------
    float or None
        ψ in degrees (−180°, +180°], or ``None`` when ψ is undefined
        (Q ∥ incident beam, or reference vector ∥ Q).
    """
    n_hkl = geometry.azimuthal_reference
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

    if geometry.azimuthal_reference is None:
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
    3. If they disagree (beyond 0.1° tolerance): returns ``[]``.
    4. If they agree: delegates to the appropriate existing solver
       (bisecting, kappa-virtual, or synthetic bisecting) and returns
       all solutions.

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
        return []  # ψ undefined for this reflection

    # Compare natural psi with target (tolerance 0.1° — generous enough
    # to handle float rounding, tight enough to be physically meaningful)
    diff = abs(natural_psi - psi_target)
    # Handle wraparound: e.g. -179.9 vs 180.1
    if diff > 180.0:
        diff = 360.0 - diff
    if diff > 0.1:
        return []  # this (h,k,l) is not accessible at the stored ψ

    # ψ is satisfied — delegate to the appropriate existing solver.
    # The psi constraint is automatically satisfied by ALL Bragg solutions.

    if mode.has_bisect:
        # psic, kappa6c: mode already has a BisectConstraint
        return _solve_bisecting(geometry, Q_phi, ttheta_deg, mode)

    if _is_kappa_virtual_mode(geometry, mode):  # pragma: no cover
        # kappa4cv, kappa4ch: kappa virtual angle mode
        return _solve_kappa_virtual(geometry, Q_phi, ttheta_deg, mode)

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
        """Compute sample rotation matrix from angle values."""
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

        # Pre-compute the Z_prefix for stages before the outer stage
        Z_pre_outer = np.eye(3)
        for i in range(outer_stage_idx):
            s = sample_stages[i]
            angle = angles_base.get(s.name, s.angle)
            Z_pre_outer = _rotation_matrix_normalized(s._axis_hat, angle) @ Z_pre_outer  # noqa: SLF001

    def _solve_inner_eulerian_fast(
        outer_deg: float,
    ) -> list[tuple[float, float]]:
        """
        Fast analytic (chi, phi) solve for Eulerian geometries.

        Returns list of (chi_deg, phi_deg) pairs.  No ForwardContext created.

        When the decomposition is degenerate (q parallel to n_chi, so chi
        is indeterminate), returns an empty list.  The caller handles
        degenerate points via `_degenerate_outers`.
        """
        # Z_prefix = R_outer @ Z_pre_outer
        R_outer = _rotation_matrix_normalized(
            outer_stage_obj._axis_hat,
            outer_deg,  # noqa: SLF001
        )
        Z_prefix = R_outer @ Z_pre_outer

        # q = Z_prefix @ Q_phi_target
        q = Z_prefix @ Q_phi

        # --- Solve for chi: A*cos(chi) + B*sin(chi) = C ---
        A = float(np.dot(n_phi, q))
        B = float(np.dot(n_phi, np.cross(n_chi, q)))
        C = float(np.dot(n_phi, v))

        R_chi_amp = math.sqrt(A * A + B * B)
        if R_chi_amp < 1e-12:
            return []  # degenerate — handled separately

        cos_arg_chi = C / R_chi_amp
        if abs(cos_arg_chi) > 1.0 + 1e-8:  # pragma: no cover
            return []  # pragma: no cover
        cos_arg_chi = max(-1.0, min(1.0, cos_arg_chi))

        chi0 = math.atan2(B, A)
        delta_chi = math.acos(cos_arg_chi)
        chi_candidates_rad = [chi0 + delta_chi, chi0 - delta_chi]

        # --- For each chi, solve for phi (2x2 system) ---
        nchi_v = float(np.dot(n_chi, v))
        n3_v = float(np.dot(n3, v))
        nchi_q = float(np.dot(n_chi, q))
        det_phi = nchi_v * nchi_v + n3_v * n3_v
        if det_phi < 1e-20:  # pragma: no cover
            return []  # pragma: no cover

        n3_q = float(np.dot(n3, q))

        results = []
        for chi_rad in chi_candidates_rad:
            c_chi = math.cos(chi_rad)
            s_chi = math.sin(chi_rad)
            rhs_n3 = c_chi * n3_q + s_chi * A
            cos_phi = (nchi_v * nchi_q + n3_v * rhs_n3) / det_phi
            sin_phi = (nchi_v * rhs_n3 - n3_v * nchi_q) / det_phi
            chi_d = math.degrees(chi_rad)
            phi_d = math.degrees(math.atan2(sin_phi, cos_phi))
            # Normalize to (-180, 180]
            chi_d = (chi_d + 180.0) % 360.0 - 180.0
            phi_d = (phi_d + 180.0) % 360.0 - 180.0
            results.append((chi_d, phi_d))

        return results

    def _find_degenerate_outers() -> list[float]:
        """
        Find outer angles where the analytic decomposition is degenerate.

        Degeneracy occurs when q = Z_prefix @ Q_phi is parallel to n_chi,
        meaning R_chi_amp = 0.  At these points chi is indeterminate and
        solutions (if they exist) must be found by scanning chi.
        """
        degenerate = []
        for i in range(720):
            outer_deg = -180.0 + i * 0.5
            R_outer = _rotation_matrix_normalized(
                outer_stage_obj._axis_hat,
                outer_deg,  # noqa: SLF001
            )
            q = R_outer @ Z_pre_outer @ Q_phi
            A = float(np.dot(n_phi, q))
            B = float(np.dot(n_phi, np.cross(n_chi, q)))
            if math.sqrt(A * A + B * B) < 0.1:
                degenerate.append(outer_deg)
        return degenerate

    def _solve_degenerate_outer(
        outer_deg: float,
    ) -> list[dict[str, float]]:
        """
        At degenerate outer angles (q parallel to n_chi), scan chi to find
        solutions that satisfy both the Bragg and Ewald conditions.

        At these points, chi is free and phi is determined by chi.  We scan
        chi, compute phi from the remaining Bragg components, and detect
        Ewald sign changes between adjacent chi values.  A bisection then
        refines the chi to the exact root.
        """
        angles = dict(angles_base)
        angles[outer_stage_name] = outer_deg

        ctx = ForwardContext(geometry)
        ctx.prepare_caching(angles, {chi_stage_name, phi_stage_name})

        R_outer = _rotation_matrix_normalized(
            outer_stage_obj._axis_hat,
            outer_deg,  # noqa: SLF001
        )
        q = R_outer @ Z_pre_outer @ Q_phi
        nchi_v = float(np.dot(n_chi, v))
        n3_v = float(np.dot(n3, v))
        nchi_q = float(np.dot(n_chi, q))
        n3_q = float(np.dot(n3, q))
        A_local = float(np.dot(n_phi, q))
        det_phi = nchi_v * nchi_v + n3_v * n3_v
        if det_phi < 1e-20:  # pragma: no cover
            return []  # pragma: no cover

        def _chi_to_trial(chi_deg_f: float) -> dict[str, float] | None:
            """Given a chi value, compute phi and return a trial dict."""
            chi_rad = math.radians(chi_deg_f)
            c_chi = math.cos(chi_rad)
            s_chi = math.sin(chi_rad)
            rhs_n3 = c_chi * n3_q + s_chi * A_local
            cos_phi = (nchi_v * nchi_q + n3_v * rhs_n3) / det_phi
            sin_phi = (nchi_v * rhs_n3 - n3_v * nchi_q) / det_phi
            phi_d = math.degrees(math.atan2(sin_phi, cos_phi))

            trial = dict(angles)
            trial[chi_stage_name] = chi_deg_f
            trial[phi_stage_name] = phi_d

            # Verify Bragg condition
            Q_computed = ctx.q_phi(trial)
            if float(np.linalg.norm(Q_computed - Q_phi)) > 1e-3:
                return None
            return trial

        # Scan chi at 2-degree intervals, detect sign changes in Ewald
        candidates = []
        prev_chi: float | None = None
        prev_ew: float = 0.0

        for chi_int in range(-180, 180, 2):
            chi_f = float(chi_int)
            trial = _chi_to_trial(chi_f)
            if trial is None:
                prev_chi = None
                continue
            ew = _ewald_residual(trial)

            if abs(ew) < 1e-3:  # pragma: no cover
                candidates.append(trial)  # pragma: no cover
            elif prev_chi is not None and prev_ew * ew < 0:
                # Sign change: bisect chi to find root
                lo, hi = prev_chi, chi_f
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

            prev_chi = chi_f
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
# Surface diffraction solvers (ReferenceConstraint modes)
# ---------------------------------------------------------------------------


def _is_surface_mode(geometry: AdHocDiffractometer, mode) -> bool:
    """Return True when mode has a surface ReferenceConstraint and surface_normal is set.

    The ``"psi"`` ReferenceConstraint is NOT a surface mode — it is handled
    by :func:`_is_psi_mode` / :func:`_solve_psi_mode` instead.
    """
    from .mode import ReferenceConstraint

    return geometry.surface_normal is not None and any(
        isinstance(c, ReferenceConstraint) and c.name != "psi" for c in mode.constraints
    )


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

    # Set the active detector stage to ttheta_deg
    det_stage = geometry.detector_stages[-1]
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
    else:  # pragma: no cover — no preset currently uses non-90° qaz
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
