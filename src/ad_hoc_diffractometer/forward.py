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
BisectingMode (four-sample-stage geometries with one detector arm)
    Classic Eulerian four-circle bisecting solution (Busing & Levy 1967,
    section "Angle settings").  Valid for psic (eta/delta bisecting),
    fourcv, fourch (omega/ttheta bisecting), and analogous geometries.

    Algorithm:
        1. Q_phi = UB @ (h, k, l)             — target in phi frame
        2. \|Q\| → ttheta via Bragg's law
        3. frozen stages set from mode.frozen_angles
        4. detector_stage = ttheta (computed from \|Q\|)
        5. sample_stage   = ttheta / 2 (bisecting)
        6. Remaining free sample stages (chi, phi or kchi, kphi) solved
           from the direction of Q_phi, choosing among the standard
           solution branches (two solutions: chi in [0°,180°] and
           chi in [-180°, 0°]).

FixedAngleMode
    The named stage is frozen at the stored value.  The remaining
    geometry-specific solver is then called with the reduced set of
    free stages.  Currently only supported when the fixed stage is
    one of the "outer" sample stages (not the detector or the bisected
    stage).

Raises
------
NotImplementedError
    If the active mode is None or is not a supported type for this geometry.
ValueError
    If wavelength or UB matrix is not set.
ValueError
    If (h, k, l) = (0, 0, 0).
ValueError
    If the requested \|Q\| exceeds the Ewald sphere (wavelength too long).

References
----------
Busing & Levy, Acta Cryst. 22, 457-464 (1967) — angle-setting equations
You, J. Appl. Cryst. 32, 614-623 (1999) — psic geometry
"""

from __future__ import annotations

import logging
import math
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .geometry import AdHocDiffractometer

logger = logging.getLogger(__name__)

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
    ValueError
        If the scattering vector magnitude exceeds the Ewald sphere
        (\|Q\| > 4π/λ, i.e. Bragg condition cannot be satisfied).
    NotImplementedError
        If no active mode is set or the active mode type is not supported
        for this geometry.
    """
    from .mode import BisectingMode
    from .mode import FixedAngleMode

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
        raise ValueError(
            f"forward(): |Q| = {Q_mag:.6g} Å⁻¹ exceeds the Ewald sphere "
            f"(max = {4.0 * math.pi / wavelength:.6g} Å⁻¹) at "
            f"λ = {wavelength} Å.  The reflection cannot be reached."
        )

    ttheta_rad = 2.0 * math.asin(sin_theta)
    ttheta_deg = math.degrees(ttheta_rad)

    # --- Dispatch to mode-specific solver ------------------------------------

    mode = geometry.mode

    if isinstance(mode, BisectingMode):
        return _solve_bisecting(geometry, Q_phi, ttheta_deg, mode)

    if isinstance(mode, FixedAngleMode):
        return _solve_fixed_angle(geometry, Q_phi, ttheta_deg, mode)

    raise NotImplementedError(
        f"forward(): mode {geometry.mode_name!r} "
        f"({type(mode).__name__}) is not supported for geometry "
        f"{geometry.name!r}.  Supported mode types: "
        "BisectingMode, FixedAngleMode."
    )


# ---------------------------------------------------------------------------
# Mode-specific solvers
# ---------------------------------------------------------------------------


def _solve_bisecting(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Bisecting-mode forward solver.

    Implements the Busing & Levy (1967) angle-setting algorithm for an
    Eulerian geometry in the bisecting condition.

    The detector stage is set to ``ttheta_deg``.
    The sample_stage (omega/eta) is set to ``ttheta_deg / 2``.
    The remaining two sample stages (chi, phi) are computed from the
    direction of Q_phi.

    Two solutions are returned corresponding to:
        - "positive chi" branch: chi ∈ [0°, 180°]
        - "negative chi" branch: chi ∈ (-180°, 0°)

    Both solutions are validated against stage limits; those that fail are
    dropped.  Cut-points from the mode are applied to each angle.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
        Target scattering vector in the phi frame (Å⁻¹).
    ttheta_deg : float
        Detector angle (2θ) in degrees.
    mode : BisectingMode

    Returns
    -------
    list of dict[str, float]
        Valid solutions (may be empty if all are out of limits).
    """
    # Identify the three free sample stages (beyond the bisected one).
    # For standard four-circle: [omega/eta, chi, phi].
    # For psic with frozen mu/nu:  [eta, chi, phi].
    sample_stages = geometry.sample_stages

    # Build the baseline angle dict: all stages at their current angles.
    angles: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }

    # Apply frozen angles from the mode.
    for name, val in mode.frozen_angles.items():
        angles[name] = val

    # Set the detector stage (first/outermost detector stage = ttheta/delta).
    detector_name = mode.detector_stage
    angles[detector_name] = ttheta_deg
    # Any additional detector stages (e.g. nu in psic) keep frozen or current.

    # The bisecting sample stage = ttheta / 2.
    sample_bisect_name = mode.sample_stage
    angles[sample_bisect_name] = ttheta_deg / 2.0

    # The remaining two free sample stages are the ones not frozen and not bisected.
    constrained = set(mode.constrained_stages)
    constrained.add(detector_name)
    free_sample = [s for s in sample_stages if s.name not in constrained]

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

    # With all constrained angles fixed, we need to find (chi, phi) such that
    # angles_to_phi_vector(geometry, **all_angles) == Q_phi.
    #
    # This is a 2D root-finding problem.  We use a numerical solver rather
    # than an analytic formula so that the result is correct for any basis
    # convention and stage axis handedness.
    #
    # Two solution branches are seeded from the analytic decomposition of
    # Q_phi relative to the chi rotation axis, giving starting guesses that
    # are typically within a few degrees of the solution.

    from .orientation import angles_to_phi_vector as _a2phi

    Q_norm = float(np.linalg.norm(Q_phi))
    q_hat = Q_phi / Q_norm

    # Normalised chi axis vector
    chi_ax = chi_stage.axis / np.linalg.norm(chi_stage.axis)

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

    # Seed (chi, phi) pairs.  Use a combination of analytically motivated
    # seeds and a coarse grid to ensure both solution branches are found even
    # in degenerate cases (Q along chi axis, or chi ≈ 0 / 90 / 180°).
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

    for chi_start, phi_start in seed_pairs:
        sol_angles = _solve_two_angles(
            geometry=geometry,
            fixed_angles=angles,
            free_stage_1=chi_stage.name,
            free_stage_2=phi_stage.name,
            start_1=chi_start,
            start_2=phi_start,
            Q_phi_target=Q_phi,
            a2phi_fn=_a2phi,
        )
        if sol_angles is None:  # pragma: no branch
            continue  # pragma: no cover

        sol = dict(angles)
        sol[chi_stage.name] = sol_angles[0]
        sol[phi_stage.name] = sol_angles[1]

        _apply_cut_points(sol, mode, geometry)

        # De-duplicate: skip if this solution is the same as an existing one
        duplicate = False
        for existing in solutions:
            if all(abs(existing.get(k, 0) - sol.get(k, 0)) < 1e-4 for k in sol):
                duplicate = True
                break

        if not duplicate and _check_limits(geometry, sol):
            solutions.append(sol)

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
        ``angles_to_phi_vector`` function.
    max_iter, tol : int, float
        Convergence parameters.
    _one_dimensional : bool
        If True, solve for only ``free_stage_1`` (1D problem).

    Returns
    -------
    (angle_1, angle_2) : tuple of float, or None if no convergence.
        When ``_one_dimensional=True``, ``angle_2`` equals ``start_2``.
    """
    n_free = 1 if _one_dimensional else 2
    x = np.array([start_1, start_2][:n_free], dtype=float)

    def residual(x):
        trial = dict(fixed_angles)
        trial[free_stage_1] = float(x[0])
        if not _one_dimensional:
            trial[free_stage_2] = float(x[1])
        Q = a2phi_fn(geometry, **trial)
        return Q - Q_phi_target

    def jacobian_fd(x, r0, h=1e-4):
        """Finite-difference Jacobian (3 × n_free)."""
        J = np.zeros((3, n_free))
        for i in range(n_free):
            xp = x.copy()
            xp[i] += h
            J[:, i] = (residual(xp) - r0) / h
        return J

    for _ in range(max_iter):
        r = residual(x)
        if np.linalg.norm(r) < tol:
            a1 = float(x[0])
            a2 = float(x[1]) if not _one_dimensional else start_2
            return (a1, a2)
        J = jacobian_fd(x, r)
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


def _solve_fixed_angle(
    geometry: AdHocDiffractometer,
    Q_phi: np.ndarray,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    FixedAngleMode forward solver.

    Freezes the named stage at its fixed value, then delegates to the
    bisecting solver treating the remaining stages as free — provided the
    geometry has a bisecting-capable structure (detector + 3 sample stages
    where one is the equivalent of omega/eta).

    This covers the common case of fixed_chi (chi frozen at e.g. 90°,
    omega/ttheta bisecting, phi free) or fixed_mu (mu=0, rest as psic
    bisecting).

    For each frozen-stage value, the remaining angles are solved via the
    bisecting equations restricted to the appropriate stages.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Q_phi : numpy.ndarray, shape (3,)
    ttheta_deg : float
    mode : FixedAngleMode

    Returns
    -------
    list of dict[str, float]
    """
    from .mode import BisectingMode

    # Find the geometry's default bisecting mode to re-use its stage names.
    bisecting_mode = None
    for m in geometry.modes.values():
        if isinstance(m, BisectingMode):
            bisecting_mode = m
            break

    if bisecting_mode is None:
        raise NotImplementedError(
            f"forward(): FixedAngleMode solver for geometry {geometry.name!r} "
            "requires a BisectingMode to be defined in the geometry's modes "
            "to identify the bisecting and detector stage names.  "
            "No BisectingMode found."
        )

    # Build an effective BisectingMode that includes the fixed stage as frozen.
    effective_frozen = dict(bisecting_mode.frozen_angles)
    effective_frozen.update(mode.frozen_angles)

    effective_mode = BisectingMode(
        sample_stage=bisecting_mode.sample_stage,
        detector_stage=bisecting_mode.detector_stage,
        frozen_angles=effective_frozen,
        cut_points=mode.cut_points or bisecting_mode.cut_points or None,
    )

    return _solve_bisecting(geometry, Q_phi, ttheta_deg, effective_mode)


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
