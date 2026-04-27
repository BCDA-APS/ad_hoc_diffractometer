# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
kappa.py — Kappa-to-Eulerian angle conversion.

Kappa diffractometers replace the Eulerian chi circle with a kappa arm
tilted at angle alpha_0 from the omega axis.  The real kappa motor angles
(komega, kappa, kphi) map to virtual Eulerian pseudoangles (omega, chi, phi)
via the relations in Walko (2016), eq. [16]:

    chi   = 2 arcsin[sin(kappa/2) · sin(alpha_0)]
    offset = arccos[cos(kappa/2) / cos(chi/2)]
    omega = komega − offset
    phi   = kphi   − offset

Both directions are implemented here as pure-NumPy functions.

Functions
---------
:func:`kappa_to_eulerian`
    Forward direction: real (komega, kappa, kphi) → virtual (omega, chi, phi).

:func:`eulerian_to_kappa`
    Inverse direction: virtual (omega, chi, phi) → real (komega, kappa, kphi).

Notes
-----
The default kappa tilt angle is 50° (Walko 2016; Enraf-Nonius convention;
ITC Vol. C §2.2.6).

The inverse has two solutions distinguished by the sign of kappa (the *branch*
parameter).  The positive branch (kappa ≥ 0) is the default; the negative
branch gives kappa ≤ 0.  The caller should choose the branch that keeps all
angles within the hardware limits.

Singularity: when ``|chi|`` approaches ±2·alpha_0, ``|kappa|`` approaches 180° and
the offset formula becomes numerically ill-conditioned.  A :class:`ValueError`
is raised when ``|chi|`` ≥ 2·alpha_0 − epsilon.

References
----------
* Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), eq. [16].
* ITC Vol. C §2.2.6 (2006).
"""

from __future__ import annotations

import itertools
import math
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

#: Virtual Eulerian pseudoangle names on kappa geometries.
KAPPA_VIRTUAL_ANGLES: frozenset[str] = frozenset({"omega", "chi", "phi"})


def kappa_to_eulerian(
    komega: float,
    kappa: float,
    kphi: float,
    alpha_deg: float = 50.0,
) -> tuple[float, float, float]:
    """
    Convert real kappa angles to virtual Eulerian pseudoangles.

    Implements Walko (2016), eq. [16] — forward direction.

    Parameters
    ----------
    komega : float
        Real kappa-omega angle in degrees.
    kappa : float
        Real kappa angle in degrees.
    kphi : float
        Real kappa-phi angle in degrees.
    alpha_deg : float, optional
        Kappa tilt angle in degrees.  Default 50°.

    Returns
    -------
    omega, chi, phi : tuple of float
        Virtual Eulerian pseudoangles in degrees.

    Raises
    ------
    ValueError
        If the combination of kappa and alpha_deg places chi outside
        the reachable range (``|kappa/2|`` · sin(alpha_0) > 1).

    Examples
    --------
    >>> kappa_to_eulerian(57.045, 134.756, 57.045, alpha_deg=50.0)
    (-0.0, 90.0, -0.0)
    """
    a0 = math.radians(alpha_deg)
    k = math.radians(kappa)

    sin_chi_half = math.sin(k / 2.0) * math.sin(a0)
    if abs(sin_chi_half) > 1.0:  # pragma: no cover
        raise ValueError(
            f"kappa_to_eulerian: kappa={kappa:.4f}°, alpha_deg={alpha_deg:.4f}° "
            f"gives |sin(chi/2)| = {abs(sin_chi_half):.6f} > 1 — unreachable chi."
        )

    chi_half = math.asin(sin_chi_half)
    chi = math.degrees(2.0 * chi_half)

    cos_chi_half = math.cos(chi_half)
    if abs(cos_chi_half) < 1e-14:  # pragma: no cover
        raise ValueError(
            f"kappa_to_eulerian: chi approaches ±180° (chi_half={math.degrees(chi_half):.4f}°) "
            f"— singularity in the offset calculation."
        )

    cos_ratio = math.cos(k / 2.0) / cos_chi_half
    cos_ratio = max(-1.0, min(1.0, cos_ratio))  # guard against float rounding
    offset = math.degrees(math.acos(cos_ratio))
    # The sign of offset matches the sign of kappa
    if kappa < 0.0:
        offset = -offset

    omega = komega - offset
    phi = kphi - offset

    return omega, chi, phi


def eulerian_to_kappa(
    omega: float,
    chi: float,
    phi: float,
    alpha_deg: float = 50.0,
    branch: int = +1,
) -> tuple[float, float, float]:
    """
    Convert virtual Eulerian pseudoangles to real kappa angles.

    Inverse of Walko (2016), eq. [16].

    Parameters
    ----------
    omega : float
        Virtual Eulerian omega in degrees.
    chi : float
        Virtual Eulerian chi in degrees.
    phi : float
        Virtual Eulerian phi in degrees.
    alpha_deg : float, optional
        Kappa tilt angle in degrees.  Default 50°.
    branch : {+1, -1}, optional
        Branch selection.  ``+1`` gives kappa ≥ 0 (default);
        ``-1`` gives kappa ≤ 0.

    Returns
    -------
    komega, kappa, kphi : tuple of float
        Real kappa motor angles in degrees.

    Raises
    ------
    ValueError
        If ``|chi|`` ≥ 2·alpha_deg (chi outside the reachable range).
    ValueError
        If branch is not +1 or -1.

    Examples
    --------
    >>> eulerian_to_kappa(0.0, 90.0, 0.0, alpha_deg=50.0)
    (57.04519..., 134.75594..., 57.04519...)
    >>> eulerian_to_kappa(0.0, 90.0, 0.0, alpha_deg=50.0, branch=-1)
    (-57.04519..., -134.75594..., -57.04519...)
    """
    if branch not in (+1, -1):
        raise ValueError(f"eulerian_to_kappa: branch must be +1 or -1; got {branch!r}.")

    a0 = math.radians(alpha_deg)
    chi_half = math.radians(chi / 2.0)

    chi_limit = 2.0 * alpha_deg
    if abs(chi) >= chi_limit - 1e-9:
        raise ValueError(
            f"eulerian_to_kappa: |chi| = {abs(chi):.4f}° must be less than "
            f"2·alpha_deg = {chi_limit:.4f}° (kappa geometry limit).  "
            "The requested chi is outside the reachable range."
        )

    # Use |chi/2| to get the magnitude of kappa, then apply branch for sign.
    sin_kappa_half = abs(math.sin(chi_half)) / math.sin(a0)
    # Clamp for float safety (already guarded above, but defensive)
    sin_kappa_half = max(-1.0, min(1.0, sin_kappa_half))
    kappa_half = math.asin(sin_kappa_half)
    kappa = math.degrees(2.0 * kappa_half) * branch

    # Recompute chi_half from the chosen kappa branch for consistency
    k = math.radians(kappa)
    sin_ch = math.sin(k / 2.0) * math.sin(a0)
    ch = math.asin(max(-1.0, min(1.0, sin_ch)))

    cos_ch = math.cos(ch)
    if abs(cos_ch) < 1e-14:  # pragma: no cover
        raise ValueError("eulerian_to_kappa: chi approaches ±180° — singularity.")

    cos_ratio = math.cos(k / 2.0) / cos_ch
    cos_ratio = max(-1.0, min(1.0, cos_ratio))
    offset = math.degrees(math.acos(cos_ratio))
    if kappa < 0.0:
        offset = -offset

    komega = omega + offset
    kphi = phi + offset

    return komega, kappa, kphi


# ---------------------------------------------------------------------------
# Kappa virtual-angle mode detection and solving
# ---------------------------------------------------------------------------


def is_kappa_virtual_mode(
    geometry: AdHocDiffractometer,
    mode,
) -> bool:
    """
    Return True when *mode* contains a virtual Eulerian angle constraint
    (omega, chi, or phi) and *geometry* is a kappa diffractometer.

    Used by :func:`~forward._solve_constraint_set` to detect when the
    kappa virtual-angle solver should be dispatched.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    mode : ConstraintSet

    Returns
    -------
    bool
    """
    from .mode import SampleConstraint

    if geometry.kappa_alpha_deg is None:
        return False
    stage_names = {s.name for s in geometry._stages.values()}  # noqa: SLF001
    if "kappa" not in stage_names:
        return False
    return any(
        isinstance(c, SampleConstraint) and c.name in KAPPA_VIRTUAL_ANGLES
        for c in mode.constraints
    )


def solve_kappa_virtual(
    geometry: AdHocDiffractometer,
    Q_phi,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Solve a kappa virtual-angle mode: find real kappa motor angles satisfying
    the virtual Eulerian angle constraints and the Bragg condition.

    The solver works in virtual Eulerian space by using the sample rotation
    Newton-Raphson numeric solver (2D free virtual angles):

    1. Reads the fixed virtual Eulerian angles (omega, chi, phi) from the mode.
    2. Applies fixed-angle rotations to reduce Q to the remaining free-angle
       subspace, then decomposes to find the free angles.
    3. Converts each Eulerian solution to real kappa motor angles via
       :func:`eulerian_to_kappa` (both branches ±1).

    Parameters
    ----------
    geometry : AdHocDiffractometer
        A kappa geometry with ``kappa_alpha_deg`` set.
    Q_phi : array-like, shape (3,)
        Scattering vector in the phi frame (Å⁻¹).
    ttheta_deg : float
        Detector angle (2θ) in degrees.
    mode : ConstraintSet
        Must contain at least one ``SampleConstraint`` naming a virtual
        Eulerian angle.

    Returns
    -------
    list of dict[str, float]
        One angles dict per solution, keyed by real stage names.
        Cut-points, limits, and deduplication are applied by the caller
        (:func:`~forward._solve_kappa_virtual`).
    """
    import numpy as np

    from .mode import SampleConstraint

    alpha_deg = geometry.kappa_alpha_deg
    sample_stages = geometry.sample_stages
    stage_names = [s.name for s in sample_stages]

    # Locate komega, kappa, kphi by position relative to the "kappa" stage
    kappa_idx = next(
        (i for i, s in enumerate(sample_stages) if s.name == "kappa"), None
    )
    if (
        kappa_idx is None or kappa_idx == 0 or kappa_idx >= len(stage_names) - 1
    ):  # pragma: no cover
        return []

    komega_name = stage_names[kappa_idx - 1]
    kphi_name = stage_names[kappa_idx + 1]
    det_stage = geometry.detector_stages[-1]

    # The kappa stage (middle) acts as the virtual chi axis in the Newton solver
    chi_stage_eq = sample_stages[kappa_idx]  # kappa ↔ virtual chi

    # Extract fixed virtual angle values from mode constraints
    fixed_virtual: dict[str, float] = {}
    for c in mode.constraints:
        if (
            isinstance(c, SampleConstraint) and c.name in KAPPA_VIRTUAL_ANGLES
        ):  # pragma: no branch
            fixed_virtual[c.name] = float(c.value)

    Q_phi_arr = np.asarray(Q_phi, dtype=float)
    Q_mag = float(np.linalg.norm(Q_phi_arr))
    if Q_mag < 1e-14:  # pragma: no cover
        return []
    Q_hat = Q_phi_arr / Q_mag

    fixed_omega = fixed_virtual.get("omega")
    fixed_chi = fixed_virtual.get("chi")
    fixed_phi = fixed_virtual.get("phi")

    # Strategy: parameterise the problem in terms of the two free virtual angles,
    # then use a numeric Newton-Raphson to find (chi_e, phi_e) [or whichever two are
    # free] such that the scattering condition is satisfied.
    #
    # For a given set of virtual Eulerian angles (omega_e, chi_e, phi_e), the real
    # kappa angles are determined by eulerian_to_kappa.  The scattering condition
    # is:  angles_to_phi_vector(geometry, komega=ko, kappa=k, kphi=kp, ttheta=tt)
    #      == Q_phi_target
    #
    # This is 3 equations in 2 unknowns (the two free virtual angles), but the
    # system is consistent when a solution exists.  We solve the 2D sub-system
    # (2 of the 3 components) and verify the third.

    from .forward import ForwardContext

    # Create ForwardContext for fast Q_phi computation
    _ctx = ForwardContext(geometry)
    # For kappa virtual, the free stages are komega, kappa, kphi, and the
    # detector stage (ttheta) — prepare caching for any remaining fixed stages
    _ctx.prepare_caching(
        {s.name: s.angle for s in list(geometry._stages.values())},  # noqa: SLF001
        {komega_name, "kappa", kphi_name, det_stage.name},
    )

    def _newton_solve(x0: np.ndarray, branch: int = +1) -> np.ndarray | None:
        """Newton-Raphson solver for the 2D free virtual angle problem."""
        x = x0.copy()
        h = 1e-4
        for _ in range(60):
            r = _q_residual_branch(x, branch)
            rn = float(np.linalg.norm(r[:2]))  # use first 2 components
            if rn < 1e-10:
                return x
            # Finite-difference Jacobian (2×2 from 3×2)
            J = np.zeros((2, 2))
            for j in range(2):
                xph = x.copy()
                xph[j] += h
                rph = _q_residual_branch(xph, branch)
                J[:, j] = (rph[:2] - r[:2]) / h
            try:
                dx = np.linalg.solve(J, -r[:2])
            except np.linalg.LinAlgError:
                return None
            x = x + np.clip(dx, -30.0, 30.0)
        r = _q_residual_branch(x, branch)
        return x if float(np.linalg.norm(r)) < 1e-6 else None

    def _q_residual_branch(free_vals: np.ndarray, branch: int) -> np.ndarray:
        if fixed_omega is not None:
            omega_e, chi_e, phi_e = fixed_omega, free_vals[0], free_vals[1]
        elif fixed_chi is not None:
            omega_e, chi_e, phi_e = free_vals[0], fixed_chi, free_vals[1]
        else:
            omega_e, chi_e, phi_e = free_vals[0], free_vals[1], fixed_phi
        try:
            ko, k, kp = eulerian_to_kappa(
                omega_e, chi_e, phi_e, alpha_deg=alpha_deg, branch=branch
            )
        except ValueError:  # pragma: no cover
            return np.array([1e6, 1e6, 1e6])
        trial = {s.name: s.angle for s in list(geometry._stages.values())}  # noqa: SLF001
        trial[komega_name] = ko
        trial["kappa"] = k
        trial[kphi_name] = kp
        trial[det_stage.name] = ttheta_deg
        return _ctx.q_phi(trial) - Q_phi_arr

    # Build seed points from Q geometry — chi_abs is a natural seed
    chi_ax = np.asarray(chi_stage_eq.axis, dtype=float)
    chi_ax /= np.linalg.norm(chi_ax)
    q_along_chi = float(np.clip(np.dot(Q_hat, chi_ax), -1.0, 1.0))
    chi_abs = math.degrees(math.asin(abs(q_along_chi)))

    if fixed_omega is not None:
        seeds = [
            np.array([chi_abs, 0.0]),
            np.array([-chi_abs, 0.0]),
            np.array([chi_abs, 90.0]),
            np.array([-chi_abs, 90.0]),
            np.array([chi_abs, 180.0]),
            np.array([-chi_abs, 180.0]),
            np.array([chi_abs, 270.0]),
            np.array([-chi_abs, 270.0]),
            np.array([0.0, 0.0]),
            np.array([0.0, 90.0]),
        ]
    elif fixed_chi is not None:
        seeds = [
            np.array([0.0, 0.0]),
            np.array([0.0, 90.0]),
            np.array([0.0, 180.0]),
            np.array([0.0, 270.0]),
            np.array([30.0, 0.0]),
            np.array([-30.0, 0.0]),
            np.array([60.0, 45.0]),
        ]
    else:  # fixed_phi — (omega, chi) are free
        seeds = []
        for o_s in range(-180, 180, 20):
            for c_s in [-chi_abs, 0.0, chi_abs, 45.0, 90.0, -45.0, -90.0]:
                seeds.append(np.array([float(o_s), float(c_s)]))

    solutions_eulerian: list[tuple[float, float, float]] = []
    seen: list[np.ndarray] = []

    # Early termination: once _MIN_SOLUTIONS unique Eulerian solutions have
    # been found, stop after _MAX_STALE consecutive stale seeds (failed,
    # high-residual, or duplicate).  Also stop immediately at _MAX_SOLUTIONS.
    _MAX_SOLUTIONS = 4
    _MIN_SOLUTIONS = 2
    _MAX_STALE = 6
    stale_count = 0

    def _should_stop() -> bool:
        """True when enough solutions found and seeds are producing duplicates."""
        return len(solutions_eulerian) >= _MIN_SOLUTIONS and stale_count >= _MAX_STALE

    for branch, seed in itertools.product((+1, -1), seeds):
        sol = _newton_solve(seed, branch)
        if sol is None:
            stale_count += 1
            if _should_stop():
                break  # pragma: no cover
            continue
        # Check residual
        r = _q_residual_branch(sol, branch)
        if float(np.linalg.norm(r)) > 1e-6:
            stale_count += 1
            if _should_stop():
                break
            continue
        # Reconstruct full virtual angles
        if fixed_omega is not None:
            omega_e, chi_e, phi_e = fixed_omega, float(sol[0]), float(sol[1])
        elif fixed_chi is not None:
            omega_e, chi_e, phi_e = float(sol[0]), fixed_chi, float(sol[1])
        else:
            omega_e, chi_e, phi_e = float(sol[0]), float(sol[1]), fixed_phi
        # Deduplicate
        is_dup = any(
            abs(omega_e - float(s[0])) < 1e-3
            and abs(chi_e - float(s[1])) < 1e-3
            and abs(phi_e - float(s[2])) < 1e-3
            for s in seen
        )
        if is_dup:
            stale_count += 1
            if _should_stop():
                break  # pragma: no cover
            continue
        seen.append(np.array([omega_e, chi_e, phi_e]))
        solutions_eulerian.append((omega_e, chi_e, phi_e))
        stale_count = 0  # reset: found a new unique solution
        if len(solutions_eulerian) >= _MAX_SOLUTIONS:
            break

    if not solutions_eulerian:  # pragma: no cover
        return []

    # Build the baseline angles dict
    base_angles: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }
    base_angles[det_stage.name] = ttheta_deg

    results: list[dict[str, float]] = []
    for omega_e, chi_e, phi_e in solutions_eulerian:
        for branch in (+1, -1):
            try:
                ko, k, kp = eulerian_to_kappa(
                    omega_e, chi_e, phi_e, alpha_deg=alpha_deg, branch=branch
                )
            except ValueError:  # pragma: no cover
                continue
            angles = dict(base_angles)
            angles[komega_name] = ko
            angles["kappa"] = k
            angles[kphi_name] = kp
            results.append(angles)

    return results
