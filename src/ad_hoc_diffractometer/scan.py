# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
scan.py — Reciprocal-space trajectory computation.

Functions
---------
:func:`hkl_trajectory` ``(geometry, trajectory, n_points, solution_key=NEAREST_ANGLES)``
    **Generator.** Yield motor-angle dicts along a reciprocal-space path.
    Use ``list(hkl_trajectory(...))`` to collect all points at once.

:func:`psi_trajectory` ``(geometry, h, k, l, psi_values, solution_key=NEAREST_ANGLES)``
    **Generator.** Yield motor-angle dicts for ψ (psi) rotation about Q.

:func:`trajectory_plan` ``(geometry, hkl_start, hkl_end, n_points, space="hkl", solution_key=NEAREST_ANGLES)``
    **Generator.** Yield motor-angle dicts between two reciprocal-space points,
    flagging limit violations and inaccessible regions.

:data:`NEAREST_ANGLES`
    Built-in solution-key callable: selects the candidate solution whose
    motor angles are closest (in least-squares sense) to the previous point.
    Suppresses branch-flip discontinuities across a trajectory.

Notes
-----
These functions *compute* trajectories; they do not move motors or
communicate with hardware.  A future execution layer (not yet planned)
will consume these outputs to drive real diffractometers.

The ordering instability of :meth:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward`
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
:meth:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward` returns multiple solutions in seed-discovery order, which is
not reproducible and can switch between the positive-chi and negative-chi
branches at adjacent trajectory points.  All three trajectory functions
accept a ``solution_key`` parameter — a callable that scores each candidate
solution against the previous point's chosen angles.  Passing
``NEAREST_ANGLES`` (the module default) prevents branch flips by always
preferring the candidate closest in motor space to the previous point.

ψ-scan algorithm
~~~~~~~~~~~~~~~~
Two distinct definitions of ψ exist in the literature.

``geometry.psi()`` implements You (1999) eqs. 10-11: ψ is the azimuthal
angle of the reference vector **n** about **Q** measured from the projection
of the fixed lab beam direction y_hat onto the phi-frame Q-perp plane.  For
a given (hkl, UB, n_hkl) this value is *constant* across all motor-angle
solutions that satisfy the Bragg condition — it is a crystal-orientation
diagnostic, not a motor-angle observable.

``psi_trajectory()`` implements the Busing & Levy (1967) *operational* ψ:
the angle through which the sample has been rotated about the scattering
vector Q relative to a chosen reference orientation.  This is the quantity
physically varied during a ψ scan on a real diffractometer.

Algorithm (Busing & Levy 1967, Section "Angle settings"):

1. Obtain Z₀ — the sample rotation matrix at the base forward() solution.
2. Compute Q̂_lab from the fixed detector position (ttheta unchanged).
3. For each ψ_target build Z(ψ) = R(Q̂_lab, −ψ_target) · Z₀.
   Left-multiplying by R(Q̂_lab, −ψ) rotates the sample about Q in the lab
   frame while the Bragg condition Z · q̂_phi = Q̂_lab is preserved.
4. Remove the contribution of any fixed outer stages and decompose the
   inner-three-stage part of Z(ψ) into motor angles analytically:

   * **Eulerian** (fourcv, fourch, psic-style, fivec-style):
     n_phi · Z · n_om = cos(χ) isolates χ; then ω and φ follow from
     projections of Z · n_om and Z^T · n_phi onto the plane perpendicular
     to their shared rotation axis.  Two χ branches are returned.

   * **Kappa** (kappa4cv, kappa4ch, kappa6c):
     n_kphi · Z · n_komega = cos(κ)·cos²(α) + sin²(α) isolates κ (where
     α = kappa_alpha_deg); komega and kphi follow by the same projection
     method.  Two κ branches (positive and negative) are returned.

5. All candidate solutions are filtered through stage limits; the
   ``solution_key`` (default: ``NEAREST_ANGLES``) selects the one closest
   in motor space to the previous ψ point, keeping motion continuous.

The BL1967 ψ is measured as the angle between the beam direction projected
perpendicular to Q in the phi frame at Z₀ versus at Z(ψ), rotating about
q̂_phi.  ψ = 0 by definition at the base forward() solution.

Trajectory dict specification for hkl_trajectory
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
The ``trajectory`` argument is a dict whose ``"type"`` key selects the
scan shape.  Supported types:

``"line"``
    Linear interpolation between two hkl points.
    Required keys: ``"start": (h, k, l)``, ``"end": (h, k, l)``.

``"radial"``
    Scan ±extent from center along a reciprocal-space direction.
    Required keys: ``"center": (h, k, l)``, ``"direction": (dh, dk, dl)``,
    ``"extent": delta``.  The direction is normalized; extent is in r.l.u.

``"transverse"``
    Scan ±extent from center perpendicular to a reference Q direction.
    Required keys: ``"center": (h, k, l)``, ``"Q_ref": (dh, dk, dl)``,
    ``"extent": delta``.  The transverse direction is the component of the
    first non-collinear standard basis vector perpendicular to ``Q_ref``
    (Gram-Schmidt).

References
----------
* Busing & Levy, Acta Cryst. 22, 457-464 (1967) —
  Four-circle geometry, ψ scan formulation, angle-setting equations.
* You, J. Appl. Cryst. 32, 614-623 (1999) —
  psic geometry, azimuthal angle definition (eqs. 10-11).
* ITC Vol. C, Sec. 2.2.6 (2006) —
  Kappa geometry, ψ scan via ω, χ, φ.
* D.A. Walko, Ref. Module Mater. Sci. Mater. Eng. (2016) —
  Kappa conversion formulae.
"""

from __future__ import annotations

import logging
import math
from collections.abc import Iterator
from typing import TYPE_CHECKING

import numpy as np

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Solution-key callables
# ---------------------------------------------------------------------------


def NEAREST_ANGLES(
    candidate: dict[str, float],
    previous: dict[str, float] | None,
) -> float:
    """
    Score a candidate motor-angle solution by its distance from the previous point.

    Returns the sum of squared angular differences between ``candidate`` and
    ``previous`` over all stages present in ``candidate``.  Lower score means
    closer to the previous point.

    When ``previous`` is ``None`` (first trajectory point) every candidate
    scores 0.0, so the first solution returned by ``forward()`` is used.

    Parameters
    ----------
    candidate : dict[str, float]
        Motor angles for one solution candidate.
    previous : dict[str, float] or None
        Motor angles chosen at the immediately preceding trajectory point.

    Returns
    -------
    float
        Non-negative score; lower is preferred.
    """
    if previous is None:
        return 0.0
    return sum((candidate.get(k, 0.0) - previous.get(k, 0.0)) ** 2 for k in candidate)


# ---------------------------------------------------------------------------
# Internal helpers — solution selection
# ---------------------------------------------------------------------------


def _pick_solution(
    solutions: list[dict[str, float]],
    previous: dict[str, float] | None,
    solution_key,
) -> dict[str, float] | None:
    """Return the best solution scored by ``solution_key``, or None if empty."""
    if not solutions:
        return None
    if solution_key is None:
        return solutions[0]
    return min(solutions, key=lambda c: solution_key(c, previous))


def _check_limits(
    geometry: AdHocDiffractometer,
    angles: dict[str, float],
) -> bool:
    """Return True if all supplied angles are within their stage limits."""
    for name, angle in angles.items():
        try:
            stage = geometry.stage(name)
        except KeyError:
            continue
        if not stage.in_limits(angle):
            return False
    return True


# ---------------------------------------------------------------------------
# Internal helpers — hkl point generation
# ---------------------------------------------------------------------------


def _hkl_points(trajectory: dict, n_points: int) -> list[tuple[float, float, float]]:
    """
    Generate ``n_points`` evenly-spaced (h, k, l) tuples from a trajectory dict.

    Raises
    ------
    ValueError
        If ``n_points < 2``, type is unknown, or a direction/Q_ref is zero.
    KeyError
        If a required trajectory key is missing.
    """
    if n_points < 2:
        raise ValueError(f"n_points must be at least 2; got {n_points}.")

    ttype = trajectory.get("type")

    if ttype == "line":
        start = np.asarray(trajectory["start"], dtype=float)
        end = np.asarray(trajectory["end"], dtype=float)
        return [
            tuple(float(v) for v in (start + t * (end - start)))
            for t in np.linspace(0.0, 1.0, n_points)
        ]

    if ttype == "radial":
        center = np.asarray(trajectory["center"], dtype=float)
        direction = np.asarray(trajectory["direction"], dtype=float)
        extent = float(trajectory["extent"])
        d_norm = np.linalg.norm(direction)
        if d_norm < 1e-14:
            raise ValueError(
                "hkl_trajectory: radial trajectory 'direction' must be non-zero."
            )
        d_hat = direction / d_norm
        return [
            tuple(float(v) for v in (center + t * d_hat))
            for t in np.linspace(-extent, extent, n_points)
        ]

    if ttype == "transverse":
        center = np.asarray(trajectory["center"], dtype=float)
        q_ref = np.asarray(trajectory["Q_ref"], dtype=float)
        extent = float(trajectory["extent"])
        q_norm = np.linalg.norm(q_ref)
        if q_norm < 1e-14:
            raise ValueError(
                "hkl_trajectory: transverse trajectory 'Q_ref' must be non-zero."
            )
        q_hat = q_ref / q_norm
        for c in (  # pragma: no branch
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ):
            perp = c - np.dot(c, q_hat) * q_hat
            if np.linalg.norm(perp) > 0.1:
                perp = perp / np.linalg.norm(perp)
                break
        return [
            tuple(float(v) for v in (center + t * perp))
            for t in np.linspace(-extent, extent, n_points)
        ]

    raise ValueError(
        f"hkl_trajectory: unknown trajectory type {ttype!r}. "
        "Supported types: 'line', 'radial', 'transverse'."
    )


# ---------------------------------------------------------------------------
# Internal helpers — Euler decomposition  (BL1967, generalized)
# ---------------------------------------------------------------------------


def _perp_ref(axis: np.ndarray) -> np.ndarray:
    """
    Return a unit vector perpendicular to ``axis``.

    Picks the first of [y, x, z] that is not nearly parallel to ``axis``,
    then Gram-Schmidt-orthogonalises it.  This avoids the conditional
    ``if abs(dot(ref, axis)) > 0.9`` scattered across the decomposition steps.
    """
    for candidate in (
        np.array([0.0, 1.0, 0.0]),
        np.array([1.0, 0.0, 0.0]),
        np.array([0.0, 0.0, 1.0]),
    ):
        v = candidate - np.dot(candidate, axis) * axis
        norm = np.linalg.norm(v)
        if norm > 0.1:
            return v / norm
    # Unreachable for any valid unit axis vector, but satisfies type checker.
    raise ValueError(
        f"Cannot find a reference perpendicular to {axis}"
    )  # pragma: no cover


def _euler_from_Z_standard(
    Z: np.ndarray,
    omega_stage: object,
    chi_stage: object,
    phi_stage: object,
) -> list[tuple[float, float, float]]:
    """
    Decompose Z = R(n_φ, φ) · R(n_χ, χ) · R(n_ω, ω) into (ω, χ, φ) in degrees.

    This is the innermost-first product order used by
    ``geometry.sample_rotation_matrix()``: the outermost stage (omega/eta)
    is applied first (rightmost), the innermost (phi) last (leftmost).

    The decomposition generalizes Busing & Levy (1967) eqs. 14-17 to
    arbitrary signed stage axes.  Two χ branches are returned
    (χ ∈ [0°, 180°] and χ ∈ [-180°, 0°]).

    Notes
    -----
    **Chi** is extracted from the element n_φ · Z · n_ω, which equals
    n_φ · R_χ · n_ω (since R_φ fixes n_φ and R_ω fixes n_ω).
    For the standard case n_ω ∥ n_φ this equals cos(χ).

    **Phi** is extracted from Z · n_ω = R_φ · (R_χ · n_ω): the angle in the
    plane ⊥ n_φ between (R_χ · n_ω) and (Z · n_ω).

    **Omega** is extracted from Z^T · n_φ = R_χ^T · n_φ: the angle in the
    plane ⊥ n_ω between (R_χ^T · n_φ) and (Z^T · n_φ), with sign reversed.

    Parameters
    ----------
    Z : ndarray, shape (3, 3)
        Sample rotation matrix to decompose.
    omega_stage, chi_stage, phi_stage : Stage
        The three innermost sample stages in stacking order (floor → sample).

    Returns
    -------
    list of (omega_deg, chi_deg, phi_deg)
        One or two solutions.
    """
    from .rotation import rotation_matrix as _Rmat

    n_om = np.asarray(omega_stage.axis, dtype=float)
    n_chi = np.asarray(chi_stage.axis, dtype=float)
    n_phi = np.asarray(phi_stage.axis, dtype=float)
    n_om_hat = n_om / np.linalg.norm(n_om)
    n_chi_hat = n_chi / np.linalg.norm(n_chi)
    n_phi_hat = n_phi / np.linalg.norm(n_phi)

    # --- Step 1: chi from n_phi · Z · n_om ----------------------------------
    elem = float(np.dot(n_phi_hat, Z @ n_om_hat))
    elem = max(-1.0, min(1.0, elem))
    chi_pos = math.degrees(math.acos(elem))
    chi_neg = -chi_pos

    results = []
    for chi_deg in (chi_pos, chi_neg):
        R_chi = _Rmat(n_chi_hat, chi_deg)

        # --- Step 2: phi from Z · n_om = R_phi · (R_chi · n_om) -------------
        Rchi_nom = R_chi @ n_om_hat
        Z_nom = Z @ n_om_hat
        e1 = _perp_ref(n_phi_hat)
        e2 = np.cross(n_phi_hat, e1)
        phi_rad = math.atan2(
            float(np.dot(Z_nom, e2)), float(np.dot(Z_nom, e1))
        ) - math.atan2(float(np.dot(Rchi_nom, e2)), float(np.dot(Rchi_nom, e1)))
        phi_rad = ((phi_rad + math.pi) % (2 * math.pi)) - math.pi
        phi_deg = math.degrees(phi_rad)

        # --- Step 3: omega from Z^T · n_phi = R_chi^T · n_phi ---------------
        Rchi_T_nphi = R_chi.T @ n_phi_hat
        ZT_nphi = Z.T @ n_phi_hat
        f1 = _perp_ref(n_om_hat)
        f2 = np.cross(n_om_hat, f1)
        omega_rad = math.atan2(
            float(np.dot(Rchi_T_nphi, f2)), float(np.dot(Rchi_T_nphi, f1))
        ) - math.atan2(float(np.dot(ZT_nphi, f2)), float(np.dot(ZT_nphi, f1)))
        omega_rad = ((omega_rad + math.pi) % (2 * math.pi)) - math.pi
        omega_deg = math.degrees(omega_rad)

        results.append((omega_deg, chi_deg, phi_deg))

    # De-duplicate degenerate case (chi = 0): both branches produce identical
    # chi — return only the first.
    if len(results) == 2 and abs(results[0][1] - results[1][1]) < 1e-8:
        return [results[0]]
    return results


def _kappa_from_Z(
    Z: np.ndarray,
    komega_stage: object,
    kappa_stage: object,
    kphi_stage: object,
    alpha_deg: float,
) -> list[tuple[float, float, float]]:
    r"""
    Decompose ``Z = R(n_kφ, kφ) · R(n_κ, κ) · R(n_kω, kω)`` into
    ``(kω, κ, kφ)`` using a Rodrigues expansion that is **agnostic to
    the specific signed axis convention** of the kappa preset.

    Derivation
    ----------
    Apply both sides to ``n_komega``:

        Z · n_kω  =  R(n_kφ, kφ) · R(n_κ, κ) · n_kω

    (since ``R(n_kω, ·) · n_kω = n_kω``).  Take the dot product with
    ``n_kphi``:

        n_kφ · Z · n_kω  =  n_kφ · R(n_κ, κ) · n_kω

    Expanding the right side via Rodrigues (with
    ``a = n_kphi · n_komega``,
    ``b = n_kphi · (n_kappa × n_komega)``,
    ``c = (n_kphi · n_kappa)·(n_kappa · n_komega)``):

        elem  =  a·cos(κ) + b·sin(κ) + c·(1 − cos(κ))
              =  (a − c)·cos(κ) + b·sin(κ) + c

    For the Walko (2016) textbook axis convention this reduces to
    ``cos(κ)·cos²(α) + sin²(α)``; for the geometry-aware kappa axis
    introduced by issue #241 (``n_kappa = cos(α)·n_komega +
    sin(α)·n_chi_eq``, with ``n_komega`` parallel to ``n_kphi``) the
    same expression evaluates to ``cos(κ)·sin²(α) + cos²(α)``.  The
    code below uses the general ``(a − c)·cos(κ) + b·sin(κ) = γ − c``
    form so that any future kappa convention with non-parallel
    komega/kphi axes is handled correctly.

    ``komega`` and ``kphi`` are then extracted by the same projection
    method as :func:`_euler_from_Z_standard`.

    Parameters
    ----------
    Z : ndarray, shape (3, 3)
    komega_stage, kappa_stage, kphi_stage : Stage
    alpha_deg : float
        Kappa tilt angle in degrees.  Retained for API compatibility;
        the decomposition itself uses only the geometric axes.

    Returns
    -------
    list of (komega_deg, kappa_deg, kphi_deg)
        Two solutions (positive and negative κ branches), or one if
        degenerate.

    References
    ----------
    * Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016) — textbook
      special case used by earlier versions of this function.
    * Issue #241 — geometry-aware re-derivation.
    """
    from .rotation import rotation_matrix as _Rmat

    n_km = np.asarray(komega_stage.axis, dtype=float) / np.linalg.norm(
        komega_stage.axis
    )
    n_kap = np.asarray(kappa_stage.axis, dtype=float) / np.linalg.norm(kappa_stage.axis)
    n_kph = np.asarray(kphi_stage.axis, dtype=float) / np.linalg.norm(kphi_stage.axis)

    # --- Step 1: κ from the Rodrigues expansion of n_kphi · Z · n_komega ----
    elem = float(np.clip(np.dot(n_kph, Z @ n_km), -1.0, 1.0))
    a_coeff = float(np.dot(n_kph, n_km))
    b_coeff = float(np.dot(n_kph, np.cross(n_kap, n_km)))
    c_coeff = float(np.dot(n_kph, n_kap)) * float(np.dot(n_kap, n_km))
    # Solve (a − c)·cos(κ) + b·sin(κ) = elem − c
    rhs = elem - c_coeff
    M = math.hypot(a_coeff - c_coeff, b_coeff)
    if M < 1e-14:  # pragma: no cover
        return []
    raw_ratio = rhs / M
    # The target Z is outside the kappa arm's reachable range when
    # |ratio| > 1.  Reject (no real κ exists) so the caller can mark
    # the trajectory point as inaccessible instead of returning a
    # spurious clamped solution.
    if raw_ratio > 1.0 + 1e-9 or raw_ratio < -1.0 - 1e-9:
        return []
    ratio = max(-1.0, min(1.0, raw_ratio))
    phi0 = math.atan2(b_coeff, a_coeff - c_coeff)
    delta = math.acos(ratio)
    kap_pos = math.degrees(phi0 + delta)
    kap_neg = math.degrees(phi0 - delta)
    # Normalise into (-180, 180]
    kap_pos = ((kap_pos + 180.0) % 360.0) - 180.0
    kap_neg = ((kap_neg + 180.0) % 360.0) - 180.0
    # Order so the smaller |κ| comes first (matches historic ordering).
    if abs(kap_pos) > abs(kap_neg):
        kap_pos, kap_neg = kap_neg, kap_pos

    results = []
    for kap_deg in (kap_pos, kap_neg):
        R_kap = _Rmat(n_kap, kap_deg)

        # --- Step 2: kphi from Z · n_km = R_kphi · (R_kap · n_km) ----------
        Rk_nom = R_kap @ n_km
        Z_nom = Z @ n_km
        e1 = _perp_ref(n_kph)
        e2 = np.cross(n_kph, e1)
        kphi_rad = math.atan2(
            float(np.dot(Z_nom, e2)), float(np.dot(Z_nom, e1))
        ) - math.atan2(float(np.dot(Rk_nom, e2)), float(np.dot(Rk_nom, e1)))
        kphi_rad = ((kphi_rad + math.pi) % (2 * math.pi)) - math.pi
        kphi_deg = math.degrees(kphi_rad)

        # --- Step 3: komega from Z^T · n_kphi = R_kap^T · n_kphi ------------
        Rk_T_nkph = R_kap.T @ n_kph
        ZT_nkph = Z.T @ n_kph
        f1 = _perp_ref(n_km)
        f2 = np.cross(n_km, f1)
        kom_rad = math.atan2(
            float(np.dot(Rk_T_nkph, f2)), float(np.dot(Rk_T_nkph, f1))
        ) - math.atan2(float(np.dot(ZT_nkph, f2)), float(np.dot(ZT_nkph, f1)))
        kom_rad = ((kom_rad + math.pi) % (2 * math.pi)) - math.pi
        kom_deg = math.degrees(kom_rad)

        results.append((kom_deg, kap_deg, kphi_deg))

    # De-duplicate (kappa = 0): both branches produce identical kappa values.
    if len(results) == 2 and abs(results[0][1] - results[1][1]) < 1e-8:
        return [results[0]]
    return results


# ---------------------------------------------------------------------------
# Internal helpers — psi-scan core
# ---------------------------------------------------------------------------


def _psi_candidates(
    geometry: AdHocDiffractometer,
    Z0: np.ndarray,
    Q_lab_hat: np.ndarray,
    y_eff: np.ndarray,
    psi_target: float,
    base: dict[str, float],
) -> list[dict[str, float]]:
    """
    Compute all in-limits motor-angle candidates for ``psi_target``.

    Implements BL1967 step 3-4: builds Z(ψ) = R(Q̂_lab, −ψ) · Z₀ then
    decomposes the inner three stages analytically.

    Parameters
    ----------
    geometry : AdHocDiffractometer
    Z0 : ndarray (3, 3)
        Sample rotation matrix at the base forward() solution.
    Q_lab_hat : ndarray (3,)
        Unit scattering vector in the lab frame (from detector position).
    y_eff : ndarray (3,)
        Effective incident-beam direction (accounts for inclination).
    psi_target : float
        Target ψ in degrees.
    base : dict[str, float]
        Base motor angles (Bragg solution); used for outer stages and ttheta.

    Returns
    -------
    list of dict[str, float]
        All limit-passing candidate motor-angle dicts.
    """
    from .rotation import rotation_matrix as _Rmat

    # Identify stage groups
    sample_stages = geometry.sample_stages
    if len(sample_stages) < 3:
        return []

    outer_stages = sample_stages[:-3]
    omega_s, chi_s, phi_s = sample_stages[-3], sample_stages[-2], sample_stages[-1]

    # Build R_outer from the base angles of any fixed outer stages
    R_outer = np.eye(3)
    for s in outer_stages:
        R_outer = _Rmat(s.axis, base.get(s.name, s.angle)) @ R_outer

    # Rotate Z0 about Q_lab by -psi_target (negative because the convention
    # is that increasing psi rotates the sample such that the beam-perpendicular
    # direction in the phi frame rotates by +psi relative to Z0)
    Z_new = _Rmat(Q_lab_hat, -psi_target) @ Z0
    Z_inner = R_outer.T @ Z_new

    # Decompose the inner three stages
    is_kappa = geometry.kappa_alpha_deg is not None
    if is_kappa:
        sols = _kappa_from_Z(Z_inner, omega_s, chi_s, phi_s, geometry.kappa_alpha_deg)
    else:
        sols = _euler_from_Z_standard(Z_inner, omega_s, chi_s, phi_s)

    candidates = []
    for a0, a1, a2 in sols:
        angles = dict(base)
        angles[omega_s.name] = a0
        angles[chi_s.name] = a1
        angles[phi_s.name] = a2
        if _check_limits(geometry, angles):
            candidates.append(angles)

    return candidates


def _measure_psi(
    geometry: AdHocDiffractometer,
    angles: dict[str, float],
    Q_lab_hat: np.ndarray,
    y_eff: np.ndarray,
    q_hat_phi: np.ndarray,
    ref_dir: np.ndarray,
) -> float:
    """
    Measure the BL1967 operational ψ for ``angles`` relative to ``ref_dir``.

    ψ is the angle from ``ref_dir`` to (Z^T · y_eff)_⊥ about q̂_phi,
    where the perpendicular projection removes the Q component.
    """
    saved: dict[str, float] = {}
    try:
        for name, angle in angles.items():
            saved[name] = geometry.stage(name).angle
            geometry.set_angle(name, float(angle))
        Z_act = geometry.sample_rotation_matrix()
    finally:
        for name, angle in saved.items():
            geometry.set_angle(name, angle)

    y_phi = Z_act.T @ y_eff
    y_perp = y_phi - np.dot(y_phi, q_hat_phi) * q_hat_phi
    new_dir = y_perp / np.linalg.norm(y_perp)
    cos_a = float(np.clip(np.dot(ref_dir, new_dir), -1.0, 1.0))
    sin_a = float(np.dot(q_hat_phi, np.cross(ref_dir, new_dir)))
    return math.degrees(math.atan2(sin_a, cos_a))


# ---------------------------------------------------------------------------
# Public functions
# ---------------------------------------------------------------------------


def hkl_trajectory(
    geometry: AdHocDiffractometer,
    trajectory: dict,
    n_points: int,
    *,
    solution_key=NEAREST_ANGLES,
) -> Iterator[dict]:
    """
    Yield motor-angle dicts along a reciprocal-space trajectory.

    For each of the ``n_points`` equally-spaced hkl positions along the
    trajectory, calls ``geometry.forward()`` to find all valid motor-angle
    solutions, then uses ``solution_key`` to select one.  Points with no
    valid solution are yielded with ``angles=None``.

    This function is a **generator**: it yields one dict per trajectory point
    and computes each point on demand.  To get all points at once use
    ``list(hkl_trajectory(...))``.

    The ``NEAREST_ANGLES`` solution key (the default) prevents branch-flip
    discontinuities: at each point it selects the candidate whose motor
    angles are closest to those of the previous point.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``wavelength`` and ``sample.UB`` set and an active mode.
    trajectory : dict
        Trajectory specification.  The ``"type"`` key selects the scan
        shape; see module docstring for supported types and required keys.
    n_points : int
        Number of evenly-spaced points to compute (≥ 2).
    solution_key : callable or None, optional
        Scoring function ``(candidate, previous) -> float``; lower wins.
        Default: ``NEAREST_ANGLES``.  Pass ``None`` to use ``forward()``'s
        raw ordering.

    Yields
    ------
    dict
        One entry per trajectory point with keys:

        ``"hkl"`` : tuple of float
            The requested Miller indices (h, k, l).
        ``"angles"`` : dict[str, float] or None
            Motor angles in degrees, or ``None`` if no valid solution.
        ``"warning"`` : str or None
            Descriptive warning when ``angles`` is ``None``, else ``None``.

    Raises
    ------
    ValueError
        If ``n_points < 2``, trajectory type is unknown, or a required
        trajectory key is missing.
    ValueError
        If ``geometry.wavelength`` is None or ``geometry.sample.UB`` is None.
    NotImplementedError
        If no active diffraction mode is set.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.wavelength = 1.5406
    >>> g.sample.lattice = ahd.Lattice(a=4.0)
    >>> ahd.ub_identity(g.sample)
    >>> g.mode_name = "bisecting"
    >>> for pt in ahd.hkl_trajectory(
    ...     g,
    ...     {"type": "line", "start": (1, 0, 0), "end": (3, 0, 0)},
    ...     n_points=5,
    ... ):
    ...     print(pt["hkl"], pt["angles"])
    """
    points = _hkl_points(trajectory, n_points)
    previous: dict[str, float] | None = None

    for hkl in points:
        h, k, l = hkl  # noqa: E741
        try:
            solutions = geometry.forward(h, k, l)
        except (ValueError, NotImplementedError) as exc:
            yield {"hkl": hkl, "angles": None, "warning": str(exc)}
            continue

        if not solutions:
            yield {
                "hkl": hkl,
                "angles": None,
                "warning": "no valid motor solution (all candidates outside limits)",
            }
            continue

        chosen = _pick_solution(solutions, previous, solution_key)
        yield {"hkl": hkl, "angles": chosen, "warning": None}
        previous = chosen


def psi_trajectory(
    geometry: AdHocDiffractometer,
    h: float,
    k: float,
    l: float,  # noqa: E741
    psi_values,
    *,
    solution_key=NEAREST_ANGLES,
) -> Iterator[dict]:
    """
    Yield motor-angle dicts for ψ rotation about the scattering vector Q.

    For the fixed reflection (h, k, l) and each target ψ in ``psi_values``,
    finds motor angles that simultaneously satisfy the Bragg condition and
    produce the requested BL1967 operational ψ.

    ψ Definition (Busing & Levy 1967)
    ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ψ is the angle through which the sample has rotated about Q relative to
    the base ``forward()`` solution.  ψ = 0 at the base solution by definition.
    Physically, ψ is varied by rotating the sample about Q (changing ω, χ, φ
    or their kappa equivalents) while keeping 2θ fixed.

    This differs from ``geometry.psi()``, which implements the You (1999)
    definition and is constant for all Bragg solutions of a given hkl.  See
    the module docstring for the full discussion.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``wavelength`` and ``sample.UB`` set, and an active
        diffraction mode.  ``kappa_alpha_deg`` must be set for kappa
        geometries (it is set automatically by the factory functions).
    h, k, l : float
        Miller indices of the fixed reflection.
    psi_values : iterable of float
        Target ψ angles in degrees relative to the base forward() solution.
    solution_key : callable or None, optional
        Scoring function ``(candidate, previous) -> float``; lower wins.
        Default: ``NEAREST_ANGLES``.

    Yields
    ------
    dict
        One entry per ψ value with keys:

        ``"psi_target"`` : float
            The requested ψ (degrees).
        ``"psi_actual"`` : float or None
            Achieved ψ verified by internal measurement, or ``None``.
        ``"angles"`` : dict[str, float] or None
            Motor angles in degrees, or ``None``.
        ``"warning"`` : str or None
            Descriptive warning when ``angles`` is ``None``, else ``None``.

    Raises
    ------
    ValueError
        If ``geometry.wavelength`` is None or ``geometry.sample.UB`` is None.
    ValueError
        If ``(h, k, l) == (0, 0, 0)``.
    NotImplementedError
        If no active diffraction mode is set.

    Notes
    -----
    The geometry must have at least three sample stages (omega/eta, chi/kappa,
    phi/kphi).  Geometries with fewer sample stages are not supported.

    References
    ----------
    Busing & Levy, Acta Cryst. 22, 457-464 (1967), angle-setting equations.
    ITC Vol. C, Sec. 2.2.6 (2006), kappa ψ scan.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.wavelength = 1.5406
    >>> g.sample.lattice = ahd.Lattice(a=4.0)
    >>> ahd.ub_identity(g.sample)
    >>> g.mode_name = "bisecting"
    >>> for pt in ahd.psi_trajectory(g, 1, 1, 0, range(-90, 91, 30)):
    ...     print(pt["psi_target"], pt["psi_actual"])
    """

    # Obtain the Bragg base solutions (raises ValueError/NotImplementedError
    # for bad preconditions — wavelength, UB, hkl=0, mode).
    base_solutions = geometry.forward(h, k, l)

    if not base_solutions:
        no_bragg_msg = (
            "no valid Bragg solution for this reflection "
            "(all candidates outside limits)"
        )
        for p in psi_values:
            yield {
                "psi_target": float(p),
                "psi_actual": None,
                "angles": None,
                "warning": no_bragg_msg,
            }
        return

    # Pick the base solution (used to compute Z0 and as the psi=0 reference)
    base = base_solutions[0]

    # Compute Z0 and Q_lab_hat from base (temporarily apply motor angles)
    saved: dict[str, float] = {}
    try:
        for name, angle in base.items():
            saved[name] = geometry.stage(name).angle
            geometry.set_angle(name, float(angle))
        Z0 = geometry.sample_rotation_matrix()
        D = geometry.detector_rotation_matrix()
    finally:
        for name, angle in saved.items():
            geometry.set_angle(name, angle)

    y_lab = np.asarray(geometry.basis["longitudinal"], dtype=float)
    y_lab = y_lab / np.linalg.norm(y_lab)
    y_eff = geometry.inclination_matrix.T @ y_lab
    Q_lab_vec = (2.0 * math.pi / geometry.wavelength) * (D @ y_eff - y_eff)
    Q_lab_hat = Q_lab_vec / np.linalg.norm(Q_lab_vec)

    q_hat_phi = geometry.sample.UB @ np.array([float(h), float(k), float(l)])
    q_hat_phi = q_hat_phi / np.linalg.norm(q_hat_phi)

    # Build reference direction: y_eff projected perp to Q in phi frame at Z0
    y_phi_0 = Z0.T @ y_eff
    y_perp_0 = y_phi_0 - np.dot(y_phi_0, q_hat_phi) * q_hat_phi
    ref_dir = y_perp_0 / np.linalg.norm(y_perp_0)

    previous: dict[str, float] | None = None

    for psi_target in psi_values:
        psi_target = float(psi_target)

        candidates = _psi_candidates(geometry, Z0, Q_lab_hat, y_eff, psi_target, base)

        if not candidates:
            yield {
                "psi_target": psi_target,
                "psi_actual": None,
                "angles": None,
                "warning": (
                    f"no motor solution achieves psi = {psi_target:.4g}° "
                    "within stage limits"
                ),
            }
            continue

        chosen = _pick_solution(candidates, previous, solution_key)
        psi_actual = _measure_psi(
            geometry, chosen, Q_lab_hat, y_eff, q_hat_phi, ref_dir
        )

        yield {
            "psi_target": psi_target,
            "psi_actual": psi_actual,
            "angles": chosen,
            "warning": None,
        }
        previous = chosen


def trajectory_plan(
    geometry: AdHocDiffractometer,
    hkl_start: tuple[float, float, float],
    hkl_end: tuple[float, float, float],
    n_points: int,
    *,
    space: str = "hkl",
    solution_key=NEAREST_ANGLES,
) -> Iterator[dict]:
    r"""
    Yield motor-angle dicts for a planned path between two reciprocal-space points.

    Interpolates between ``hkl_start`` and ``hkl_end`` in ``n_points``
    equally-spaced steps, calls ``geometry.forward()`` at each point, and
    flags portions that are inaccessible (no valid motor solution) or that
    exceed stage limits.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Must have ``wavelength`` and ``sample.UB`` set and an active mode.
    hkl_start, hkl_end : tuple of float
        Start and end Miller indices (h, k, l).
    n_points : int
        Number of evenly-spaced points including endpoints (≥ 2).
    space : {"hkl", "Q"}, optional
        Interpolation space.

        ``"hkl"`` (default)
            Linearly interpolate Miller indices.  Natural for L-scans,
            H-scans, and crystal-axis-aligned trajectories.  Equal steps in
            hkl are *not* equal steps in \|Q\| for non-cubic lattices.

        ``"Q"``
            Linearly interpolate in reciprocal Cartesian coordinates
            (Qx, Qy, Qz = UB @ hkl).  Equal steps in Q regardless of
            lattice distortion; intermediate hkl values are generally
            non-integer.  Endpoints are exact by construction.

    solution_key : callable or None, optional
        Scoring function ``(candidate, previous) -> float``; lower wins.
        Default: ``NEAREST_ANGLES``.

    Yields
    ------
    dict
        One entry per trajectory point with keys:

        ``"hkl"`` : tuple of float
            The (h, k, l) point.
        ``"angles"`` : dict[str, float] or None
            Motor angles, or ``None`` if inaccessible.
        ``"accessible"`` : bool
            ``True`` if a valid in-limits solution was found.
        ``"warnings"`` : list of str
            Empty when accessible; otherwise one or more warning strings.

    Raises
    ------
    ValueError
        If ``n_points < 2``.
    ValueError
        If ``space`` is not ``"hkl"`` or ``"Q"``.
    ValueError
        If ``space="Q"`` and ``geometry.sample.UB`` is None or singular.
    ValueError
        If ``geometry.wavelength`` is None or ``geometry.sample.UB`` is None.
    NotImplementedError
        If no active diffraction mode is set.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.wavelength = 1.5406
    >>> g.sample.lattice = ahd.Lattice(a=4.0)
    >>> ahd.ub_identity(g.sample)
    >>> g.mode_name = "bisecting"
    >>> plan = list(ahd.trajectory_plan(g, (1, 0, 0), (3, 0, 0), n_points=11))
    >>> accessible = [pt for pt in plan if pt["accessible"]]
    >>> print(f"{len(accessible)}/{len(plan)} points accessible")
    """
    if n_points < 2:
        raise ValueError(f"n_points must be at least 2; got {n_points}.")
    if space not in ("hkl", "Q"):
        raise ValueError(f"trajectory_plan: space must be 'hkl' or 'Q'; got {space!r}.")

    hkl_start_arr = np.asarray(hkl_start, dtype=float)
    hkl_end_arr = np.asarray(hkl_end, dtype=float)

    if space == "hkl":
        hkl_points = [
            tuple(float(v) for v in (hkl_start_arr + t * (hkl_end_arr - hkl_start_arr)))
            for t in np.linspace(0.0, 1.0, n_points)
        ]
    else:  # space == "Q"
        sample = geometry.sample
        if sample.UB is None:
            raise ValueError(
                "trajectory_plan(space='Q') requires geometry.sample.UB to be set."
            )
        try:
            UB_inv = np.linalg.inv(sample.UB)
        except np.linalg.LinAlgError as exc:
            raise ValueError(
                "trajectory_plan(space='Q'): UB matrix is singular."
            ) from exc

        Q_start = sample.UB @ hkl_start_arr
        Q_end = sample.UB @ hkl_end_arr
        hkl_points = []
        for t in np.linspace(0.0, 1.0, n_points):
            Q_interp = Q_start + t * (Q_end - Q_start)
            hkl_points.append(tuple(float(v) for v in (UB_inv @ Q_interp)))
        # Force exact endpoints
        hkl_points[0] = tuple(float(v) for v in hkl_start_arr)
        hkl_points[-1] = tuple(float(v) for v in hkl_end_arr)

    previous: dict[str, float] | None = None

    for hkl in hkl_points:
        h, k, l = hkl  # noqa: E741
        warnings_list: list[str] = []

        try:
            solutions = geometry.forward(h, k, l)
        except (ValueError, NotImplementedError) as exc:
            yield {
                "hkl": hkl,
                "angles": None,
                "accessible": False,
                "warnings": [str(exc)],
            }
            continue

        if not solutions:
            yield {
                "hkl": hkl,
                "angles": None,
                "accessible": False,
                "warnings": ["no valid motor solution (all candidates outside limits)"],
            }
            continue

        chosen = _pick_solution(solutions, previous, solution_key)

        try:
            geometry.check_limits(**chosen)
        except ValueError as exc:
            warnings_list.append(str(exc))

        yield {
            "hkl": hkl,
            "angles": chosen,
            "accessible": len(warnings_list) == 0,
            "warnings": warnings_list,
        }
        if not warnings_list:
            previous = chosen
