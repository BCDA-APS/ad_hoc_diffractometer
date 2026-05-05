# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
kappa.py — Kappa-to-Eulerian angle conversion.

Kappa diffractometers replace the Eulerian chi circle with a kappa arm
tilted at angle ``alpha_0`` from the omega axis.  This module provides
**two layers** of conversion between the real motor triple
``(komega, kappa, kphi)`` and the equivalent virtual Eulerian
pseudoangles ``(omega, chi, phi)``:

1. **Geometry-aware (recommended)**.  The functions
   :func:`eulerian_to_kappa_axes` and :func:`kappa_to_eulerian_axes`
   take a :class:`KappaPseudoAngleConvention` that names the *actual*
   stage axis vectors of the kappa demo geometry (``n_komega``, ``n_kappa``,
   ``n_kphi``) and the *equivalent Eulerian* chi axis (``n_chi_eq``).
   The decomposition is derived directly from these axes and is correct
   for any signed-axis combination — including the ``ad_hoc_diffractometer``
   demo geometries that use mixed handedness or a horizontal scattering plane.

   These are the conversions used internally by the kappa virtual-angle
   solver.

2. **Walko (2016) eq. [16] textbook formula**.  The functions
   :func:`kappa_to_eulerian` and :func:`eulerian_to_kappa` implement
   Walko's published closed-form relations literally:

       chi   = 2 arcsin[sin(kappa/2) · sin(alpha_0)]
       offset = arccos[cos(kappa/2) / cos(chi/2)]
       omega = komega − offset
       phi   = kphi   − offset

   These are correct **only for the axis convention assumed in
   Walko's derivation** (omega about the transverse axis, chi about
   the longitudinal axis, phi about the transverse axis, all
   right-handed in the textbook sense).  They are retained as
   reference implementations of the published formula and for
   backward compatibility with users who study Walko's algebra
   directly.  They are **not** used inside the kappa virtual-angle
   solver, because no demo geometry in this package matches the textbook
   convention exactly: ``kappa4cv`` uses BL ``-TRANSVERSE`` for omega
   (left-handed); ``kappa4ch`` uses ``-VERTICAL`` (horizontal
   scattering plane); ``kappa6c`` uses You ``-TRANSVERSE`` (mixed
   handedness with a horizontal mu base).  See issue #241 for the
   diagnostic that uncovered this divergence.

Both directions are implemented as pure-NumPy / standard-library code.

Functions
---------
:func:`kappa_to_eulerian_axes`
    Geometry-aware forward direction: real (komega, kappa, kphi) →
    virtual (omega, chi, phi) for arbitrary signed stage axes.

:func:`eulerian_to_kappa_axes`
    Geometry-aware inverse direction: virtual (omega, chi, phi) →
    real (komega, kappa, kphi) for arbitrary signed stage axes.

:func:`kappa_to_eulerian`
    Walko (2016) eq. [16] forward direction (textbook frame only).

:func:`eulerian_to_kappa`
    Walko (2016) eq. [16] inverse direction (textbook frame only).

:func:`is_kappa_virtual_mode`
    Detect whether a mode contains a virtual Eulerian-angle constraint.

:func:`solve_kappa_virtual`
    Solve the virtual-angle forward problem for a kappa geometry.

Classes
-------
:class:`KappaPseudoAngleConvention`
    Container holding the four signed unit axes that fully specify
    the geometry-aware pseudoangle decomposition for a given kappa
    demo geometry.

Notes
-----
The default kappa tilt angle is 50° (Walko 2016; Enraf-Nonius
convention; ITC Vol. C §2.2.6).

The geometry-aware decomposition has two solutions distinguished by
the sign of kappa (the *branch* parameter).  The positive branch
(``kappa ≥ 0``) is the default; the negative branch gives
``kappa ≤ 0``.  The caller should choose the branch that keeps all
angles within the hardware limits.

The decomposition fails (``ValueError``) when the requested Eulerian
orientation cannot be represented by any kappa motor triple — this
is the kappa-arm reachability limit and is geometry-dependent (in
the textbook convention it reduces to ``|chi| < 2·alpha_0``).

References
----------
* Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), eq. [16] —
  closed-form pseudoangle relations in the textbook frame.
* ITC Vol. C §2.2.6 (2006).
* Issue #241 — diagnostic and resolution that motivated the
  geometry-aware layer.
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import TYPE_CHECKING

import numpy as np

from .rotation import _rotation_matrix_normalized

if TYPE_CHECKING:  # pragma: no cover
    from .diffractometer import AdHocDiffractometer

#: Virtual Eulerian pseudoangle names on kappa geometries.
KAPPA_VIRTUAL_ANGLES: frozenset[str] = frozenset({"omega", "chi", "phi"})


# ---------------------------------------------------------------------------
# Geometry-aware pseudoangle convention
# ---------------------------------------------------------------------------


@dataclass(frozen=True)
class KappaPseudoAngleConvention:
    """
    Per-geometry declaration of the four signed axes that determine
    the kappa ↔ Eulerian pseudoangle decomposition.

    A kappa geometry encodes its three sample stages
    (``komega``, ``kappa``, ``kphi``) with specific signed axis vectors
    in the lab frame.  An *equivalent Eulerian triple*
    (``omega``, ``chi``, ``phi``) shares the outer (``omega``) and
    inner (``phi``) axes with the kappa stack — by construction — and
    introduces a virtual ``chi`` axis perpendicular to ``omega`` that
    plays the role of the Eulerian chi rotation.

    The four axes here are sufficient to derive the conversion in
    closed form from the rotation-matrix identity

        R(n_komega, κω) · R(n_kappa, κ) · R(n_kphi, κφ)
            =  R(n_komega, ω) · R(n_chi_eq, χ) · R(n_kphi, φ).

    This identity holds for **arbitrary signed stage axes** as long
    as ``n_kappa`` lies in the plane spanned by ``n_komega`` and
    ``n_chi_eq`` (which is the geometric definition of a kappa arm).

    Parameters
    ----------
    n_komega : numpy.ndarray, shape (3,)
        Unit axis of the outer kappa stage (= virtual ω axis).
    n_kappa : numpy.ndarray, shape (3,)
        Unit axis of the kappa arm.
    n_kphi : numpy.ndarray, shape (3,)
        Unit axis of the inner kappa stage (= virtual φ axis).
    n_chi_eq : numpy.ndarray, shape (3,)
        Unit axis of the equivalent Eulerian χ rotation, perpendicular
        to ``n_komega`` and coplanar with ``n_komega`` and ``n_kappa``.

    Notes
    -----
    The convention object is immutable (``frozen=True``) and arrays
    are stored as a copy so the caller cannot accidentally mutate
    them.  Use :func:`make_kappa_pseudo_angle_convention` to construct
    one from physical-direction names + a basis dict.

    References
    ----------
    * Issue #241 — derivation and resolution.
    """

    n_komega: np.ndarray
    n_kappa: np.ndarray
    n_kphi: np.ndarray
    n_chi_eq: np.ndarray

    def __post_init__(self) -> None:
        # Normalise and freeze every axis.  ``object.__setattr__`` is
        # required because the dataclass is frozen.
        for name in ("n_komega", "n_kappa", "n_kphi", "n_chi_eq"):
            v = np.asarray(getattr(self, name), dtype=float).copy()
            n = float(np.linalg.norm(v))
            if n < 1e-14:
                raise ValueError(f"KappaPseudoAngleConvention: {name} has zero norm.")
            v = v / n
            v.setflags(write=False)
            object.__setattr__(self, name, v)


def make_kappa_pseudo_angle_convention(
    n_komega: np.ndarray,
    n_kappa: np.ndarray,
    n_kphi: np.ndarray,
    n_chi_eq: np.ndarray,
) -> KappaPseudoAngleConvention:
    """
    Construct a :class:`KappaPseudoAngleConvention` from four axes.

    Thin wrapper that exists so callers do not need to import the
    dataclass directly; the validation lives in
    :meth:`KappaPseudoAngleConvention.__post_init__`.

    Parameters
    ----------
    n_komega, n_kappa, n_kphi, n_chi_eq : array-like, shape (3,)
        Signed axis vectors.  Will be normalized internally.

    Returns
    -------
    KappaPseudoAngleConvention
    """
    return KappaPseudoAngleConvention(
        n_komega=n_komega,
        n_kappa=n_kappa,
        n_kphi=n_kphi,
        n_chi_eq=n_chi_eq,
    )


def kappa_axis_from_eulerian(
    n_komega: np.ndarray,
    n_chi_eq: np.ndarray,
    alpha_deg: float,
) -> np.ndarray:
    r"""
    Compute the kappa rotation axis as a tilt of the omega axis toward
    the equivalent Eulerian chi axis.

    The canonical definition of a kappa diffractometer (Walko 2016;
    ITC Vol. C §2.2.6) is that the kappa axis lies in the plane
    spanned by the outer ``omega`` axis and the equivalent chi axis,
    tilted by ``α`` from ``n_komega`` toward ``n_chi_eq``:

        n_kappa = n_komega · cos(α) + n_chi_eq · sin(α)

    This formula is **geometry-aware**: the resulting axis is correct
    for any signed convention of ``n_komega`` and ``n_chi_eq``,
    including the BL ``-TRANSVERSE`` orientation used by the
    ``kappa4cv`` demo geometry and the BL ``-VERTICAL`` orientation used by
    the ``kappa4ch`` demo geometry.

    The historic helper :func:`~ad_hoc_diffractometer.axes.kappa_axis`
    instead returns ``vertical·cos(α) + transverse·sin(α)``
    unconditionally, which assumes ``n_komega = +vertical`` and is
    therefore wrong for any demo geometry whose ``n_komega`` is not literally
    ``+vertical``.  See issue #241 for the diagnostic that uncovered
    this.

    Parameters
    ----------
    n_komega : array-like, shape (3,)
        Outer kappa stage axis.  Need not be normalized.
    n_chi_eq : array-like, shape (3,)
        Equivalent Eulerian chi axis (perpendicular to ``n_komega``).
        Need not be normalized.
    alpha_deg : float
        Kappa tilt angle in degrees, measured from ``n_komega``
        toward ``n_chi_eq``.  Typical value 50° (Enraf-Nonius
        convention).

    Returns
    -------
    n_kappa : numpy.ndarray, shape (3,)
        Unit vector of the kappa axis in the lab frame.

    Raises
    ------
    ValueError
        If ``n_komega`` and ``n_chi_eq`` are not perpendicular within
        ``1e-9``.

    References
    ----------
    * D. A. Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016) — kappa
      arm geometry.
    * ITC Vol. C, Sec. 2.2.6 (2006) — Enraf-Nonius α convention.
    * Issue #241 — divergence between the canonical α definition and
      the historic vertical/transverse hard-coded helper.
    """
    n_om = np.asarray(n_komega, dtype=float)
    n_om = n_om / np.linalg.norm(n_om)
    n_ch = np.asarray(n_chi_eq, dtype=float)
    n_ch = n_ch / np.linalg.norm(n_ch)
    if abs(float(np.dot(n_om, n_ch))) > 1e-9:
        raise ValueError(
            "kappa_axis_from_eulerian: n_komega and n_chi_eq must be "
            "perpendicular; got dot product "
            f"{float(np.dot(n_om, n_ch)):.6e}."
        )
    a = math.radians(alpha_deg)
    return math.cos(a) * n_om + math.sin(a) * n_ch


# ---------------------------------------------------------------------------
# Geometry-aware kappa ↔ Eulerian decomposition (issue #241)
# ---------------------------------------------------------------------------


def _eulerian_rotation(
    convention: KappaPseudoAngleConvention,
    omega_deg: float,
    chi_deg: float,
    phi_deg: float,
) -> np.ndarray:
    """Compose the equivalent-Eulerian sample rotation matrix.

    Stage order matches :func:`~orientation._compute_q_phi`: the
    innermost stage (``phi``) is the leftmost matrix in the product,
    so that ``Z · v`` rotates ``v`` by the outermost stage first when
    applied right-to-left.
    """
    R_om = _rotation_matrix_normalized(convention.n_komega, omega_deg)
    R_ch = _rotation_matrix_normalized(convention.n_chi_eq, chi_deg)
    R_ph = _rotation_matrix_normalized(convention.n_kphi, phi_deg)
    return R_ph @ R_ch @ R_om


def _kappa_rotation(
    convention: KappaPseudoAngleConvention,
    komega_deg: float,
    kappa_deg: float,
    kphi_deg: float,
) -> np.ndarray:
    """Compose the real kappa-stack sample rotation matrix.

    Stage order matches :func:`~orientation._compute_q_phi`: innermost
    stage (``kphi``) is the leftmost matrix in the product.
    """
    R_om = _rotation_matrix_normalized(convention.n_komega, komega_deg)
    R_ka = _rotation_matrix_normalized(convention.n_kappa, kappa_deg)
    R_ph = _rotation_matrix_normalized(convention.n_kphi, kphi_deg)
    return R_ph @ R_ka @ R_om


def _angle_about_axis(v_from: np.ndarray, v_to: np.ndarray, axis: np.ndarray) -> float:
    """
    Return the signed rotation angle (deg) about ``axis`` that takes
    ``v_from`` to ``v_to``.  Both vectors are assumed to lie on the
    same cone about ``axis`` (i.e. ``v_from · axis == v_to · axis``).
    """
    # Project onto plane perpendicular to axis
    p_from = v_from - np.dot(v_from, axis) * axis
    p_to = v_to - np.dot(v_to, axis) * axis
    n_from = float(np.linalg.norm(p_from))
    n_to = float(np.linalg.norm(p_to))
    if n_from < 1e-12 or n_to < 1e-12:
        # Degenerate: the vectors are parallel to the axis.  Any
        # rotation about the axis maps one to the other.  Return 0.
        return 0.0
    p_from = p_from / n_from
    p_to = p_to / n_to
    cos_t = float(np.clip(np.dot(p_from, p_to), -1.0, 1.0))
    sin_t = float(np.dot(axis, np.cross(p_from, p_to)))
    return math.degrees(math.atan2(sin_t, cos_t))


def _solve_alpha_cos_beta_sin_eq_gamma(
    alpha: float, beta: float, gamma: float
) -> list[float]:
    """
    Solve ``alpha · cos(θ) + beta · sin(θ) = gamma`` for θ in degrees.

    Returns the two solutions (or one solution duplicated, when the
    discriminant is zero) in the range (-180, 180].  Raises
    :class:`ValueError` when the equation has no real solution
    (``|gamma| > sqrt(alpha² + beta²)``).
    """
    M = math.hypot(alpha, beta)
    if M < 1e-14:
        # Degenerate: equation is 0 = gamma.  Any θ solves if gamma==0.
        if abs(gamma) < 1e-12:  # pragma: no cover
            return [0.0, 180.0]
        raise ValueError(
            "alpha·cos + beta·sin = gamma is unsolvable: "
            f"alpha={alpha}, beta={beta}, gamma={gamma} (LHS amplitude is zero)."
        )
    ratio = gamma / M
    # Allow tiny excursions from the unit interval that arise from
    # floating-point error.
    if ratio > 1.0 + 1e-9 or ratio < -1.0 - 1e-9:
        raise ValueError(
            "alpha·cos(θ) + beta·sin(θ) = gamma has no real solution: "
            f"|gamma|/sqrt(alpha²+beta²) = {abs(ratio):.6f} > 1."
        )
    ratio = max(-1.0, min(1.0, ratio))
    phi0 = math.atan2(beta, alpha)
    delta = math.acos(ratio)
    return [
        math.degrees(phi0 + delta),
        math.degrees(phi0 - delta),
    ]


def eulerian_to_kappa_axes(
    omega_deg: float,
    chi_deg: float,
    phi_deg: float,
    convention: KappaPseudoAngleConvention,
    branch: int = +1,
) -> tuple[float, float, float]:
    """
    Convert virtual Eulerian pseudoangles to real kappa motor angles
    using the geometry's actual signed stage axes.

    Inverts the rotation-matrix identity (with stage order matching
    :func:`~orientation._compute_q_phi` — innermost on the left)

        R(n_kphi, κφ) · R(n_kappa, κ) · R(n_komega, κω)
            =  R(n_kphi, φ) · R(n_chi_eq, χ) · R(n_komega, ω).

    The decomposition is fully analytic and selects between the two
    kappa-branch solutions via the ``branch`` parameter.

    Algorithm
    ---------
    Let ``R_eul = R(n_kphi, φ) · R(n_chi_eq, χ) · R(n_komega, ω)``
    be the target rotation.  Apply both sides to ``n_komega``:

        R(n_kphi, κφ) · R(n_kappa, κ) · n_komega  =  R_eul · n_komega

    (since ``R(n_komega, ·) · n_komega = n_komega``).  Take the dot
    product with ``n_kphi``:

        n_kphi · R(n_kappa, κ) · n_komega  =  n_kphi · R_eul · n_komega.

    Expanding ``R(n_kappa, κ) · n_komega`` via the Rodrigues formula
    yields a single trigonometric equation in κ of the form
    ``alpha·cos(κ) + beta·sin(κ) = gamma`` with two solutions
    (the two kappa branches).

    Once κ is fixed, ``κφ`` is the rotation about ``n_kphi`` that
    takes ``R(n_kappa, κ) · n_komega`` onto ``R_eul · n_komega``
    (these vectors lie on the same cone about ``n_kphi`` by
    construction).  Finally ``κω`` is the residual rotation about
    ``n_komega`` extracted from
    ``R(n_komega, κω) = R(n_kappa, -κ) · R(n_kphi, -κφ) · R_eul``.

    Parameters
    ----------
    omega_deg, chi_deg, phi_deg : float
        Virtual Eulerian pseudoangles in degrees.
    convention : KappaPseudoAngleConvention
        Per-geometry axis declaration; see
        :func:`make_kappa_pseudo_angle_convention`.
    branch : {+1, -1}, optional
        Branch selection.  Each branch produces a distinct kappa
        triple satisfying the same Eulerian rotation.  ``+1`` (default)
        selects the branch with the larger ``cos(κ)`` (equivalently
        the smaller ``|κ|``); ``-1`` selects the other.

    Returns
    -------
    komega, kappa, kphi : tuple of float
        Real kappa motor angles in degrees, in the cut-point range
        (-180, 180].

    Raises
    ------
    ValueError
        If ``branch`` is not ``+1`` or ``-1``.
    ValueError
        If the requested Eulerian orientation is not reachable by
        any kappa triple in the given convention (the equation in κ
        has no real solution).

    Notes
    -----
    For ``kappa4cv`` (BL convention), the convention is

        n_komega = -TRANSVERSE,   n_kappa = kappa_axis(α),
        n_kphi   = -TRANSVERSE,   n_chi_eq = +LONGITUDINAL.

    Confirming the geometry-aware result against Walko (2016) eq. [16]
    requires using the textbook convention exactly:

        n_komega = +VERTICAL,     n_kappa  = +VERTICAL·cos(α) + +TRANSVERSE·sin(α),
        n_kphi   = +VERTICAL,     n_chi_eq = +LONGITUDINAL.

    See :func:`eulerian_to_kappa` for the textbook closed form.

    References
    ----------
    * Issue #241 — derivation and motivation.
    * Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016) — textbook
      special case.
    """
    if branch not in (+1, -1):
        raise ValueError(
            f"eulerian_to_kappa_axes: branch must be +1 or -1; got {branch!r}."
        )

    # Special case: χ = 0 reduces the kappa arm to the identity rotation
    # (κ = 0).  When ``n_komega`` and ``n_kphi`` are parallel — as they
    # are in every demo geometry shipped with the package — the Z matrix
    # ``R(n_kphi, φ) · R(n_komega, ω)`` collapses to a single rotation
    # ``R(n_komega, ω + φ·sign)`` and the (κω, κφ) split is
    # one-parameter degenerate.  Pick the natural assignment
    # ``(κω, κφ) = (ω, φ)`` so the round-trip
    # ``kappa → eulerian → kappa`` is the identity at this point.
    if abs(chi_deg) < 1e-12:
        return _wrap_180(omega_deg), 0.0, _wrap_180(phi_deg)

    R_eul = _eulerian_rotation(convention, omega_deg, chi_deg, phi_deg)

    n_om = convention.n_komega
    n_ka = convention.n_kappa
    n_ph = convention.n_kphi

    # Step 1 — solve for κ from
    #   n_ph · R(n_ka, κ) · n_om  =  n_ph · R_eul · n_om
    v_target = R_eul @ n_om
    gamma = float(np.dot(n_ph, v_target))

    # Rodrigues expansion of R(n_ka, κ) · n_om:
    #   R(n_ka, κ) · n_om = n_om·cos(κ) + (n_ka × n_om)·sin(κ)
    #                       + n_ka·(n_ka·n_om)·(1 − cos(κ))
    A = float(np.dot(n_ph, n_om))
    B = float(np.dot(n_ph, np.cross(n_ka, n_om)))
    C = float(np.dot(n_ph, n_ka)) * float(np.dot(n_ka, n_om))
    # → (A − C)·cos(κ) + B·sin(κ) = γ − C
    kappa_candidates = _solve_alpha_cos_beta_sin_eq_gamma(A - C, B, gamma - C)

    # Order branches: +1 is the candidate with the larger cos(κ)
    # (equivalently the smaller |κ|).  This matches the historic
    # default that the +1 branch yields ``kappa = 0`` when χ = 0.
    kappa_candidates.sort(key=lambda k: -math.cos(math.radians(k)))
    if branch == +1:
        kappa_deg = kappa_candidates[0]
    else:
        kappa_deg = kappa_candidates[-1]

    # Degeneracy at κ = 0: when the kappa rotation is the identity,
    # ``R(n_kphi, κφ) · R(n_komega, κω) = R(n_kphi, φ) · R(n_komega, ω)``
    # has a one-parameter family of solutions whenever ``n_kphi`` is
    # parallel (or anti-parallel) to ``n_komega`` — and *all* kappa
    # demo geometries in this package have that property by design (komega
    # and kphi share the same physical motor axis).  Pick the
    # representative that maps ``ω → κω`` and ``φ → κφ`` directly so
    # the round-trip with :func:`kappa_to_eulerian_axes` is identity.
    parallel_axes = abs(abs(float(np.dot(n_om, n_ph))) - 1.0) < 1e-12
    if abs(kappa_deg) < 1e-12 and parallel_axes:  # pragma: no cover
        return _wrap_180(omega_deg), 0.0, _wrap_180(phi_deg)

    # Step 2 — recover κφ as the signed rotation about n_kphi taking
    # R(n_kappa, κ) · n_komega onto R_eul · n_komega.
    R_ka = _rotation_matrix_normalized(n_ka, kappa_deg)
    w = R_ka @ n_om
    kphi_deg = _angle_about_axis(w, v_target, n_ph)

    # Step 3 — recover κω as the residual rotation about n_komega.
    # From R(n_kphi, κφ) · R(n_kappa, κ) · R(n_komega, κω) = R_eul,
    # we have R(n_komega, κω) = R(n_kappa, -κ) · R(n_kphi, -κφ) · R_eul.
    R_ph_inv = _rotation_matrix_normalized(n_ph, -kphi_deg)
    R_ka_inv = _rotation_matrix_normalized(n_ka, -kappa_deg)
    R_residual = R_ka_inv @ R_ph_inv @ R_eul
    # R_residual should be a pure rotation about n_komega.  Recover
    # the angle from the action on any vector perpendicular to n_om.
    if abs(n_om[0]) < 0.9:
        ref = np.array([1.0, 0.0, 0.0])
    else:
        ref = np.array([0.0, 1.0, 0.0])
    perp = ref - np.dot(ref, n_om) * n_om
    perp = perp / np.linalg.norm(perp)
    perp_rot = R_residual @ perp
    komega_deg = _angle_about_axis(perp, perp_rot, n_om)

    return _wrap_180(komega_deg), _wrap_180(kappa_deg), _wrap_180(kphi_deg)


def kappa_to_eulerian_axes(
    komega_deg: float,
    kappa_deg: float,
    kphi_deg: float,
    convention: KappaPseudoAngleConvention,
) -> tuple[float, float, float]:
    """
    Convert real kappa motor angles to virtual Eulerian pseudoangles
    using the geometry's actual signed stage axes.

    Inverts :func:`eulerian_to_kappa_axes`.  Decomposes the kappa
    sample-rotation matrix (with stage order matching
    :func:`~orientation._compute_q_phi`) as

        R(n_kphi, κφ) · R(n_kappa, κ) · R(n_komega, κω)
            =  R(n_kphi, φ) · R(n_chi_eq, χ) · R(n_komega, ω)

    and returns ``(ω, χ, φ)``.  The decomposition has two χ branches;
    this function returns the branch consistent with the input kappa
    sign.

    Algorithm
    ---------
    Let ``R = R(n_kphi, κφ) · R(n_kappa, κ) · R(n_komega, κω)`` be
    the realized rotation.  By the same reasoning as in
    :func:`eulerian_to_kappa_axes` (applied to the Eulerian-equivalent
    triple):

    1. Solve for χ from
       ``n_kphi · R(n_chi_eq, χ) · n_komega  =  n_kphi · R · n_komega``.

    2. Recover φ as the signed rotation about ``n_kphi`` taking
       ``R(n_chi_eq, χ) · n_komega`` onto ``R · n_komega``.

    3. Recover ω from the residual rotation about ``n_komega``.

    Parameters
    ----------
    komega_deg, kappa_deg, kphi_deg : float
        Real kappa motor angles in degrees.
    convention : KappaPseudoAngleConvention

    Returns
    -------
    omega, chi, phi : tuple of float
        Virtual Eulerian pseudoangles in degrees, in the cut-point
        range (-180, 180].  The sign of χ is chosen consistently with
        the sign of ``kappa_deg``.

    Notes
    -----
    The chi-sign convention here matches the Walko (2016) sign
    convention: positive κ gives positive χ in the textbook frame.
    For other conventions (e.g. left-handed komega), the sign carries
    through naturally because the decomposition uses the geometry's
    own signed axes.
    """
    # Special case: κ = 0 ⇒ χ = 0 ⇒ degenerate (ω, φ) split.  In every
    # demo geometry shipped with the package ``n_komega`` and ``n_kphi`` are
    # parallel, so any pair satisfying ``ω + φ = κω + κφ`` is a valid
    # decomposition.  The natural assignment ``ω = κω, φ = κφ`` makes
    # the round-trip ``eulerian → kappa → eulerian`` an identity at
    # this point.
    if abs(kappa_deg) < 1e-12:
        return _wrap_180(komega_deg), 0.0, _wrap_180(kphi_deg)

    R = _kappa_rotation(convention, komega_deg, kappa_deg, kphi_deg)

    n_om = convention.n_komega
    n_ch = convention.n_chi_eq
    n_ph = convention.n_kphi

    # Step 1 — solve for χ from
    #   n_ph · R(n_ch, χ) · n_om  =  n_ph · R · n_om
    v_target = R @ n_om
    gamma = float(np.dot(n_ph, v_target))
    A = float(np.dot(n_ph, n_om))
    B = float(np.dot(n_ph, np.cross(n_ch, n_om)))
    C = float(np.dot(n_ph, n_ch)) * float(np.dot(n_ch, n_om))
    chi_candidates = _solve_alpha_cos_beta_sin_eq_gamma(A - C, B, gamma - C)

    # The two χ candidates differ by reflection through the cone-axis
    # plane.  Pick the one that, when fed back through
    # ``eulerian_to_kappa_axes``, reproduces the input motor triple
    # under the branch implied by ``sign(kappa_deg)``.
    best = None
    best_err = float("inf")
    target_branch = +1 if kappa_deg >= 0 else -1
    for chi_deg in chi_candidates:
        # Recover φ from the n_kphi cone identity.
        R_ch = _rotation_matrix_normalized(n_ch, chi_deg)
        w = R_ch @ n_om
        phi_deg = _angle_about_axis(w, v_target, n_ph)
        # Recover ω from the residual rotation about n_komega.
        R_ph_inv = _rotation_matrix_normalized(n_ph, -phi_deg)
        R_ch_inv = _rotation_matrix_normalized(n_ch, -chi_deg)
        R_residual = R_ch_inv @ R_ph_inv @ R
        if abs(n_om[0]) < 0.9:
            ref = np.array([1.0, 0.0, 0.0])
        else:
            ref = np.array([0.0, 1.0, 0.0])
        perp = ref - np.dot(ref, n_om) * n_om
        perp = perp / np.linalg.norm(perp)
        perp_rot = R_residual @ perp
        omega_deg = _angle_about_axis(perp, perp_rot, n_om)
        # Round-trip error: forward through eulerian_to_kappa_axes.
        try:
            ko, k, kp = eulerian_to_kappa_axes(
                omega_deg, chi_deg, phi_deg, convention, branch=target_branch
            )
        except ValueError:  # pragma: no cover
            continue
        err = (
            (ko - _wrap_180(komega_deg)) ** 2
            + (k - _wrap_180(kappa_deg)) ** 2
            + (kp - _wrap_180(kphi_deg)) ** 2
        )
        if err < best_err:
            best_err = err
            best = (omega_deg, chi_deg, phi_deg)

    if best is None:  # pragma: no cover
        raise ValueError(
            "kappa_to_eulerian_axes: no consistent χ branch found for "
            f"(komega={komega_deg}, kappa={kappa_deg}, kphi={kphi_deg})."
        )
    omega_deg, chi_deg, phi_deg = best
    return _wrap_180(omega_deg), _wrap_180(chi_deg), _wrap_180(phi_deg)


def _wrap_180(angle_deg: float) -> float:
    """Wrap an angle in degrees to the half-open interval (-180, 180]."""
    a = (float(angle_deg) + 180.0) % 360.0 - 180.0
    # The modulo yields [-180, 180); shift the −180 endpoint to +180
    # so the wrap is symmetric around 0.
    if a == -180.0:
        return 180.0
    return a


# ---------------------------------------------------------------------------
# Walko (2016) eq. [16] textbook formula — backward-compatible
# ---------------------------------------------------------------------------


def kappa_to_eulerian(
    komega: float,
    kappa: float,
    kphi: float,
    alpha_deg: float = 50.0,
) -> tuple[float, float, float]:
    """
    Convert real kappa angles to virtual Eulerian pseudoangles using
    the closed-form Walko (2016) eq. [16] formula.

    .. warning::

       This function implements Walko's published algebra **literally**
       and is correct only for the axis convention assumed in that
       paper (omega and phi about the transverse axis with the
       textbook handedness; chi about the longitudinal axis).  No
       demo geometry in this package matches that convention exactly.  For
       any geometry-dependent calculation, use
       :func:`kappa_to_eulerian_axes` with the demo geometry's
       :class:`KappaPseudoAngleConvention` instead.

       Issue #241 records the diagnostic that uncovered this
       limitation.

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
        Virtual Eulerian pseudoangles in degrees in Walko's frame.

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
            f"kappa_to_eulerian: chi approaches ±180° "
            f"(chi_half={math.degrees(chi_half):.4f}°) — singularity."
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
    Convert virtual Eulerian pseudoangles to real kappa angles using
    the inverse of the closed-form Walko (2016) eq. [16] formula.

    .. warning::

       This function implements Walko's published algebra **literally**
       and is correct only for the textbook axis convention; see the
       warning in :func:`kappa_to_eulerian` and issue #241.  For any
       geometry-dependent calculation, use
       :func:`eulerian_to_kappa_axes` with the demo geometry's
       :class:`KappaPseudoAngleConvention` instead.

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
        Real kappa motor angles in degrees in Walko's frame.

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


def _build_eulerian_equivalent_geometry(
    geometry: AdHocDiffractometer,
) -> tuple[AdHocDiffractometer, str, str, str]:
    """
    Build a synthetic Eulerian-equivalent geometry from a kappa
    geometry.

    The synthetic geometry has identical outer and inner stages but
    replaces the kappa triple ``(komega, kappa, kphi)`` with the
    Eulerian triple ``(omega, chi, phi)`` whose axes come from the
    kappa demo geometry's :class:`KappaPseudoAngleConvention`.  The result
    is a standard Eulerian geometry that the existing solvers handle
    natively.

    All non-kappa stage angles (and limits) are copied verbatim from
    the source geometry.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Kappa geometry with ``kappa_pseudo_angle_convention`` set.

    Returns
    -------
    eul_geom : AdHocDiffractometer
        Synthetic Eulerian equivalent.  Has no modes / sample / UB
        attached; intended for direct use with the low-level forward
        dispatcher.
    komega_name, kphi_name : str
        Names of the outermost and innermost kappa stages on the
        original geometry.  Returned for convenience so callers can
        substitute Eulerian solutions back into the kappa motor dict.
    """
    from .diffractometer import AdHocDiffractometer
    from .stage import Stage

    convention = geometry.kappa_pseudo_angle_convention
    sample_stages = geometry.sample_stages
    stage_names = [s.name for s in sample_stages]
    kappa_idx = next(i for i, s in enumerate(sample_stages) if s.name == "kappa")

    komega_name = stage_names[kappa_idx - 1]
    kphi_name = stage_names[kappa_idx + 1]

    eul_stages: list[Stage] = []
    for i, s in enumerate(sample_stages):
        if i == kappa_idx - 1:
            eul_stages.append(
                Stage(
                    "omega",
                    convention.n_komega,
                    parent=s.parent,
                    role="sample",
                    angle=s.angle,
                    limits=s.limits,
                )
            )
        elif i == kappa_idx:
            eul_stages.append(
                Stage(
                    "chi",
                    convention.n_chi_eq,
                    parent="omega",
                    role="sample",
                    angle=0.0,
                    limits=s.limits,
                )
            )
        elif i == kappa_idx + 1:
            eul_stages.append(
                Stage(
                    "phi",
                    convention.n_kphi,
                    parent="chi",
                    role="sample",
                    angle=s.angle,
                    limits=s.limits,
                )
            )
        else:
            # Outer / inner non-kappa sample stages (e.g. mu on kappa6c)
            eul_stages.append(
                Stage(
                    s.name,
                    s.axis,
                    parent=s.parent,
                    role="sample",
                    angle=s.angle,
                    limits=s.limits,
                )
            )

    # Detector stages copied verbatim
    for s in geometry.detector_stages:
        eul_stages.append(
            Stage(
                s.name,
                s.axis,
                parent=s.parent,
                role="detector",
                angle=s.angle,
                limits=s.limits,
            )
        )

    eul_geom = AdHocDiffractometer(
        name=f"{geometry.name}__eulerian_equivalent",
        stages=eul_stages,
        basis=dict(geometry.basis),
        wavelength=geometry.wavelength,
    )
    return eul_geom, komega_name, kphi_name


def _translate_mode_to_eulerian(
    mode,
    det_stage_name: str,
):
    """Rewrite a kappa virtual-angle mode onto the Eulerian-equivalent
    geometry.

    Virtual-angle ``SampleConstraint`` entries on ``omega``/``chi``/
    ``phi`` already speak the Eulerian language and are passed through.
    All other constraints (fixed real outer/inner stages, detector
    constraints, bisect constraints) are preserved verbatim.  A
    ``VirtualBisectConstraint`` is downgraded to a plain
    ``BisectConstraint`` on the equivalent ``omega`` stage so the
    standard analytic Eulerian bisect solver is dispatched.
    """
    from .mode import BisectConstraint
    from .mode import ConstraintSet
    from .mode import SampleConstraint
    from .mode import VirtualBisectConstraint

    new_constraints = []
    for c in mode.constraints:
        if isinstance(c, VirtualBisectConstraint):
            new_constraints.append(
                BisectConstraint(sample_stage="omega", detector_stage=c.detector_stage)
            )
        elif isinstance(c, SampleConstraint) and c.name in KAPPA_VIRTUAL_ANGLES:
            new_constraints.append(SampleConstraint(c.name, c.value))
        else:
            new_constraints.append(c)
    return ConstraintSet(
        new_constraints,
        computed=["omega", "chi", "phi", det_stage_name],
    )


def solve_kappa_virtual(
    geometry: AdHocDiffractometer,
    Q_phi,
    ttheta_deg: float,
    mode,
) -> list[dict[str, float]]:
    """
    Solve a kappa virtual-angle mode: find real kappa motor angles
    satisfying the virtual Eulerian-angle constraints and the Bragg
    condition.

    The solver uses the **geometry-aware** decomposition
    (:func:`eulerian_to_kappa_axes`) introduced by issue #241.  The
    flow is:

    1. Build a synthetic Eulerian-equivalent geometry whose three
       sample stages have axes ``(n_komega, n_chi_eq, n_kphi)`` taken
       from the kappa demo geometry's
       :class:`KappaPseudoAngleConvention`.
    2. Translate the kappa virtual-angle mode onto the equivalent
       geometry (a ``VirtualBisectConstraint`` becomes a plain
       ``BisectConstraint`` on the synthetic ``omega`` stage; fixed
       virtual ``SampleConstraint`` entries pass through verbatim).
    3. Dispatch into the standard Eulerian forward solver
       (``_solve_constraint_set``), which handles the bisect or
       fixed-angle variant analytically.
    4. Convert each Eulerian motor triple to kappa motors via
       :func:`eulerian_to_kappa_axes` for both branches ±1, with
       deduplication.

    Each step is closed-form; no Newton iteration is required.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        Kappa geometry with ``kappa_pseudo_angle_convention`` set.
    Q_phi : array-like, shape (3,)
        Scattering vector in the phi frame (Å⁻¹).
    ttheta_deg : float
        Detector angle (2θ) in degrees.
    mode : ConstraintSet
        Either a virtual-bisect mode (``VirtualBisectConstraint``) or
        a fixed-virtual-angle mode (``SampleConstraint`` on one of
        ``omega``, ``chi``, ``phi``).

    Returns
    -------
    list of dict[str, float]
        One angles dict per solution, keyed by real stage names.
        Cut-points, limits, and deduplication on the kappa side are
        applied by the caller in :mod:`~ad_hoc_diffractometer.forward`.
    """
    from .mode import DetectorConstraint
    from .mode import SampleConstraint

    convention = geometry.kappa_pseudo_angle_convention
    if convention is None:  # pragma: no cover
        return []

    det_stage = geometry.detector_stages[-1]

    Q_phi_arr = np.asarray(Q_phi, dtype=float)

    eul_geom, komega_name, kphi_name = _build_eulerian_equivalent_geometry(geometry)

    # Apply the original mode's fixed sample/detector constraints to
    # the equivalent geometry's live stage angles, so the synthetic
    # forward computation sees the right baseline.
    for c in mode.constraints:
        if isinstance(c, SampleConstraint) and c.name in eul_geom._stages:  # noqa: SLF001
            eul_geom.stage(c.name).angle = float(c.value)
    if isinstance(getattr(mode, "detector_constraint", None), DetectorConstraint):
        dc = mode.detector_constraint
        # The qaz-detector branch is never reached from a kappa
        # virtual-angle mode in any shipped demo geometry (qaz modes route
        # through ``_solve_qaz_mode``, not through here).
        if (
            not getattr(dc, "is_qaz", False) and dc.name in eul_geom._stages  # noqa: SLF001
        ):  # pragma: no branch
            eul_geom.stage(dc.name).angle = float(dc.value)

    eul_mode = _translate_mode_to_eulerian(mode, det_stage.name)

    # Install the synthetic mode and dispatch.  The equivalent geometry
    # has no UB / sample attached, but _solve_constraint_set works in
    # the phi-frame and only needs Q_phi + ttheta + mode.
    eul_geom._modes._data["__synthetic__"] = eul_mode  # noqa: SLF001
    eul_geom._mode_name = "__synthetic__"  # noqa: SLF001

    from .forward import _solve_constraint_set

    eul_solutions = _solve_constraint_set(eul_geom, Q_phi_arr, ttheta_deg, eul_mode)

    # Build the kappa-side baseline angles.
    base_angles: dict[str, float] = {
        s.name: s.angle
        for s in list(geometry._stages.values())  # noqa: SLF001
    }
    base_angles[det_stage.name] = ttheta_deg
    for c in mode.constraints:
        if isinstance(c, SampleConstraint) and c.name in geometry._stages:  # noqa: SLF001
            base_angles[c.name] = float(c.value)
    if isinstance(getattr(mode, "detector_constraint", None), DetectorConstraint):
        dc = mode.detector_constraint
        # See note above: the qaz branch is unreachable from a
        # virtual-angle mode in shipped demo geometries.
        if (
            not getattr(dc, "is_qaz", False) and dc.name in geometry._stages  # noqa: SLF001
        ):  # pragma: no branch
            base_angles[dc.name] = float(dc.value)

    results: list[dict[str, float]] = []
    seen: list[tuple[float, float, float]] = []
    for eul_sol in eul_solutions:
        omega_e = float(eul_sol.get("omega", 0.0))
        chi_e = float(eul_sol.get("chi", 0.0))
        phi_e = float(eul_sol.get("phi", 0.0))
        for branch in (+1, -1):
            try:
                ko, k, kp = eulerian_to_kappa_axes(
                    omega_e, chi_e, phi_e, convention, branch=branch
                )
            except ValueError:
                continue
            sig = (round(ko, 6), round(k, 6), round(kp, 6))
            if any(
                abs(sig[0] - s[0]) < 1e-6
                and abs(sig[1] - s[1]) < 1e-6
                and abs(sig[2] - s[2]) < 1e-6
                for s in seen
            ):
                continue
            seen.append(sig)
            angles = dict(base_angles)
            angles[komega_name] = ko
            angles["kappa"] = k
            angles[kphi_name] = kp
            results.append(angles)

    return results
