# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
refinement.py — Lattice and orientation refinement from multiple reflections.

Functions
---------
refine_lattice_bl1967(sample, reflections, ...)
    Simultaneously refine cell parameters and orientation (U matrix) from
    N ≥ 3 reflections by iterative least-squares (Busing & Levy 1967,
    §"Refinement of lattice and orientation parameters").  Uses finite-
    difference Jacobians; no scipy dependency.

refine_lattice_simplex(sample, reflections, ...)
    Derivative-free refinement of cell parameters (and optionally orientation)
    using the Nelder-Mead simplex method.  Falls back to a pure-numpy
    implementation when scipy is not available.

Parameter-space notes
---------------------
Both functions accept a ``refine_all`` keyword argument (default ``False``).

* ``refine_all=False`` — refine only the **free parameters** for the
  crystal system currently deduced from ``sample.lattice``.  Symmetry
  constraints are maintained at every iteration: constrained parameters
  are re-derived from the free ones by constructing a new ``Lattice``
  object.  For example, for a cubic crystal only ``a`` is refined;
  ``b``, ``c``, and all angles are kept equal to ``a`` and 90°.

* ``refine_all=True`` — treat all six cell parameters as independent,
  regardless of the crystal system of the starting lattice.  Use this
  when you expect the refined lattice to break the nominal symmetry
  (e.g. due to strain or twinning), or when the crystal system is
  genuinely unknown.

The same ``refine_all`` argument applies to both ``refine_lattice_bl1967``
and ``refine_lattice_simplex``.

References
----------
* Busing & Levy, Acta Cryst. 22, 457-464 (1967),
  §"Refinement of lattice and orientation parameters".
* Nelder & Mead, The Computer Journal 7(4), 308-313 (1965).
"""

from __future__ import annotations

import logging
import math

import numpy as np

from .lattice import _SYSTEM_FREE_PARAMS
from .lattice import Lattice
from .orientation import angles_to_phi_vector

logger = logging.getLogger(__name__)

# All six cell parameters, in the order stored internally.
_ALL_CELL_PARAMS = ("a", "b", "c", "alpha", "beta", "gamma")
_ORIENT_LABELS = ("phi1", "phi2", "phi3")

# ---------------------------------------------------------------------------
# Crystal-system constraint helpers
# ---------------------------------------------------------------------------


def _free_params_for_system(system: str) -> tuple[str, ...]:
    """Return the free parameter names for the given crystal system."""
    return _SYSTEM_FREE_PARAMS.get(system, _SYSTEM_FREE_PARAMS["triclinic"])


# Fixed (constrained) parameters that must be supplied to Lattice.__init__
# to disambiguate the crystal system from the free parameters alone.
# For example, hexagonal needs gamma=120 explicitly; without it Lattice
# infers tetragonal from (a, c) alone.
_SYSTEM_HINT_PARAMS: dict[str, dict] = {
    "cubic": {},
    "tetragonal": {},
    "orthorhombic": {},
    "hexagonal": {"gamma": 120.0},
    "trigonal": {},
    "monoclinic": {},
    "triclinic": {},
}


def _full_params_from_free(free_vals: dict, system: str) -> dict:
    """
    Expand the free-parameter dict to a full six-parameter dict by applying
    the symmetry constraints of ``system``.

    Constructs a temporary ``Lattice`` from the free values plus any
    system-specific hint parameters (e.g. ``gamma=120`` for hexagonal)
    needed for the ``Lattice`` class to recognize the correct system, then
    reads back all six parameters.

    Parameters
    ----------
    free_vals : dict
        Mapping from free parameter name → value for the given system.
    system : str
        Crystal system name (key in ``_SYSTEM_FREE_PARAMS``).

    Returns
    -------
    dict with keys 'a','b','c','alpha','beta','gamma'.
    """
    kw: dict = dict(_SYSTEM_HINT_PARAMS.get(system, {}))
    for pname in ("a", "b", "c", "alpha", "beta", "gamma"):
        if pname in free_vals:
            kw[pname] = free_vals[pname]

    lat = Lattice(**kw)
    return {
        "a": lat.a,
        "b": lat.b,
        "c": lat.c,
        "alpha": lat.alpha,
        "beta": lat.beta,
        "gamma": lat.gamma,
    }


def _active_cell_params(lat: Lattice, refine_all: bool) -> tuple[str, ...]:
    """
    Return the parameter names that will be varied during refinement.

    Parameters
    ----------
    lat : Lattice
        Current sample lattice (used to read crystal system when
        ``refine_all=False``).
    refine_all : bool
        If True, all six parameters are free.

    Returns
    -------
    tuple of str — the active parameter names.
    """
    if refine_all:
        return _ALL_CELL_PARAMS
    return _free_params_for_system(lat.system)


# ---------------------------------------------------------------------------
# Reflection resolution
# ---------------------------------------------------------------------------


def _resolve_reflections(sample, reflections: list) -> list:
    """
    Return a list of Reflection objects from a mix of names and objects.
    """
    from .reflection import Reflection

    result = []
    for r in reflections:
        if isinstance(r, str):
            result.append(sample.reflections[r])
        elif isinstance(r, Reflection):
            result.append(r)
        else:
            raise TypeError(
                f"Each reflection must be a Reflection object or a name string; "
                f"got {type(r).__name__!r}."
            )
    return result


# ---------------------------------------------------------------------------
# B-matrix and Jacobian helpers
# ---------------------------------------------------------------------------


def _B_from_full_params(params: dict) -> np.ndarray:
    """Return the B matrix for a complete six-parameter dict."""
    return Lattice(
        a=params["a"],
        b=params["b"],
        c=params["c"],
        alpha=params["alpha"],
        beta=params["beta"],
        gamma=params["gamma"],
    ).B


def _dB_dp_finite_diff(
    full_params: dict,
    free_params: dict,
    pname: str,
    system: str,
    refine_all: bool,
    eps: float = 1e-6,
) -> np.ndarray:
    """
    Finite-difference derivative of the B matrix with respect to free parameter
    ``pname``, enforcing symmetry constraints when ``refine_all=False``.

    Parameters
    ----------
    full_params : dict — current full six-parameter set.
    free_params : dict — current free-parameter values.
    pname : str — the free parameter being differentiated.
    system : str — crystal system.
    refine_all : bool
    eps : float — finite-difference step.

    Returns
    -------
    dB : numpy.ndarray, shape (3, 3)
    """
    fp_plus = dict(free_params)
    fp_plus[pname] += eps
    fp_minus = dict(free_params)
    fp_minus[pname] -= eps

    if refine_all:
        # No constraints: perturb directly in full-param space
        p_plus = dict(full_params)
        p_plus[pname] += eps
        p_minus = dict(full_params)
        p_minus[pname] -= eps
    else:
        p_plus = _full_params_from_free(fp_plus, system)
        p_minus = _full_params_from_free(fp_minus, system)

    Bp = _B_from_full_params(p_plus)
    Bm = _B_from_full_params(p_minus)
    return (Bp - Bm) / (2.0 * eps)


def _rotation_from_vector(v: np.ndarray) -> np.ndarray:
    """Rotation matrix for a small rotation vector v (Rodrigues formula)."""
    angle = np.linalg.norm(v)
    if angle < 1e-15:
        return np.eye(3)
    from .rotation import rotation_matrix

    return rotation_matrix(v / angle, math.degrees(angle))


# ---------------------------------------------------------------------------
# Core least-squares building block
# ---------------------------------------------------------------------------


def _build_jacobian_and_residuals(
    geometry,
    refl_list: list,
    UB: np.ndarray,
    full_params: dict,
    free_params: dict,
    active_cell: tuple[str, ...],
    system: str,
    refine_all: bool,
    refine_cell: bool,
    refine_orientation: bool,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Build Jacobian J (3N × n_params) and residual vector r (3N,).

    Cell-param columns use finite-difference on B; orientation columns use
    the analytical cross-product formula.
    """
    n = len(refl_list)
    n_cell = len(active_cell) if refine_cell else 0
    n_orient = 3 if refine_orientation else 0
    n_params = n_cell + n_orient

    J = np.zeros((3 * n, n_params))
    r = np.zeros(3 * n)

    # U from current UB and B
    try:
        U = UB @ np.linalg.inv(_B_from_full_params(full_params))
    except np.linalg.LinAlgError:  # pragma: no cover
        U = np.eye(3)

    for i, refl in enumerate(refl_list):
        h = np.asarray(refl.hkl, dtype=float)
        h_phi_obs = angles_to_phi_vector(geometry, **refl.angles)
        h_phi_calc = UB @ h
        r[3 * i : 3 * i + 3] = h_phi_obs - h_phi_calc

        col = 0
        if refine_cell:
            for pname in active_cell:
                dBdp = _dB_dp_finite_diff(
                    full_params, free_params, pname, system, refine_all
                )
                J[3 * i : 3 * i + 3, col] = (U @ dBdp) @ h
                col += 1

        if refine_orientation:
            for j in range(3):
                e_j = np.zeros(3)
                e_j[j] = 1.0
                J[3 * i : 3 * i + 3, col] = np.cross(e_j, h_phi_calc)
                col += 1

    return J, r


def _apply_delta(
    full_params: dict,
    free_params: dict,
    UB: np.ndarray,
    delta_p: np.ndarray,
    active_cell: tuple[str, ...],
    system: str,
    refine_all: bool,
    refine_cell: bool,
    refine_orientation: bool,
) -> tuple[dict, dict, np.ndarray]:
    """
    Apply parameter update delta_p; return (new_full_params, new_free_params, new_UB).
    """
    new_free = dict(free_params)
    col = 0

    if refine_cell:
        for pname in active_cell:
            new_free[pname] = free_params[pname] + float(delta_p[col])
            col += 1
        # Clamp free params to physical range
        for pname in ("a", "b", "c"):
            if pname in new_free:
                new_free[pname] = max(new_free[pname], 1e-4)
        for pname in ("alpha", "beta", "gamma"):
            if pname in new_free:
                new_free[pname] = float(np.clip(new_free[pname], 1.1, 178.9))

    # Expand free → full
    if refine_cell:
        if refine_all:
            new_full = {k: new_free.get(k, full_params[k]) for k in _ALL_CELL_PARAMS}
        else:
            new_full = _full_params_from_free(new_free, system)
    else:
        new_full = dict(full_params)

    # Update UB
    B_new = _B_from_full_params(new_full)
    try:
        B_old = _B_from_full_params(full_params)
        U_old = UB @ np.linalg.inv(B_old)
    except np.linalg.LinAlgError:  # pragma: no cover
        U_old = np.eye(3)

    if refine_orientation:
        phi_vec = delta_p[col : col + 3]
        dR = _rotation_from_vector(phi_vec)
        U_new = dR @ U_old
    else:
        U_new = U_old

    new_UB = U_new @ B_new
    return new_full, new_free, new_UB


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def refine_lattice_bl1967(
    sample,
    reflections,
    refine_cell: bool = True,
    refine_orientation: bool = True,
    refine_all: bool = False,
    max_iter: int = 100,
    tol: float | None = None,
) -> dict:
    """
    Refine lattice parameters and/or orientation from N≥3 reflections.

    Uses the Busing & Levy (1967) least-squares approach (§"Refinement of
    lattice and orientation parameters").  Starting from the current
    ``sample.UB`` and ``sample.lattice``, iteratively solves the linearized
    normal equations to minimize the sum of squared residuals between
    observed and calculated phi-frame scattering vectors.

    Observed vector for reflection i:
        ``hiφ_obs = angles_to_phi_vector(geometry, **ri.angles)``
    Calculated vector:
        ``hiφ_calc = UB @ hi``
    Residual:
        ``Δhiφ = hiφ_obs − hiφ_calc``

    The Jacobian of ``hiφ_calc`` with respect to cell parameters is
    computed by finite difference on the B matrix; the Jacobian with
    respect to the three orientation parameters (infinitesimal rotation
    components) uses the analytical cross-product formula.

    Parameters
    ----------
    sample : Sample
        The sample whose ``lattice``, ``U``, and ``UB`` are updated
        in-place.  ``sample.parent`` must be a geometry with
        ``wavelength`` set (needed by ``angles_to_phi_vector``).
        ``sample.UB`` must be set (call one of the ``ub_from_*``
        functions first).
    reflections : list of Reflection or str
        N ≥ 3 reflections.  Each may be a ``Reflection`` object or a
        name string in ``sample.reflections``.
    refine_cell : bool, optional
        If True (default), refine the cell parameters.
    refine_orientation : bool, optional
        If True (default), refine the three orientation parameters
        (infinitesimal rotation components).
    refine_all : bool, optional
        Controls which cell parameters are treated as free:

        * ``False`` (default) — refine only the parameters that are
          free in the crystal system of the *current* ``sample.lattice``
          (e.g. only ``a`` for cubic, ``a`` and ``c`` for tetragonal,
          etc.).  Symmetry constraints are enforced at every iteration.
        * ``True`` — treat all six parameters (a, b, c, α, β, γ) as
          independent, regardless of the starting crystal system.  Use
          this when the refined cell is expected to break the nominal
          symmetry, or when the crystal system is unknown.

    max_iter : int, optional
        Maximum number of iterations (default 100).
    tol : float or None, optional
        Convergence tolerance on the RMS residual change between
        iterations.  If None, defaults to 1e-10 Å⁻¹.

    Returns
    -------
    result : dict
        ``lattice``    — refined :class:`Lattice` (also stored on ``sample``)
        ``UB``         — refined UB matrix (3×3 ndarray)
        ``U``          — derived U = UB @ B⁻¹ (3×3 ndarray or None)
        ``residuals``  — per-reflection residual vectors, shape (N, 3) Å⁻¹
        ``rms``        — root-mean-square residual (scalar, Å⁻¹)
        ``converged``  — True if converged within ``max_iter`` iterations
        ``n_iter``     — number of iterations performed

    Raises
    ------
    ValueError
        If ``sample.parent`` is None.
    ValueError
        If ``sample.UB`` is None (call a ``ub_from_*`` function first).
    ValueError
        If fewer than 3 reflections are supplied.
    ValueError
        If neither ``refine_cell`` nor ``refine_orientation`` is True.
    TypeError
        If any element of ``reflections`` is not a Reflection or str.

    Notes
    -----
    UB is updated in-place at each iteration.  U is derived at the end
    as ``U = UB @ B⁻¹``.  Both ``sample.U`` and ``sample.UB`` are set
    in-place before returning, and ``sample.lattice`` is replaced by
    the refined ``Lattice`` object.

    If the least-squares system is ill-conditioned (e.g. reflections are
    nearly coplanar), the normal equations are solved with
    ``numpy.linalg.lstsq`` (rcond=1e-15), which returns the minimum-norm
    solution.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.fourcv()
    >>> g.wavelength = 1.5406
    >>> g.add_sample("sapphire", ahd.Lattice(a=4.758, c=12.991, gamma=120))
    >>> g.sample = "sapphire"
    >>> # ... add reflections, compute initial UB ...
    >>> result = ahd.refine_lattice_bl1967(
    ...     g.sample, ["r1", "r2", "r3"],
    ...     refine_all=False,   # keep hexagonal constraints
    ... )
    >>> result["converged"]
    True
    >>> result["lattice"].a    # refined a (= b enforced by hexagonal symmetry)
    ...

    References
    ----------
    Busing & Levy, Acta Cryst. 22, 457-464 (1967),
    §"Refinement of lattice and orientation parameters".
    """
    if sample.parent is None:
        raise ValueError(
            "refine_lattice_bl1967 requires sample.parent to be set "
            "(an AdHocDiffractometer with wavelength)."
        )
    if sample.UB is None:
        raise ValueError(
            "refine_lattice_bl1967 requires sample.UB to be set. "
            "Call ub_from_two_reflections_bl1967() or similar first."
        )
    if not refine_cell and not refine_orientation:
        raise ValueError(
            "At least one of refine_cell or refine_orientation must be True."
        )

    refl_list = _resolve_reflections(sample, reflections)
    if len(refl_list) < 3:
        raise ValueError(
            f"refine_lattice_bl1967 requires at least 3 reflections; "
            f"got {len(refl_list)}."
        )

    if tol is None:
        tol = 1e-10

    geometry = sample.parent
    lat = sample.lattice
    system = lat.system
    active_cell = _active_cell_params(lat, refine_all)

    # Build current full and free parameter dicts
    full_params = {
        "a": lat.a,
        "b": lat.b,
        "c": lat.c,
        "alpha": lat.alpha,
        "beta": lat.beta,
        "gamma": lat.gamma,
    }
    free_params = {k: full_params[k] for k in active_cell}

    UB = sample.UB.copy()
    rms_prev = np.inf
    converged = False
    n_iter = 0

    for iteration in range(max_iter):
        n_iter = iteration + 1
        J, r = _build_jacobian_and_residuals(
            geometry,
            refl_list,
            UB,
            full_params,
            free_params,
            active_cell,
            system,
            refine_all,
            refine_cell,
            refine_orientation,
        )

        rms = float(np.sqrt(np.mean(r**2)))
        if abs(rms_prev - rms) < tol:
            converged = True
            break
        rms_prev = rms

        delta_p, _, _, _ = np.linalg.lstsq(J, r, rcond=1e-15)
        full_params, free_params, UB = _apply_delta(
            full_params,
            free_params,
            UB,
            delta_p,
            active_cell,
            system,
            refine_all,
            refine_cell,
            refine_orientation,
        )

    # Final residuals
    _, r_final = _build_jacobian_and_residuals(
        geometry,
        refl_list,
        UB,
        full_params,
        free_params,
        active_cell,
        system,
        refine_all,
        refine_cell,
        refine_orientation,
    )
    n = len(refl_list)
    residuals = r_final.reshape(n, 3)
    rms_final = float(np.sqrt(np.mean(r_final**2)))

    # Build refined Lattice
    try:
        refined_lattice = Lattice(
            a=full_params["a"],
            b=full_params["b"],
            c=full_params["c"],
            alpha=full_params["alpha"],
            beta=full_params["beta"],
            gamma=full_params["gamma"],
        )
    except ValueError:  # pragma: no cover
        refined_lattice = sample.lattice

    B_refined = refined_lattice.B
    try:
        U_refined = UB @ np.linalg.inv(B_refined)
    except np.linalg.LinAlgError:  # pragma: no cover
        U_refined = None

    sample.lattice = refined_lattice
    sample.UB = UB
    sample.U = U_refined

    return {
        "lattice": refined_lattice,
        "UB": UB,
        "U": U_refined,
        "residuals": residuals,
        "rms": rms_final,
        "converged": converged,
        "n_iter": n_iter,
    }


def refine_lattice_simplex(
    sample,
    reflections,
    refine_cell: bool = True,
    refine_orientation: bool = False,
    refine_all: bool = False,
    max_iter: int = 1000,
    tol: float | None = None,
) -> dict:
    """
    Derivative-free lattice refinement using the Nelder-Mead simplex method.

    Minimises the scalar cost function::

        f(params) = Σ_i |hiφ_obs_i − UB(params) @ hi|²

    where ``UB(params) = U @ B(params)`` and U is held fixed at its current
    value (unless ``refine_orientation=True``).

    Uses ``scipy.optimize.minimize`` with ``method='Nelder-Mead'`` when
    scipy is available, and falls back to a built-in pure-numpy
    Nelder-Mead implementation otherwise.

    Parameters
    ----------
    sample : Sample
        The sample whose ``lattice``, ``U``, and ``UB`` are updated
        in-place.  ``sample.parent`` must be a geometry with
        ``wavelength`` set.  ``sample.UB`` must already be set.
    reflections : list of Reflection or str
        N ≥ 3 reflections.
    refine_cell : bool, optional
        If True (default), refine the cell parameters.
    refine_orientation : bool, optional
        If True, also perturb the orientation parameters (three rotation
        components) during the simplex search.  Default False.
    refine_all : bool, optional
        Controls which cell parameters are treated as free:

        * ``False`` (default) — refine only the free parameters for the
          crystal system of the *current* ``sample.lattice``.  Symmetry
          constraints are enforced by reconstructing a ``Lattice`` at
          each function evaluation.
        * ``True`` — treat all six parameters as independent.

    max_iter : int, optional
        Maximum number of function evaluations (default 1000).
    tol : float or None, optional
        Convergence tolerance on the cost function value.
        If None, defaults to 1e-12.

    Returns
    -------
    result : dict
        Same keys as ``refine_lattice_bl1967``:
        ``lattice``, ``UB``, ``U``, ``residuals``, ``rms``,
        ``converged``, ``n_iter``.

    Raises
    ------
    ValueError
        If ``sample.parent`` is None or ``sample.UB`` is None.
    ValueError
        If fewer than 3 reflections are supplied.
    ValueError
        If neither ``refine_cell`` nor ``refine_orientation`` is True.
    TypeError
        If any element of ``reflections`` is not a Reflection or str.

    Notes
    -----
    The default ``refine_orientation=False`` leaves U fixed; only the
    cell parameters (and thus B, and thus UB = U @ B) are optimized.
    This is the typical use case: use ``ub_from_two_reflections_bl1967``
    to get U, then use the simplex to find the best cell parameters.

    When scipy is not available, a minimal pure-numpy Nelder-Mead
    implementation is used.

    Examples
    --------
    >>> result = ahd.refine_lattice_simplex(
    ...     g.sample, ["r1", "r2", "r3"],
    ...     refine_all=False,   # keep hexagonal constraints on a=b, gamma=120
    ... )
    >>> result["converged"]
    True

    References
    ----------
    Nelder & Mead, The Computer Journal 7(4), 308-313 (1965).
    """
    if sample.parent is None:
        raise ValueError("refine_lattice_simplex requires sample.parent to be set.")
    if sample.UB is None:
        raise ValueError(
            "refine_lattice_simplex requires sample.UB to be set. "
            "Call ub_from_two_reflections_bl1967() or similar first."
        )
    if not refine_cell and not refine_orientation:
        raise ValueError(
            "At least one of refine_cell or refine_orientation must be True."
        )

    refl_list = _resolve_reflections(sample, reflections)
    if len(refl_list) < 3:
        raise ValueError(
            f"refine_lattice_simplex requires at least 3 reflections; "
            f"got {len(refl_list)}."
        )

    if tol is None:
        tol = 1e-12

    geometry = sample.parent
    lat = sample.lattice
    system = lat.system
    active_cell = _active_cell_params(lat, refine_all)

    full_params0 = {
        "a": lat.a,
        "b": lat.b,
        "c": lat.c,
        "alpha": lat.alpha,
        "beta": lat.beta,
        "gamma": lat.gamma,
    }
    free_params0 = {k: full_params0[k] for k in active_cell}

    UB0 = sample.UB.copy()
    B0 = _B_from_full_params(full_params0)
    try:
        U_fixed = UB0 @ np.linalg.inv(B0)
    except np.linalg.LinAlgError:  # pragma: no cover
        U_fixed = np.eye(3)

    # Build initial parameter vector for the simplex
    p0_cell = np.array([free_params0[k] for k in active_cell])
    if refine_cell and refine_orientation:
        p0 = np.concatenate([p0_cell, np.zeros(3)])
    elif refine_cell:
        p0 = p0_cell.copy()
    else:
        p0 = np.zeros(3)

    # Pre-compute observations
    h_phi_obs = np.array(
        [angles_to_phi_vector(geometry, **r.angles) for r in refl_list]
    )
    h_vecs = np.array([r.hkl for r in refl_list], dtype=float)

    def cost(p: np.ndarray) -> float:
        if refine_cell and refine_orientation:
            p_cell = p[: len(active_cell)]
            phi_vec = p[len(active_cell) : len(active_cell) + 3]
        elif refine_cell:
            p_cell = p
            phi_vec = np.zeros(3)
        else:
            p_cell = p0_cell
            phi_vec = p

        free_vals = {k: float(p_cell[i]) for i, k in enumerate(active_cell)}
        # Clamp lengths and angles
        for pname in ("a", "b", "c"):
            if pname in free_vals:
                free_vals[pname] = max(free_vals[pname], 1e-4)
        for pname in ("alpha", "beta", "gamma"):
            if pname in free_vals:
                free_vals[pname] = float(np.clip(free_vals[pname], 1.1, 178.9))

        if refine_all:
            full = {k: free_vals.get(k, full_params0[k]) for k in _ALL_CELL_PARAMS}
        else:
            try:
                full = _full_params_from_free(free_vals, system)
            except Exception:  # pragma: no cover
                return 1e30

        try:
            B = _B_from_full_params(full)
        except Exception:  # pragma: no cover
            return 1e30

        dR = _rotation_from_vector(phi_vec)
        U = dR @ U_fixed
        UB = U @ B
        diff = h_phi_obs - (UB @ h_vecs.T).T
        return float(np.sum(diff**2))

    # Minimise
    n_iter = 0
    converged = False
    p_opt = p0.copy()

    try:
        from scipy.optimize import minimize  # type: ignore[import]

        options = {
            "maxiter": max_iter,
            "xatol": tol**0.5,
            "fatol": tol,
        }
        res = minimize(cost, p0, method="Nelder-Mead", options=options)
        p_opt = res.x
        n_iter = res.nit
        converged = res.success

    except ImportError:
        p_opt, n_iter, converged = _nelder_mead_numpy(cost, p0, max_iter, tol)

    # Reconstruct refined UB from optimal parameters
    if refine_cell and refine_orientation:
        p_cell_opt = p_opt[: len(active_cell)]
        phi_vec_opt = p_opt[len(active_cell) : len(active_cell) + 3]
    elif refine_cell:
        p_cell_opt = p_opt
        phi_vec_opt = np.zeros(3)
    else:
        p_cell_opt = p0_cell
        phi_vec_opt = p_opt

    free_opt = {k: float(p_cell_opt[i]) for i, k in enumerate(active_cell)}
    for pname in ("a", "b", "c"):
        if pname in free_opt:
            free_opt[pname] = max(free_opt[pname], 1e-4)
    for pname in ("alpha", "beta", "gamma"):
        if pname in free_opt:
            free_opt[pname] = float(np.clip(free_opt[pname], 1.1, 178.9))

    if refine_all:
        full_opt = {k: free_opt.get(k, full_params0[k]) for k in _ALL_CELL_PARAMS}
    else:
        try:
            full_opt = _full_params_from_free(free_opt, system)
        except Exception:  # pragma: no cover
            full_opt = full_params0

    try:
        refined_lattice = Lattice(
            a=full_opt["a"],
            b=full_opt["b"],
            c=full_opt["c"],
            alpha=full_opt["alpha"],
            beta=full_opt["beta"],
            gamma=full_opt["gamma"],
        )
    except ValueError:  # pragma: no cover
        refined_lattice = sample.lattice

    B_opt = refined_lattice.B
    dR_opt = _rotation_from_vector(phi_vec_opt)
    U_opt = dR_opt @ U_fixed
    UB_opt = U_opt @ B_opt

    diff_final = h_phi_obs - (UB_opt @ h_vecs.T).T
    residuals = diff_final
    rms_final = float(np.sqrt(np.mean(diff_final**2)))

    try:
        U_refined = UB_opt @ np.linalg.inv(B_opt)
    except np.linalg.LinAlgError:  # pragma: no cover
        U_refined = None

    sample.lattice = refined_lattice
    sample.UB = UB_opt
    sample.U = U_refined

    return {
        "lattice": refined_lattice,
        "UB": UB_opt,
        "U": U_refined,
        "residuals": residuals,
        "rms": rms_final,
        "converged": converged,
        "n_iter": n_iter,
    }


# ---------------------------------------------------------------------------
# Pure-numpy Nelder-Mead fallback
# ---------------------------------------------------------------------------


def _nelder_mead_numpy(
    f,
    x0: np.ndarray,
    max_iter: int = 1000,
    tol: float = 1e-12,
) -> tuple[np.ndarray, int, bool]:
    """
    Minimal pure-numpy Nelder-Mead simplex minimisation.

    Implements the standard reflection/expansion/contraction/shrink
    operations (Nelder & Mead 1965).
    """
    n = len(x0)
    simplex = np.zeros((n + 1, n))
    simplex[0] = x0
    step = np.where(np.abs(x0) > 1e-8, 0.05 * np.abs(x0), 0.00025)
    for i in range(n):
        simplex[i + 1] = x0.copy()
        simplex[i + 1, i] += step[i]

    fvals = np.array([f(simplex[i]) for i in range(n + 1)])

    alpha, gamma_nm, rho, sigma = 1.0, 2.0, 0.5, 0.5
    converged = False
    n_iter = 0

    for _i in range(1, max_iter + 1):
        n_iter = _i
        order = np.argsort(fvals)
        simplex = simplex[order]
        fvals = fvals[order]

        size = np.max(np.abs(simplex[1:] - simplex[0]))
        if size < tol**0.5 and fvals[0] < tol:
            converged = True
            break

        centroid = np.mean(simplex[:-1], axis=0)
        x_r = centroid + alpha * (centroid - simplex[-1])
        f_r = f(x_r)

        if fvals[0] <= f_r < fvals[-2]:
            simplex[-1] = x_r
            fvals[-1] = f_r
        elif f_r < fvals[0]:
            x_e = centroid + gamma_nm * (x_r - centroid)
            f_e = f(x_e)
            if f_e < f_r:
                simplex[-1] = x_e
                fvals[-1] = f_e
            else:
                simplex[-1] = x_r
                fvals[-1] = f_r
        else:
            x_c = centroid + rho * (simplex[-1] - centroid)
            f_c = f(x_c)
            if f_c < fvals[-1]:
                simplex[-1] = x_c
                fvals[-1] = f_c
            else:
                simplex[1:] = simplex[0] + sigma * (simplex[1:] - simplex[0])
                fvals[1:] = np.array([f(simplex[i]) for i in range(1, n + 1)])

    if not converged:
        logger.debug(
            "Nelder-Mead did not converge after %d iterations (tol=%.2g).",
            n_iter,
            tol,
        )
    return simplex[0], n_iter, converged
