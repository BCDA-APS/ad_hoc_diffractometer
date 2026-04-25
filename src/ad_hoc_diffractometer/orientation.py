# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
r"""
orientation.py — U and UB matrix computation from orienting reflections.

Functions
---------
angles_to_phi_vector(geometry, \*\*motor_angles)
    Convert a set of motor angles to the scattering vector expressed in the
    phi-axis (innermost sample-stage) frame.  This is the foundational
    computation needed for U and UB matrix determination.

ub_from_one_reflection(sample, reflection, reference_hkl, reference_stage)
    Compute a provisional U and UB from one reflection using the Rodrigues
    rotation that takes the crystal direction ``B @ reference_hkl`` to the
    lab direction given by ``reference_stage``.  Sets ``sample.U`` and
    ``sample.UB`` in-place; returns UB.

ub_from_two_reflections_bl1967(sample, r1, r2)
    Compute U and UB from two orienting reflections using the Busing & Levy
    (1967) algorithm (eqs. 23-27).  Sets ``sample.U`` and ``sample.UB``
    in-place; returns UB.

ub_from_three_reflections_bl1967(sample, r1, r2, r3)
    Compute UB directly from three reflections using the Busing & Levy (1967)
    direct method (eqs. 29-31), without prior knowledge of the lattice.
    Sets ``sample.UB`` in-place; also sets ``sample.U`` if a lattice B is
    available.  Returns UB.

ub_identity(sample)
    Set U = I, UB = B; return UB.  The crudest assumption.

Future functions (separate issues):
    ub_from_three_reflections_bl1967 — Busing & Levy 1967, eqs. 29-31  (#6)

References
----------
* Busing & Levy, Acta Cryst. 22, 457-464 (1967)
* You, J. Appl. Cryst. 32, 614-623 (1999)
"""

from __future__ import annotations

import logging

import numpy as np

from .rotation import rotation_matrix
from .stage import Stage

logger = logging.getLogger(__name__)


def angles_to_phi_vector(geometry, **motor_angles: float) -> np.ndarray:
    """
    Convert a set of motor angles to the scattering vector in the phi frame.

    The "phi frame" is the coordinate system seen from the innermost sample
    stage — the frame in which crystal reflections are expressed when
    computing the orientation (U) matrix.

    Algorithm (Busing & Levy 1967, section "The phi-axis frame"):

    1. Temporarily set the supplied motor angles on their stages (preserving
       the original values so the geometry is restored afterwards).
    2. Compute the total sample rotation matrix ``Z`` (product of all sample
       stage rotation matrices, floor-most first).
    3. Compute the total detector rotation matrix ``D``.
    4. The incident-beam unit vector in the lab frame is ``ŷ`` (longitudinal
       direction, ``geometry.basis["longitudinal"]``).
    5. The scattered-beam unit vector in the lab frame is ``D @ ŷ``.
    6. The scattering vector in the lab frame is::

           Q_lab = (2π / λ) * (D @ ŷ - ŷ)

    7. Rotate Q_lab back through the sample stack::

           Q_phi = Z⁻¹ @ Q_lab = Zᵀ @ Q_lab

       (Z is orthogonal, so Z⁻¹ = Zᵀ.)

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer geometry.  Must have ``wavelength`` set (not None).
    **motor_angles : float
        Motor angles in degrees, keyed by stage name.  All stages present in
        the geometry may be supplied; stages not supplied keep their current
        ``angle`` attribute.  Only sample and detector stages affect the
        result; other stages (if any) are ignored.

    Returns
    -------
    Q_phi : numpy.ndarray, shape (3,)
        Scattering vector in the phi frame, in units of Å⁻¹.

    Raises
    ------
    KeyError
        If a supplied stage name does not exist in the geometry.
    ValueError
        If ``geometry.wavelength`` is None.
    ValueError
        If the geometry has no sample stages.
    ValueError
        If the geometry has no detector stages.

    Notes
    -----
    The function modifies stage angles temporarily and restores them
    afterwards, even if an exception is raised.  It is therefore safe to
    call inside a ``try`` block or from multiple threads as long as each
    call uses a separate geometry instance.

    The scattering vector Q_phi is independent of which sample stage is
    designated the "phi" axis; it is expressed in the frame of the *last*
    sample stage in the stacking order (the one closest to the sample).

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.psic()
    >>> g.wavelength = 1.5406          # Cu Kα in Å
    >>> Q_phi = ahd.angles_to_phi_vector(
    ...     g,
    ...     mu=0, eta=20.97, chi=90, phi=0, nu=0, delta=41.94,
    ... )
    >>> Q_phi  # scattering vector for sapphire (006) in phi frame
    array([...])

    References
    ----------
    Busing & Levy, Acta Cryst. 22, 457-464 (1967) — phi-axis frame
    You, J. Appl. Cryst. 32, 614-623 (1999) — psic geometry conventions
    """
    if geometry.wavelength is None:
        raise ValueError(
            "geometry.wavelength must be set before calling angles_to_phi_vector. "
            "Set it with e.g. geometry.wavelength = 1.5406."
        )
    if not geometry.sample_stages:
        raise ValueError(
            f"Geometry {geometry.name!r} has no sample stages; "
            "cannot compute a phi-frame vector."
        )
    if not geometry.detector_stages:
        raise ValueError(
            f"Geometry {geometry.name!r} has no detector stages; "
            "cannot compute a phi-frame vector."
        )

    # Validate all supplied stage names up-front (raises KeyError for unknown)
    for name in motor_angles:
        geometry.stage(name)  # raises KeyError if not found

    # Save current angles and apply the requested ones
    saved: dict[str, float] = {}
    try:
        for name, angle in motor_angles.items():
            saved[name] = geometry.stage(name).angle
            geometry.set_angle(name, float(angle))

        # Sample rotation matrix Z and detector rotation matrix D
        Z = geometry.sample_rotation_matrix()
        D = geometry.detector_rotation_matrix()

    finally:
        # Restore original angles unconditionally
        for name, angle in saved.items():
            geometry.set_angle(name, angle)

    # Incident-beam direction: longitudinal basis vector (normalized)
    y_hat = np.asarray(geometry.basis["longitudinal"], dtype=float)
    y_norm = np.linalg.norm(y_hat)
    y_hat = y_hat / y_norm

    # Apply diffractometer inclination: when the instrument is mounted at a
    # non-zero angle relative to the beam, the effective beam direction in the
    # diffractometer frame is R_inc.T @ ŷ.  For a zero inclination (R_inc = I)
    # this reduces to the standard y_hat.
    R_inc = geometry.inclination_matrix
    y_eff = R_inc.T @ y_hat

    # Scattering vector in lab frame: Q_lab = (2π/λ) * (D @ ŷ_eff - ŷ_eff)
    two_pi_over_lambda = 2.0 * np.pi / geometry.wavelength
    Q_lab = two_pi_over_lambda * (D @ y_eff - y_eff)

    # Rotate into phi frame: Q_phi = Z^T @ Q_lab  (Z is orthogonal)
    Q_phi = Z.T @ Q_lab

    return Q_phi


def ub_from_one_reflection(
    sample,
    reflection,
    reference_hkl: tuple[float, float, float] = (0.0, 0.0, 1.0),
    reference_stage=None,
) -> np.ndarray:
    """
    Compute a provisional U and UB from one reflection (Rodrigues method).

    A common first step in a diffractometer alignment session: assume that
    a known high-symmetry crystal direction (``reference_hkl``) is parallel
    to a specific diffractometer axis (``reference_stage``).  This gives a
    "fake" UB sufficient to predict where to scan for a second reflection.

    The algorithm:

    1. Compute the crystal-frame direction: ``Bh = B @ reference_hkl``
       (normalized to ``Bh_hat``).
    2. Extract the lab-frame direction from ``reference_stage`` (normalized
       to ``r_hat``).
    3. Find the minimal rotation (Rodrigues) that takes ``Bh_hat`` to
       ``r_hat``:
       ``axis  = cross(Bh_hat, r_hat)``
       ``angle = arccos(clip(dot(Bh_hat, r_hat), -1, 1))``
       ``U     = rotation_matrix(axis, degrees(angle))``
    4. ``UB = U @ B``
    5. Store ``sample.U = U`` and ``sample.UB = UB``; return ``UB``.

    Edge cases:

    - **Parallel** (``angle ≈ 0``): ``Bh_hat`` already points along
      ``r_hat``; ``U = I``.
    - **Anti-parallel** (``angle ≈ π``): choose an arbitrary perpendicular
      axis (the first vector from ``[XHAT, YHAT, ZHAT]`` not parallel to
      ``Bh_hat``); rotate 180° about it.

    .. note::
       The result is approximate.  If ``reference_hkl`` is not truly
       parallel to ``reference_stage.axis`` (e.g. χ = 89.32° rather than
       90.00°), predicted angles for subsequent reflections will be slightly
       wrong.  Refine with ``ub_from_two_reflections_bl1967()`` once a
       second reflection is measured.

    Parameters
    ----------
    sample : Sample
        The sample whose ``U`` and ``UB`` attributes are updated in-place.
        ``sample.lattice`` must be set.  If ``sample.parent`` is set, it
        is used to resolve a string ``reference_stage``.
    reflection : Reflection or str
        A ``Reflection`` object or the name of a reflection in
        ``sample.reflections``.
    reference_hkl : tuple of float, optional
        Miller indices of the crystal direction assumed to be aligned with
        ``reference_stage``.  Default ``(0, 0, 1)`` (c-axis).
    reference_stage : Stage, str, or None, optional
        The diffractometer axis assumed to be parallel to ``reference_hkl``.

        - ``Stage`` object: ``stage.axis`` is used directly (recommended;
          the sign convention is already encoded in the Stage).
        - ``str``: looked up as ``sample.parent.stage(name)``.
        - ``None`` and ``sample.parent`` is set: defaults to
          ``sample.parent.stage("phi")``.
        - ``None`` and ``sample.parent`` is ``None``: raises ``ValueError``.

    Returns
    -------
    UB : numpy.ndarray, shape (3, 3)
        Sets ``sample.U`` and ``sample.UB`` in-place before returning.

    Raises
    ------
    KeyError
        If ``reflection`` is a string not found in ``sample.reflections``.
    ValueError
        If ``reference_stage`` cannot be resolved (no parent, no stage).
    ValueError
        If ``reference_hkl`` maps to the zero vector under B.

    Examples
    --------
    >>> g = psic()
    >>> g.add_sample("sapphire", Lattice(a=4.758, c=12.991))
    >>> g.sample = "sapphire"
    >>> g.add_reflection("r1", hkl=(0, 0, 6),
    ...                  angles={"mu": 0, "eta": 20.97, "chi": 90,
    ...                          "phi": 0, "nu": 0, "delta": 41.94})
    >>> g.sample.reflections.setor0("r1")
    >>> UB = ub_from_one_reflection(
    ...     g.sample, "r1",
    ...     reference_hkl=(0, 0, 1),
    ...     reference_stage=g.stage("phi"),
    ... )
    """
    from .reflection import Reflection

    # --- Resolve reflection ------------------------------------------------
    if isinstance(reflection, str):
        reflection = sample.reflections[reflection]
    if not isinstance(reflection, Reflection):
        raise TypeError(
            f"reflection must be a Reflection or a name string; "
            f"got {type(reflection).__name__!r}."
        )

    # --- Resolve reference_stage → lab-frame axis vector ------------------
    if reference_stage is None:
        if sample.parent is None:
            raise ValueError(
                "reference_stage is required when sample has no parent geometry. "
                "Pass a Stage object or stage name, e.g. geometry.stage('phi')."
            )
        reference_stage = sample.parent.stage("phi")

    if isinstance(reference_stage, str):
        if sample.parent is None:
            raise ValueError(
                f"Cannot look up stage {reference_stage!r}: "
                f"sample has no parent geometry."
            )
        reference_stage = sample.parent.stage(reference_stage)

    if isinstance(reference_stage, Stage):
        r_vec = np.asarray(reference_stage.axis, dtype=float)
    else:
        r_vec = np.asarray(reference_stage, dtype=float)

    r_norm = np.linalg.norm(r_vec)
    if r_norm < 1e-14:
        raise ValueError("reference_stage axis vector is zero.")
    r_hat = r_vec / r_norm

    # --- Crystal-frame direction from reference_hkl and B ------------------
    B = sample.lattice.B
    Bh = B @ np.asarray(reference_hkl, dtype=float)
    Bh_norm = np.linalg.norm(Bh)
    if Bh_norm < 1e-14:
        raise ValueError(
            f"reference_hkl {reference_hkl} maps to the zero vector under B. "
            f"Choose a non-zero Miller index triple."
        )
    Bh_hat = Bh / Bh_norm

    # --- Rodrigues rotation from Bh_hat to r_hat ---------------------------
    cos_angle = float(np.clip(np.dot(Bh_hat, r_hat), -1.0, 1.0))
    angle_rad = np.arccos(cos_angle)

    if abs(angle_rad) < 1e-10:
        # Already parallel — U = I
        U = np.eye(3)
    elif abs(angle_rad - np.pi) < 1e-10:
        # Anti-parallel — choose any perpendicular axis
        for candidate in (  # pragma: no branch
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ):
            ax = np.cross(Bh_hat, candidate)
            if np.linalg.norm(ax) > 1e-10:
                break
        # At least one standard basis vector is guaranteed non-parallel to
        # any unit vector, so the loop always breaks.
        U = rotation_matrix(ax, 180.0)
    else:
        ax = np.cross(Bh_hat, r_hat)
        U = rotation_matrix(ax, np.degrees(angle_rad))

    UB = U @ B
    sample.U = U
    sample.UB = UB
    return UB


def _gram_schmidt_triple(v1: np.ndarray, v2: np.ndarray) -> np.ndarray:
    """
    Build a right-handed orthonormal 3×3 matrix from two linearly independent
    vectors using Gram-Schmidt orthogonalisation.

    The columns of the returned matrix ``T`` are::

        t1 = v1 / |v1|
        t3 = t1 × v2 / |t1 × v2|
        t2 = t3 × t1

    so that ``t1 ∥ v1``, ``t2`` lies in the plane of ``v1`` and ``v2``, and
    ``t3`` is perpendicular to that plane.  The triple ``(t1, t2, t3)`` is
    right-handed and orthonormal.

    Parameters
    ----------
    v1 : numpy.ndarray, shape (3,)
        Primary vector (must be non-zero).
    v2 : numpy.ndarray, shape (3,)
        Secondary vector (must not be parallel to ``v1``).

    Returns
    -------
    T : numpy.ndarray, shape (3, 3)
        Columns are ``[t1, t2, t3]``, forming a right-handed orthonormal basis.

    Raises
    ------
    ValueError
        If ``v1`` is the zero vector or ``v1`` and ``v2`` are parallel
        (cross product is zero).

    Notes
    -----
    This is the Gram-Schmidt construction used in Busing & Levy (1967) to
    build the orthonormal triples ``Tc`` (crystal frame) and ``Tφ`` (phi
    frame) for the two-reflection orientation algorithm (eqs. 23-27).
    """
    v1 = np.asarray(v1, dtype=float)
    v2 = np.asarray(v2, dtype=float)

    n1 = np.linalg.norm(v1)
    if n1 < 1e-14:
        raise ValueError(
            "_gram_schmidt_triple: v1 is the zero vector; cannot build an "
            "orthonormal basis."
        )
    t1 = v1 / n1

    cross = np.cross(t1, v2)
    n_cross = np.linalg.norm(cross)
    if n_cross < 1e-14:
        raise ValueError(
            "_gram_schmidt_triple: v1 and v2 are parallel (or v2 is zero); "
            "cannot build an orthonormal basis."
        )
    t3 = cross / n_cross
    t2 = np.cross(t3, t1)

    return np.column_stack([t1, t2, t3])


def ub_from_two_reflections_bl1967(
    sample,
    r1=None,
    r2=None,
) -> np.ndarray:
    """
    Compute U and UB from two orienting reflections (Busing & Levy 1967, eqs. 23-27).

    Given two reflections with known hkl and measured motor angles, and a
    known lattice (B matrix), this function computes the orientation matrix U
    and then UB = U @ B, storing both on the sample.

    Algorithm (BL1967 eqs. 23-27):

    1. For each reflection, call ``angles_to_phi_vector()`` to get the
       scattering vector in the phi frame: ``u1φ``, ``u2φ``.
    2. From hkl and the lattice B matrix: ``h1c = B @ h1``, ``h2c = B @ h2``.
    3. Build orthonormal triple ``Tc`` in the crystal frame via Gram-Schmidt:
       ``t1c ∥ h1c``, ``t2c`` in the plane of ``h1c`` and ``h2c``,
       ``t3c = t1c × t2c``.
    4. Build the matching triple ``Tφ`` in the phi frame from ``u1φ``, ``u2φ``.
    5. Compute ``U = Tφ @ Tc.T``  (eq. 27; Tc is orthogonal so Tc⁻¹ = Tc.T).
    6. Compute ``UB = U @ B``.
    7. Store ``sample.U = U``, ``sample.UB = UB``; return ``UB``.

    Parameters
    ----------
    sample : Sample
        The sample whose ``U`` and ``UB`` attributes are updated in-place.
        ``sample.lattice`` must be set.  ``sample.parent`` must be a geometry
        with ``wavelength`` set (it is used to call ``angles_to_phi_vector``).
    r1 : Reflection, str, or None
        Primary orienting reflection.  If ``None``, defaults to
        ``sample.reflections.orienting_reflections[0]``.
    r2 : Reflection, str, or None
        Secondary orienting reflection.  If ``None``, defaults to
        ``sample.reflections.orienting_reflections[1]``.

    Returns
    -------
    UB : numpy.ndarray, shape (3, 3)
        Sets ``sample.U`` (first) and ``sample.UB = sample.U @ B`` in-place
        before returning.

    Raises
    ------
    KeyError
        If ``r1`` or ``r2`` is a string not found in ``sample.reflections``.
    ValueError
        If ``r1`` or ``r2`` is ``None`` and the required orienting reflection
        has not been designated (``setor0``/``setor1`` not called).
    ValueError
        If ``sample.parent`` is ``None`` (needed to call
        ``angles_to_phi_vector``).
    ValueError
        If ``sample.parent.wavelength`` is ``None``.
    ValueError
        If the two reflections are parallel in the crystal frame (h1c and
        h2c collinear) or in the phi frame (u1φ and u2φ collinear).
    TypeError
        If ``r1`` or ``r2`` is not a ``Reflection``, string, or ``None``.

    Notes
    -----
    U is computed first (``sample.U = Tφ @ Tc.T``), then UB is derived from
    it (``sample.UB = sample.U @ B``).

    The wavelength used for ``angles_to_phi_vector`` is taken from
    ``sample.parent.wavelength``.  If a reflection carries its own
    ``wavelength`` attribute, that is *not* used here; the geometry's
    wavelength governs the conversion from motor angles to Q_phi.

    References
    ----------
    Busing & Levy, Acta Cryst. 22, 457-464 (1967), eqs. 23-27.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> g = ahd.psic()
    >>> g.wavelength = 1.5406
    >>> g.add_sample("sapphire", ahd.Lattice(a=4.758, c=12.991))
    >>> g.sample = "sapphire"
    >>> g.add_reflection("r1", hkl=(0, 0, 6),
    ...     angles={"mu": 0, "eta": 20.97, "chi": 90, "phi": 0,
    ...             "nu": 0, "delta": 41.94})
    >>> g.add_reflection("r2", hkl=(1, 0, 4),
    ...     angles={"mu": 0, "eta": 23.72, "chi": 57.04, "phi": 0,
    ...             "nu": 0, "delta": 48.13})
    >>> g.sample.reflections.setor0("r1")
    >>> g.sample.reflections.setor1("r2")
    >>> UB = ahd.ub_from_two_reflections_bl1967(g.sample)
    """
    from .reflection import Reflection

    # --- Require a parent geometry (needed for angles_to_phi_vector) ----------
    if sample.parent is None:
        raise ValueError(
            "ub_from_two_reflections_bl1967 requires sample.parent to be set "
            "(an AdHocDiffractometer with wavelength).  Attach the sample to a "
            "geometry before calling this function."
        )
    geometry = sample.parent

    # --- Resolve r1 -----------------------------------------------------------
    if r1 is None:
        ors = sample.reflections.orienting_reflections
        if len(ors) < 1 or ors[0] is None:
            raise ValueError(
                "r1 is None and no primary orienting reflection (or1) has been "
                "designated.  Call sample.reflections.setor0() first, or pass "
                "r1 explicitly."
            )
        r1 = ors[0]
    elif isinstance(r1, str):
        r1 = sample.reflections[r1]
    if not isinstance(r1, Reflection):
        raise TypeError(
            f"r1 must be a Reflection, a name string, or None; "
            f"got {type(r1).__name__!r}."
        )

    # --- Resolve r2 -----------------------------------------------------------
    if r2 is None:
        ors = sample.reflections.orienting_reflections
        if len(ors) < 2:
            raise ValueError(
                "r2 is None and no secondary orienting reflection (or2) has been "
                "designated.  Call sample.reflections.setor1() first, or pass "
                "r2 explicitly."
            )
        r2 = ors[1]
    elif isinstance(r2, str):
        r2 = sample.reflections[r2]
    if not isinstance(r2, Reflection):
        raise TypeError(
            f"r2 must be a Reflection, a name string, or None; "
            f"got {type(r2).__name__!r}."
        )

    # --- Phi-frame scattering vectors from motor angles -----------------------
    u1_phi = angles_to_phi_vector(geometry, **r1.angles)
    u2_phi = angles_to_phi_vector(geometry, **r2.angles)

    # --- Crystal-frame vectors from hkl and B ---------------------------------
    B = sample.lattice.B
    h1c = B @ np.asarray(r1.hkl, dtype=float)
    h2c = B @ np.asarray(r2.hkl, dtype=float)

    # --- Build orthonormal triples via Gram-Schmidt ---------------------------
    # Tc: crystal frame.  t1c ∥ h1c, t2c in plane(h1c, h2c), t3c = t1c × t2c.
    try:
        Tc = _gram_schmidt_triple(h1c, h2c)
    except ValueError as exc:
        raise ValueError(
            "The two reflections are parallel in the crystal frame "
            f"(h1c = B @ {r1.hkl} and h2c = B @ {r2.hkl} are collinear). "
            "Choose two reflections that are not parallel."
        ) from exc

    # Tphi: phi frame.  t1φ ∥ u1φ, t2φ in plane(u1φ, u2φ), t3φ = t1φ × t2φ.
    try:
        Tphi = _gram_schmidt_triple(u1_phi, u2_phi)
    except ValueError as exc:
        raise ValueError(
            "The two reflections are parallel in the phi frame "
            f"(Q_phi vectors for {r1.name!r} and {r2.name!r} are collinear). "
            "Choose two reflections that are not parallel."
        ) from exc

    # --- BL1967 eq. 27: U = Tφ @ Tc.T  (Tc orthogonal → Tc⁻¹ = Tc.T) ------
    U = Tphi @ Tc.T
    UB = U @ B

    sample.U = U
    sample.UB = UB
    return UB


def ub_from_three_reflections_bl1967(
    sample,
    r1,
    r2,
    r3,
) -> np.ndarray:
    """
    Compute UB directly from three reflections (Busing & Levy 1967, eqs. 29-31).

    This method requires no prior knowledge of the lattice: it computes UB
    directly by matrix inversion from three measured reflections.  If a
    lattice B matrix is available on the sample, U is also derived.

    Algorithm (BL1967 eqs. 28-31):

    1. For each reflection i compute the phi-frame scattering vector
       ``hiφ = angles_to_phi_vector(geometry, **ri.angles)`` (eq. 28 gives
       the magnitude; ``angles_to_phi_vector`` already carries the full
       vector in Å⁻¹).
    2. Stack as column matrices::

           Hφ = [h1φ | h2φ | h3φ]    (3×3, columns are phi-frame vectors)
           H  = [h1  | h2  | h3 ]    (3×3, columns are Miller-index triples)

    3. ``UB = Hφ @ inv(H)``  (eq. 31).
    4. If ``sample.lattice`` is set: ``U = UB @ inv(B)`` (derived from UB).
    5. Store ``sample.UB = UB`` (first) and ``sample.U = U``; return ``UB``.

    Parameters
    ----------
    sample : Sample
        The sample whose ``UB`` (and ``U``) attributes are updated in-place.
        ``sample.parent`` must be a geometry with ``wavelength`` set.
    r1, r2, r3 : Reflection or str
        Three orienting reflections.  Each may be a ``Reflection`` object or
        the name of a reflection in ``sample.reflections``.

    Returns
    -------
    UB : numpy.ndarray, shape (3, 3)
        ``sample.UB`` is set first (directly, via eq. 31), then
        ``sample.U = UB @ inv(B)`` is derived.  Both are set in-place.

    Raises
    ------
    KeyError
        If any of ``r1``, ``r2``, ``r3`` is a string not found in
        ``sample.reflections``.
    TypeError
        If any argument is not a ``Reflection`` or a string.
    ValueError
        If ``sample.parent`` is ``None``.
    ValueError
        If the three Miller-index column matrix ``H`` is singular
        (``|det(H)| < tol``), i.e. the three hkl vectors are coplanar.
    ValueError
        If ``sample.parent.wavelength`` is ``None``.

    Warns
    -----
    UserWarning
        If ``det(H) < 0``, the hkl triples form a left-handed system.
        The computation proceeds but the sign convention may give U with
        ``det(U) = -1``; consider swapping r1 and r2 to make det(H) > 0.

    Notes
    -----
    UB is computed first (``sample.UB = Hφ @ H⁻¹``).  U is then derived
    as ``sample.U = UB @ B⁻¹``.  This is the opposite order from
    ``ub_from_two_reflections_bl1967``, where U is computed first.

    The method does not require a known lattice: ``H`` is formed from the
    raw hkl indices, not from ``B @ hkl``.  However, if ``sample.lattice``
    is the package default (cubic, a = 1 Å) rather than a measured lattice,
    the derived U will not be physically meaningful.

    If ``det(H)`` is exactly zero (degenerate reflections), ``numpy.linalg.inv``
    will raise ``LinAlgError``, which is caught and re-raised as ``ValueError``.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> import math
    >>> g = ahd.psic()
    >>> g.wavelength = 2 * math.pi
    >>> g.sample.lattice = ahd.Lattice(a=2 * math.pi)
    >>> g.add_reflection("r1", hkl=(1, 0, 0),
    ...     angles={"mu": 0, "eta": 30, "chi": 0, "phi": 0, "nu": 0, "delta": 60})
    >>> g.add_reflection("r2", hkl=(0, 1, 0),
    ...     angles={"mu": 0, "eta": 30, "chi": 0, "phi": 90, "nu": 0, "delta": 60})
    >>> g.add_reflection("r3", hkl=(0, 0, 1),
    ...     angles={"mu": 0, "eta": 30, "chi": 90, "phi": 30, "nu": 0, "delta": 60})
    >>> UB = ahd.ub_from_three_reflections_bl1967(g.sample, "r1", "r2", "r3")

    References
    ----------
    Busing & Levy, Acta Cryst. 22, 457-464 (1967), eqs. 28-31.
    """
    import warnings

    from .reflection import Reflection

    # --- Require a parent geometry -------------------------------------------
    if sample.parent is None:
        raise ValueError(
            "ub_from_three_reflections_bl1967 requires sample.parent to be set "
            "(an AdHocDiffractometer with wavelength).  Attach the sample to a "
            "geometry before calling this function."
        )
    geometry = sample.parent

    # --- Resolve reflections -------------------------------------------------
    def _resolve(r, label: str) -> Reflection:
        if isinstance(r, str):
            r = sample.reflections[r]
        if not isinstance(r, Reflection):
            raise TypeError(
                f"{label} must be a Reflection or a name string; "
                f"got {type(r).__name__!r}."
            )
        return r

    r1 = _resolve(r1, "r1")
    r2 = _resolve(r2, "r2")
    r3 = _resolve(r3, "r3")

    # --- Phi-frame scattering vectors ----------------------------------------
    h1_phi = angles_to_phi_vector(geometry, **r1.angles)
    h2_phi = angles_to_phi_vector(geometry, **r2.angles)
    h3_phi = angles_to_phi_vector(geometry, **r3.angles)

    # --- Column matrices Hφ and H -------------------------------------------
    H_phi = np.column_stack([h1_phi, h2_phi, h3_phi])  # 3×3
    H = np.column_stack(
        [
            np.asarray(r1.hkl, dtype=float),
            np.asarray(r2.hkl, dtype=float),
            np.asarray(r3.hkl, dtype=float),
        ]
    )  # 3×3

    # --- Check H is non-singular --------------------------------------------
    det_H = float(np.linalg.det(H))
    tol = 1e-10
    if abs(det_H) < tol:
        raise ValueError(
            "The three reflections are coplanar in reciprocal space: "
            f"det(H) = {det_H:.3g} ≈ 0.  Choose three reflections whose "
            "Miller-index vectors are linearly independent."
        )

    if det_H < 0:
        msg = (
            f"det(H) = {det_H:.6g} < 0: the three hkl triples form a "
            "left-handed system.  U may have det(U) = -1.  Consider "
            "swapping r1 and r2 to restore a right-handed system."
        )
        logger.warning(msg)
        warnings.warn(msg, UserWarning, stacklevel=2)

    # --- BL1967 eq. 31: UB = Hφ @ H⁻¹ -------------------------------------
    try:
        H_inv = np.linalg.inv(H)
    except np.linalg.LinAlgError as exc:
        raise ValueError(
            "Cannot invert the Miller-index matrix H; the three reflections "
            "are linearly dependent (coplanar in reciprocal space)."
        ) from exc

    UB = H_phi @ H_inv

    # --- Derive U = UB @ B⁻¹ -----------------------------------------------
    B = sample.lattice.B
    try:
        B_inv = np.linalg.inv(B)
    except np.linalg.LinAlgError:
        U = None
    else:
        U = UB @ B_inv

    # --- Store in-place (UB first, then U) ----------------------------------
    sample.UB = UB
    sample.U = U
    return UB


def ub_identity(sample) -> np.ndarray:
    """
    Set U = I (identity) and UB = B; return UB.

    The crudest orientation assumption: crystal axes are aligned with the
    lab axes.  Useful as a starting point when no reflection information
    is available at all.

    Parameters
    ----------
    sample : Sample
        The sample whose ``U`` and ``UB`` attributes are updated in-place.

    Returns
    -------
    UB : numpy.ndarray, shape (3, 3)
    """
    B = sample.lattice.B
    sample.U = np.eye(3)
    sample.UB = B.copy()
    return sample.UB
