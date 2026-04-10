"""
orientation.py — U and UB matrix computation from orienting reflections.

Functions
---------
angles_to_phi_vector(geometry, **motor_angles)
    Convert a set of motor angles to the scattering vector expressed in the
    phi-axis (innermost sample-stage) frame.  This is the foundational
    computation needed for U and UB matrix determination.

ub_from_one_reflection(sample, reflection, reference_hkl, reference_stage)
    Compute a provisional U and UB from one reflection using the Rodrigues
    rotation that takes the crystal direction ``B @ reference_hkl`` to the
    lab direction given by ``reference_stage``.  Sets ``sample.U`` and
    ``sample.UB`` in-place; returns UB.

ub_identity(sample)
    Set U = I, UB = B; return UB.  The crudest assumption.

Future functions (separate issues):
    ub_from_two_reflections_bl1967   — Busing & Levy 1967, eqs. 23-27  (#5)
    ub_from_three_reflections_bl1967 — Busing & Levy 1967, eqs. 29-31  (#6)

References
----------
Busing & Levy, Acta Cryst. 22, 457-464 (1967)
You, J. Appl. Cryst. 32, 614-623 (1999)
"""

from __future__ import annotations

import numpy as np

from .rotation import rotation_matrix
from .stage import Stage


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

    # Incident-beam direction: longitudinal basis vector (normalised)
    y_hat = np.asarray(geometry.basis["longitudinal"], dtype=float)
    y_norm = np.linalg.norm(y_hat)
    y_hat = y_hat / y_norm

    # Scattering vector in lab frame: Q_lab = (2π/λ) * (D @ ŷ - ŷ)
    two_pi_over_lambda = 2.0 * np.pi / geometry.wavelength
    Q_lab = two_pi_over_lambda * (D @ y_hat - y_hat)

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
       (normalised to ``Bh_hat``).
    2. Extract the lab-frame direction from ``reference_stage`` (normalised
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
    >>> g.sample.reflections.setor1("r1")
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
        for candidate in (
            np.array([1.0, 0.0, 0.0]),
            np.array([0.0, 1.0, 0.0]),
            np.array([0.0, 0.0, 1.0]),
        ):
            ax = np.cross(Bh_hat, candidate)
            if np.linalg.norm(ax) > 1e-10:
                break
        U = rotation_matrix(ax, 180.0)
    else:
        ax = np.cross(Bh_hat, r_hat)
        U = rotation_matrix(ax, np.degrees(angle_rad))

    UB = U @ B
    sample.U = U
    sample.UB = UB
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
