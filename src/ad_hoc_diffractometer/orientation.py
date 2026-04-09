"""
orientation.py — U and UB matrix computation from orienting reflections.

Functions
---------
ub_identity(sample)
    Set U = I, UB = B; return UB.  The crudest assumption.

ub_from_one_reflection(sample, reflection, reference_hkl, reference_stage)
    Compute a provisional U and UB from one reflection using the Rodrigues
    rotation that takes the crystal direction ``B @ reference_hkl`` to the
    lab direction given by ``reference_stage``.  Sets ``sample.U`` and
    ``sample.UB`` in-place; returns UB.

Future functions (separate issues):
    ub_from_two_reflections_bl1967   — Busing & Levy 1967, eqs. 23-27  (#5)
    ub_from_three_reflections_bl1967 — Busing & Levy 1967, eqs. 29-31  (#6)

References
----------
Busing & Levy, Acta Cryst. 22, 457-464 (1967)
"""

from __future__ import annotations

import numpy as np

from .rotation import rotation_matrix
from .stage import Stage


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
