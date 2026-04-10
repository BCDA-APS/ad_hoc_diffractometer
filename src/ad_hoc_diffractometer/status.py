"""
status.py — wh() and pa() status display commands.

These commands print the diffractometer status in a format modelled on
the SPEC ``wh`` (where) and ``pa`` (parameter) commands.

Functions
---------
wh(geometry)
    Print a terse one-screen summary: current reciprocal-space position
    (HKL), wavelength, and current motor angles.  Analogous to SPEC's
    ``wh`` command.

pa(geometry)
    Print a verbose parameter listing: geometry name, orienting
    reflections, lattice constants (real and reciprocal space),
    wavelength.  Analogous to SPEC's ``pa`` command.

References
----------
Align4Pete.log — real SPEC session (7-ID-C fourc, Dec 2020)
"""

from __future__ import annotations

import math


def wh(geometry) -> str:
    """
    Return a terse status string showing the current diffractometer position.

    Analogous to the SPEC ``wh`` (where) command.  Shows the current
    reciprocal-space position (HKL) computed from the live motor angles
    and the active sample's UB matrix, the wavelength, and a table of
    current motor positions.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer to query.  ``geometry.wavelength`` and
        ``geometry.sample.UB`` must both be set for the HKL calculation;
        if either is absent the HKL line shows ``"not available"``.

    Returns
    -------
    str
        Multi-line status string, ready to print.

    Notes
    -----
    The function does not modify any state; it only reads the current
    motor angles via ``stage.angle`` for each stage in the geometry.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> from ad_hoc_diffractometer.status import wh
    >>> g = ahd.fourcv()
    >>> g.wavelength = 1.5406
    >>> print(wh(g))
    H K L =  0  0  0
    Lambda = 1.5406
    <BLANKLINE>
      TwoTheta      Theta        Chi        Phi
         0.000      0.000      0.000      0.000

    References
    ----------
    Align4Pete.log — ``wh`` command outputs, 7-ID-C fourc session, Dec 2020.
    """
    lines: list[str] = []

    # --- HKL position -------------------------------------------------------
    hkl_str = "not available"
    try:
        current_angles = {s.name: s.angle for s in geometry._stages.values()}
        hkl = geometry.inverse(current_angles)
        hkl_str = "  {:g}  {:g}  {:g}".format(*[_clean_zero(v) for v in hkl])
    except Exception:
        pass

    lines.append(f"H K L = {hkl_str}")

    # --- Wavelength ---------------------------------------------------------
    lam = geometry.wavelength
    lam_str = f"{lam:g}" if lam is not None else "not set"
    lines.append(f"Lambda = {lam_str}")

    # --- Motor angle table --------------------------------------------------
    lines.append("")
    stage_names = list(geometry._stages.keys())
    # Header: right-align each name in a 10-char field
    header = "".join(f"{_spec_motor_name(n):>10s}" for n in stage_names)
    lines.append(header)
    values = "".join(f"{geometry._stages[n].angle:>10.3f}" for n in stage_names)
    lines.append(values)

    return "\n".join(lines)


def pa(geometry) -> str:
    """
    Return a verbose parameter listing for the diffractometer.

    Analogous to the SPEC ``pa`` (parameter) command.  Shows the geometry
    name, the two designated orienting reflections (if set), the lattice
    constants in real and reciprocal space, and the wavelength.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        The diffractometer to query.

    Returns
    -------
    str
        Multi-line parameter string, ready to print.

    Notes
    -----
    Fields that have no equivalent in the package yet (diffraction mode,
    sector, azimuthal reference, cut points) are shown with placeholder
    values matching the SPEC fourc defaults to aid visual comparison with
    real SPEC output.

    Examples
    --------
    >>> import ad_hoc_diffractometer as ahd
    >>> from ad_hoc_diffractometer.status import pa
    >>> g = ahd.fourcv()
    >>> print(pa(g))                        # doctest: +ELLIPSIS
    Four-Circle Geometry (fourcv)
    <BLANKLINE>
      Primary Reflection: not set
    ...

    References
    ----------
    Align4Pete.log — ``pa`` command outputs, 7-ID-C fourc session, Dec 2020.
    """
    lines: list[str] = []

    # --- Header -------------------------------------------------------------
    lines.append(f"Geometry: {geometry.name}")
    lines.append("")

    # --- Orienting reflections ----------------------------------------------
    sample = geometry.sample
    ors = sample.reflections.orienting_reflections

    def _refl_block(label: str, refl) -> list[str]:
        if refl is None:
            return [f"  {label}: not set"]
        block = []
        lam_r = refl.wavelength if refl.wavelength is not None else geometry.wavelength
        lam_str = f"{lam_r:g}" if lam_r is not None else "not set"
        block.append(f"  {label} (at lambda {lam_str}):")
        # Motor angles — show in SPEC order if possible
        ang_str = "  ".join(f"{v:g}" for v in refl.angles.values())
        ang_keys = "  ".join(refl.angles.keys())
        block.append(f"    {ang_keys} = {ang_str}")
        h, k, l = refl.hkl  # noqa: E741
        block.append(
            f"    H K L = {_clean_zero(h):g}  {_clean_zero(k):g}  {_clean_zero(l):g}"
        )
        return block

    or1 = ors[0] if len(ors) >= 1 else None
    or2 = ors[1] if len(ors) >= 2 else None
    lines.extend(_refl_block("Primary Reflection", or1))
    lines.append("")
    lines.extend(_refl_block("Secondary Reflection", or2))
    lines.append("")

    # --- Lattice constants --------------------------------------------------
    lat = sample.lattice
    a, b, c = lat.a, lat.b, lat.c
    alpha, beta, gamma = lat.alpha, lat.beta, lat.gamma

    lines.append("  Lattice Constants (lengths / angles):")
    lines.append(f"      real space = {a:g} {b:g} {c:g} / {alpha:g} {beta:g} {gamma:g}")

    # Reciprocal space from Lattice.reciprocal_lattice_vectors
    import numpy as np

    rvecs = lat.reciprocal_lattice_vectors  # tuple of 3 arrays (with 2π factor)
    a_star = float(np.linalg.norm(rvecs[0]))
    b_star = float(np.linalg.norm(rvecs[1]))
    c_star = float(np.linalg.norm(rvecs[2]))

    def _angle_between(u, v):
        cos_a = float(
            np.clip(
                np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)),
                -1.0,
                1.0,
            )
        )
        return math.degrees(math.acos(cos_a))

    r1, r2, r3 = (
        np.asarray(rvecs[0]),
        np.asarray(rvecs[1]),
        np.asarray(rvecs[2]),
    )
    alpha_star = _angle_between(r2, r3)
    beta_star = _angle_between(r1, r3)
    gamma_star = _angle_between(r1, r2)

    lines.append(
        f"    reciprocal space = {a_star:.4g} {b_star:.4g} {c_star:.4g}"
        f" / {alpha_star:.4g} {beta_star:.4g} {gamma_star:.4g}"
    )
    lines.append("")

    # --- Wavelength ---------------------------------------------------------
    lam = geometry.wavelength
    lam_str = f"{lam:g}" if lam is not None else "not set"
    lines.append(f"  Lambda = {lam_str}")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _clean_zero(v: float, atol: float = 1e-10) -> float:
    """Return 0.0 if abs(v) < atol, else v.  Avoids -0.0 or 1e-16 noise."""
    return 0.0 if abs(v) < atol else v


def _spec_motor_name(name: str) -> str:
    """
    Map internal stage names to SPEC-style column headers.

    SPEC uses abbreviated names in motor angle tables (e.g. ``tth`` for
    ``two_theta``, ``th`` for ``omega`` / ``theta``).  For unrecognised
    names the internal name is used as-is.
    """
    _MAP = {
        "two_theta": "TwoTheta",
        "omega": "Theta",
        "chi": "Chi",
        "phi": "Phi",
        "mu": "Mu",
        "eta": "Eta",
        "nu": "Nu",
        "delta": "Delta",
    }
    return _MAP.get(name, name)
