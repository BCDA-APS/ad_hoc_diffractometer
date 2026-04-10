"""
spec.py — SPEC #G1 line parser and emitter for fourc geometry.

The SPEC diffractometer control software stores the complete geometry state
in ``#G`` header lines at the start of every scan.  This module handles the
``#G1`` line for the **fourc** geometry only (Busing & Levy convention,
Eulerian four-circle: omega, chi, phi, two_theta).

.. warning::
   The ``#G1`` field layout is **geometry-specific**.  The implementation
   here applies **only** to the ``fourc`` geometry (``fourcv`` / ``fourch``
   in this package).  Other geometries (sixc, psic, kappa, surf, …) use
   different field counts and meanings.

Functions
---------
parse_fourc_g1(line)
    Parse a SPEC ``#G1`` line for the fourc geometry and return a
    ``FourcG1`` named-tuple.

emit_fourc_g1(g1)
    Format a ``FourcG1`` named-tuple as a SPEC ``#G1`` line string.

g1_to_sample(g1, geometry)
    Apply a parsed ``FourcG1`` to a geometry: set the lattice, add the two
    orienting reflections, designate them as or1/or2.

sample_to_g1(geometry)
    Emit the current sample state of a ``fourcv`` / ``fourch`` geometry as
    a ``FourcG1`` (and thus as a ``#G1`` string).

References
----------
SPEC #G1 format documented in:
    ``references/2020-12-13-fourcc-alignment-7-id-c/spec_G_lines.md``
Busing & Levy, Acta Cryst. 22, 457-464 (1967).
"""

from __future__ import annotations

from typing import NamedTuple

from .lattice import Lattice
from .reflection import Reflection

# ---------------------------------------------------------------------------
# Data container
# ---------------------------------------------------------------------------


class FourcG1(NamedTuple):
    """
    Parsed contents of a SPEC ``#G1`` line for the fourc geometry.

    All 34 fields are represented.  The two ``unused_*`` tuples hold the
    zero padding that SPEC writes but does not interpret.

    Parameters
    ----------
    a, b, c : float
        Direct-lattice parameters (Å).
    alpha, beta, gamma : float
        Direct-lattice angles (degrees).
    a_star, b_star, c_star : float
        Reciprocal-lattice parameters (Å⁻¹, with 2π factor).
    alpha_star, beta_star, gamma_star : float
        Reciprocal-lattice angles (degrees).
    or1_h, or1_k, or1_l : float
        Primary orienting reflection Miller indices.
    or2_h, or2_k, or2_l : float
        Secondary orienting reflection Miller indices.
    or1_two_theta, or1_omega, or1_chi, or1_phi : float
        Primary reflection motor angles (degrees).
        Note: SPEC fourc uses ``2θ, θ(=omega), χ, φ`` — **not** ``two_theta``.
    unused1 : tuple of float
        Two unused zeros following the primary angles (indices 22–23).
    or2_two_theta, or2_omega, or2_chi, or2_phi : float
        Secondary reflection motor angles (degrees).
    unused2 : tuple of float
        Two unused zeros following the secondary angles (indices 28–29).
    lambda1 : float
        Wavelength for primary reflection (Å).
    lambda2 : float
        Wavelength for secondary reflection (Å).
    unused3 : tuple of float
        Two unused zeros at the end (indices 32–33).
    """

    # Direct lattice (fields 0-5)
    a: float
    b: float
    c: float
    alpha: float
    beta: float
    gamma: float
    # Reciprocal lattice (fields 6-11)
    a_star: float
    b_star: float
    c_star: float
    alpha_star: float
    beta_star: float
    gamma_star: float
    # Primary orienting reflection hkl (fields 12-14)
    or1_h: float
    or1_k: float
    or1_l: float
    # Secondary orienting reflection hkl (fields 15-17)
    or2_h: float
    or2_k: float
    or2_l: float
    # Primary angles: 2θ, θ(omega), χ, φ (fields 18-21)
    or1_two_theta: float
    or1_omega: float
    or1_chi: float
    or1_phi: float
    # Unused zeros (fields 22-23)
    unused1: tuple
    # Secondary angles: 2θ, θ(omega), χ, φ (fields 24-27)
    or2_two_theta: float
    or2_omega: float
    or2_chi: float
    or2_phi: float
    # Unused zeros (fields 28-29)
    unused2: tuple
    # Wavelengths (fields 30-31)
    lambda1: float
    lambda2: float
    # Unused zeros (fields 32-33)
    unused3: tuple


# ---------------------------------------------------------------------------
# Parse
# ---------------------------------------------------------------------------

_FOURC_G1_NFIELDS = 34


def parse_fourc_g1(line: str) -> FourcG1:
    """
    Parse a SPEC ``#G1`` line for the fourc geometry.

    Accepts the raw line exactly as it appears in a SPEC data file, with or
    without the leading ``#G1`` tag.

    Parameters
    ----------
    line : str
        A SPEC ``#G1`` line, e.g.::

            #G1 4.785 4.785 12.991 90 90 120 1.516 1.516 0.484 90 90 60 \\
                0 0 6  1 0 0  41.94 20.97 90 0  0 0  60 30 0 0  0 0  1.5498 1.5498  0 0

    Returns
    -------
    FourcG1
        Parsed named-tuple with all 34 fields.

    Raises
    ------
    ValueError
        If the line does not contain exactly 34 numeric tokens (after
        stripping the optional ``#G1`` prefix).

    Notes
    -----
    The function is **geometry-agnostic** in the sense that it does not
    validate whether the lattice or wavelength values are physically
    reasonable.  Callers that need physical validation should check the
    returned ``FourcG1`` fields themselves.

    Examples
    --------
    >>> line = (
    ...     "#G1 4.785 4.785 12.991 90 90 120 "
    ...     "1.516237713 1.516237713 0.483656786 90 90 60 "
    ...     "0 0 6  1 0 0  41.94188 20.97 90 0  0 0  60 30 0 0  0 0  "
    ...     "1.549802558 1.549802558  0 0"
    ... )
    >>> g1 = parse_fourc_g1(line)
    >>> g1.a
    4.785
    >>> g1.or1_h, g1.or1_k, g1.or1_l
    (0.0, 0.0, 6.0)
    """
    # Strip optional leading tag
    stripped = line.strip()
    if stripped.startswith("#G1"):
        stripped = stripped[len("#G1") :].strip()

    tokens = stripped.split()
    if len(tokens) != _FOURC_G1_NFIELDS:
        raise ValueError(
            f"SPEC #G1 line for fourc must have {_FOURC_G1_NFIELDS} numeric fields; "
            f"got {len(tokens)}.  "
            "Ensure this is a fourc geometry #G1 line, not another geometry."
        )

    f = [float(t) for t in tokens]

    return FourcG1(
        # Direct lattice
        a=f[0],
        b=f[1],
        c=f[2],
        alpha=f[3],
        beta=f[4],
        gamma=f[5],
        # Reciprocal lattice
        a_star=f[6],
        b_star=f[7],
        c_star=f[8],
        alpha_star=f[9],
        beta_star=f[10],
        gamma_star=f[11],
        # Primary reflection hkl
        or1_h=f[12],
        or1_k=f[13],
        or1_l=f[14],
        # Secondary reflection hkl
        or2_h=f[15],
        or2_k=f[16],
        or2_l=f[17],
        # Primary angles
        or1_two_theta=f[18],
        or1_omega=f[19],
        or1_chi=f[20],
        or1_phi=f[21],
        unused1=(f[22], f[23]),
        # Secondary angles
        or2_two_theta=f[24],
        or2_omega=f[25],
        or2_chi=f[26],
        or2_phi=f[27],
        unused2=(f[28], f[29]),
        # Wavelengths
        lambda1=f[30],
        lambda2=f[31],
        unused3=(f[32], f[33]),
    )


# ---------------------------------------------------------------------------
# Emit
# ---------------------------------------------------------------------------


def emit_fourc_g1(g1: FourcG1) -> str:
    """
    Format a ``FourcG1`` named-tuple as a SPEC ``#G1`` line string.

    The output is a single-line string starting with ``#G1`` followed by
    34 space-separated floating-point values in the SPEC fourc field order.

    Parameters
    ----------
    g1 : FourcG1
        Parsed or programmatically constructed ``FourcG1`` container.

    Returns
    -------
    str
        A SPEC ``#G1`` line, e.g.::

            #G1 4.785 4.785 12.991 90.0 90.0 120.0 ...

    Examples
    --------
    >>> g1 = parse_fourc_g1(
    ...     "#G1 4.785 4.785 12.991 90 90 120 "
    ...     "1.516237713 1.516237713 0.483656786 90 90 60 "
    ...     "0 0 6  1 0 0  41.94188 20.97 90 0  0 0  60 30 0 0  0 0  "
    ...     "1.549802558 1.549802558  0 0"
    ... )
    >>> line = emit_fourc_g1(g1)
    >>> line.startswith("#G1")
    True
    >>> len(line.split()) - 1  # 34 fields after the tag
    34
    """

    def _fmt(v: float) -> str:
        # Use repr-like formatting: no trailing zeros but enough precision
        return f"{v:g}" if v == int(v) and abs(v) < 1e15 else f"{v:.10g}"

    fields: list[float] = [
        g1.a,
        g1.b,
        g1.c,
        g1.alpha,
        g1.beta,
        g1.gamma,
        g1.a_star,
        g1.b_star,
        g1.c_star,
        g1.alpha_star,
        g1.beta_star,
        g1.gamma_star,
        g1.or1_h,
        g1.or1_k,
        g1.or1_l,
        g1.or2_h,
        g1.or2_k,
        g1.or2_l,
        g1.or1_two_theta,
        g1.or1_omega,
        g1.or1_chi,
        g1.or1_phi,
        g1.unused1[0],
        g1.unused1[1],
        g1.or2_two_theta,
        g1.or2_omega,
        g1.or2_chi,
        g1.or2_phi,
        g1.unused2[0],
        g1.unused2[1],
        g1.lambda1,
        g1.lambda2,
        g1.unused3[0],
        g1.unused3[1],
    ]
    return "#G1 " + " ".join(_fmt(v) for v in fields)


# ---------------------------------------------------------------------------
# Conversion helpers
# ---------------------------------------------------------------------------


def g1_to_sample(g1: FourcG1, geometry) -> None:
    """
    Apply a parsed ``FourcG1`` to a ``fourcv`` / ``fourch`` geometry.

    Sets the active sample's lattice from the direct-lattice parameters in
    ``g1``, adds the two orienting reflections (using fourc motor-angle names
    ``omega``, ``chi``, ``phi``, ``two_theta``), designates them as or1 and
    or2, and sets the geometry wavelength to ``g1.lambda1``.

    Parameters
    ----------
    g1 : FourcG1
        Parsed #G1 data.
    geometry : AdHocDiffractometer
        A ``fourcv`` or ``fourch`` geometry instance.  Must already have an
        active sample (the default ``"test"`` sample is fine).

    Notes
    -----
    - If a reflection named ``"or1"`` or ``"or2"`` already exists on the
      sample, it is removed and replaced.
    - The wavelength is set from ``g1.lambda1``; ``g1.lambda2`` is stored on
      the secondary reflection object but is not set on the geometry (SPEC
      always uses the same wavelength for both in fourc).
    - The fourc motor-angle keys are ``omega``, ``chi``, ``phi``,
      ``two_theta``.  SPEC's ``#G1`` stores ``2θ, θ, χ, φ`` for each
      reflection; ``θ`` is mapped to ``omega``.

    Examples
    --------
    >>> from ad_hoc_diffractometer import fourcv, Lattice
    >>> from ad_hoc_diffractometer.spec import parse_fourc_g1, g1_to_sample
    >>> g = fourcv()
    >>> line = (
    ...     "#G1 4.785 4.785 12.991 90 90 120 "
    ...     "1.516237713 1.516237713 0.483656786 90 90 60 "
    ...     "0 0 6  1 0 0  41.94188 20.97 90 0  0 0  60 30 0 0  0 0  "
    ...     "1.549802558 1.549802558  0 0"
    ... )
    >>> g1 = parse_fourc_g1(line)
    >>> g1_to_sample(g1, g)
    >>> g.wavelength
    1.549802558
    >>> g.sample.lattice.a
    4.785
    """
    sample = geometry.sample

    # Set geometry wavelength from primary wavelength
    geometry.wavelength = g1.lambda1

    # Update sample lattice
    sample.lattice = Lattice(
        a=g1.a,
        b=g1.b,
        c=g1.c,
        alpha=g1.alpha,
        beta=g1.beta,
        gamma=g1.gamma,
    )

    # Remove pre-existing or1/or2 reflections if they exist
    for name in ("or1", "or2"):
        if name in sample.reflections:
            sample.reflections.remove(name)

    # fourc angle mapping: SPEC stores (2θ, θ, χ, φ); θ maps to omega
    or1_angles = {
        "two_theta": g1.or1_two_theta,
        "omega": g1.or1_omega,
        "chi": g1.or1_chi,
        "phi": g1.or1_phi,
    }
    or2_angles = {
        "two_theta": g1.or2_two_theta,
        "omega": g1.or2_omega,
        "chi": g1.or2_chi,
        "phi": g1.or2_phi,
    }

    sample.reflections.add(
        "or1",
        hkl=(g1.or1_h, g1.or1_k, g1.or1_l),
        angles=or1_angles,
        wavelength=g1.lambda1,
    )
    sample.reflections.add(
        "or2",
        hkl=(g1.or2_h, g1.or2_k, g1.or2_l),
        angles=or2_angles,
        wavelength=g1.lambda2,
    )
    sample.reflections.setor1("or1")
    sample.reflections.setor2("or2")


def sample_to_g1(geometry) -> FourcG1:
    """
    Build a ``FourcG1`` from the active sample state of a fourc geometry.

    Reads the active sample's lattice, the two designated orienting
    reflections (or1 and or2), and the geometry wavelength to produce a
    ``FourcG1`` suitable for writing as a SPEC ``#G1`` line.

    Parameters
    ----------
    geometry : AdHocDiffractometer
        A ``fourcv`` or ``fourch`` geometry instance with:

        - ``geometry.wavelength`` set.
        - ``geometry.sample.lattice`` set.
        - ``geometry.sample.reflections.orienting_reflections`` containing
          at least one reflection (two for the full #G1 round-trip).

    Returns
    -------
    FourcG1
        Named-tuple ready for ``emit_fourc_g1()``.

    Raises
    ------
    ValueError
        If ``geometry.wavelength`` is None.
    ValueError
        If no orienting reflections are designated (or1 is required).

    Notes
    -----
    - Reciprocal-lattice parameters are derived from the lattice ``B`` matrix
      using the package's standard computation.
    - If only one orienting reflection is designated (or2 not set), the
      secondary reflection fields are filled with zeros.
    - The fourc motor-angle keys ``omega``, ``chi``, ``phi``, ``two_theta``
      are mapped back to SPEC's ``θ, χ, φ, 2θ`` field order.

    Examples
    --------
    >>> from ad_hoc_diffractometer import fourcv
    >>> from ad_hoc_diffractometer.spec import parse_fourc_g1, g1_to_sample, sample_to_g1
    >>> g = fourcv()
    >>> line = (
    ...     "#G1 4.785 4.785 12.991 90 90 120 "
    ...     "1.516237713 1.516237713 0.483656786 90 90 60 "
    ...     "0 0 6  1 0 0  41.94188 20.97 90 0  0 0  60 30 0 0  0 0  "
    ...     "1.549802558 1.549802558  0 0"
    ... )
    >>> g1_parsed = parse_fourc_g1(line)
    >>> g1_to_sample(g1_parsed, g)
    >>> g1_emitted = sample_to_g1(g)
    >>> g1_emitted.a == g1_parsed.a
    True
    """
    import math

    import numpy as np

    if geometry.wavelength is None:
        raise ValueError("sample_to_g1 requires geometry.wavelength to be set.")

    ors = geometry.sample.reflections.orienting_reflections
    if len(ors) < 1:
        raise ValueError(
            "sample_to_g1 requires at least one designated orienting "
            "reflection (call sample.reflections.setor1() first)."
        )

    sample = geometry.sample
    lat = sample.lattice

    # Direct lattice parameters
    a, b, c = lat.a, lat.b, lat.c
    alpha, beta, gamma = lat.alpha, lat.beta, lat.gamma

    # Reciprocal lattice parameters from the Lattice.reciprocal_lattice_vectors
    # property, which returns the three reciprocal basis vectors as a tuple of
    # numpy arrays (already including the 2π factor per the I16 convention).
    rvecs = lat.reciprocal_lattice_vectors  # tuple of 3 arrays
    r1, r2, r3 = np.asarray(rvecs[0]), np.asarray(rvecs[1]), np.asarray(rvecs[2])

    a_star = float(np.linalg.norm(r1))
    b_star = float(np.linalg.norm(r2))
    c_star = float(np.linalg.norm(r3))

    def _angle_between(u, v):
        cos_a = float(
            np.clip(np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v)), -1.0, 1.0)
        )
        return math.degrees(math.acos(cos_a))

    alpha_star = _angle_between(r2, r3)
    beta_star = _angle_between(r1, r3)
    gamma_star = _angle_between(r1, r2)

    # Primary orienting reflection
    or1: Reflection = ors[0]
    or1_h, or1_k, or1_l = or1.hkl
    or1_two_theta = or1.angles.get("two_theta", 0.0)
    or1_omega = or1.angles.get("omega", 0.0)
    or1_chi = or1.angles.get("chi", 0.0)
    or1_phi = or1.angles.get("phi", 0.0)
    lam1 = or1.wavelength if or1.wavelength is not None else geometry.wavelength

    # Secondary orienting reflection (zeros if not set)
    if len(ors) >= 2:
        or2: Reflection = ors[1]
        or2_h, or2_k, or2_l = or2.hkl
        or2_two_theta = or2.angles.get("two_theta", 0.0)
        or2_omega = or2.angles.get("omega", 0.0)
        or2_chi = or2.angles.get("chi", 0.0)
        or2_phi = or2.angles.get("phi", 0.0)
        lam2 = or2.wavelength if or2.wavelength is not None else geometry.wavelength
    else:
        or2_h = or2_k = or2_l = 0.0
        or2_two_theta = or2_omega = or2_chi = or2_phi = 0.0
        lam2 = lam1

    return FourcG1(
        a=a,
        b=b,
        c=c,
        alpha=alpha,
        beta=beta,
        gamma=gamma,
        a_star=a_star,
        b_star=b_star,
        c_star=c_star,
        alpha_star=alpha_star,
        beta_star=beta_star,
        gamma_star=gamma_star,
        or1_h=or1_h,
        or1_k=or1_k,
        or1_l=or1_l,
        or2_h=or2_h,
        or2_k=or2_k,
        or2_l=or2_l,
        or1_two_theta=or1_two_theta,
        or1_omega=or1_omega,
        or1_chi=or1_chi,
        or1_phi=or1_phi,
        unused1=(0.0, 0.0),
        or2_two_theta=or2_two_theta,
        or2_omega=or2_omega,
        or2_chi=or2_chi,
        or2_phi=or2_phi,
        unused2=(0.0, 0.0),
        lambda1=lam1,
        lambda2=lam2,
        unused3=(0.0, 0.0),
    )
