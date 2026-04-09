"""
ad_hoc_diffractometer — multi-circle diffractometer geometry calculations.

Suggested import alias:  import ad_hoc_diffractometer as ahd

Based on:
  Busing & Levy, Acta Cryst. 22, 457-464 (1967)
  H. You, J. Appl. Cryst. 32, 614-623 (1999). DOI: 10.1107/S0021889899001223
  M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
  D.A. Walko, Reference Module in Materials Science and Materials Engineering (2016).
"""

import numpy as np

# ---------------------------------------------------------------------------
# Standard basis vectors in the You (1999) / problem2.md convention:
#   xHat = vertical, yHat = longitudinal (beam), zHat = lateral
# ---------------------------------------------------------------------------

XHAT = np.array([1.0, 0.0, 0.0])
YHAT = np.array([0.0, 1.0, 0.0])
ZHAT = np.array([0.0, 0.0, 1.0])


def rotation_matrix(axis, angle_deg):
    """
    Compute a right-handed rotation matrix about a unit axis vector.

    Uses the Rodrigues formula:
        R = I cos(theta) + (1 - cos(theta))(n x n) + sin(theta)[n]_x

    For a left-handed rotation (e.g. eta, phi, delta in You 1999), pass the
    negated axis: rotation_matrix(-ZHAT, angle_deg).

    Parameters
    ----------
    axis : numpy.ndarray, shape (3,)
        Unit vector defining the rotation axis.  Need not be normalised;
        it will be normalised internally.
    angle_deg : float
        Rotation angle in degrees (right-handed sense).

    Returns
    -------
    R : numpy.ndarray, shape (3, 3)
        Rotation matrix.
    """
    n = np.asarray(axis, dtype=float)
    n = n / np.linalg.norm(n)
    theta = np.deg2rad(angle_deg)
    c = np.cos(theta)
    s = np.sin(theta)
    nx, ny, nz = n
    skew = np.array(
        [
            [0, -nz, ny],
            [nz, 0, -nx],
            [-ny, nx, 0],
        ]
    )
    return c * np.eye(3) + (1 - c) * np.outer(n, n) + s * skew


# ---------------------------------------------------------------------------
# Diffractometer stage and geometry description
# ---------------------------------------------------------------------------


class Stage:
    """
    One rotary stage of a diffractometer.

    Each stage is characterised by:
      - a name (e.g. 'mu', 'eta', 'chi')
      - a rotation axis expressed as a signed Cartesian vector in the lab frame
        (e.g. +XHAT for right-handed vertical, -ZHAT for left-handed lateral)
      - a parent stage name (the stage on which this one sits), or None if it
        sits on the floor / fixed lab frame

    The sign of the axis vector encodes the handedness of the rotation:
      +nHat  =>  right-handed rotation about nHat
      -nHat  =>  left-handed rotation about nHat
                 (equivalent to right-handed rotation about nHat with negated angle)

    Parameters
    ----------
    name : str
        Human-readable name for this stage (e.g. 'mu', 'S2-1').
    axis : numpy.ndarray, shape (3,)
        Signed rotation axis vector in the lab frame.
    parent : str or None
        Name of the stage on which this stage is mounted, or None if it
        sits directly on the lab frame (floor).
    role : str, optional
        'sample' or 'detector', for bookkeeping.  Default is 'sample'.
    angle : float, optional
        Current angle setting in degrees.  Default is 0.0.
    """

    def __init__(self, name, axis, parent=None, role="sample", angle=0.0):
        self.name = name
        self.axis = np.asarray(axis, dtype=float)
        self.parent = parent
        self.role = role
        self.angle = angle

    def rotation_matrix(self):
        """Return the 3x3 rotation matrix for the current angle setting."""
        return rotation_matrix(self.axis, self.angle)

    def __repr__(self):
        return (
            f"Stage(name={self.name!r}, axis={self.axis}, "
            f"parent={self.parent!r}, role={self.role!r}, angle={self.angle})"
        )


class AdHocDiffractometer:
    """
    Description of a diffractometer as an ordered collection of rotary stages.

    The geometry is defined by a list of Stage objects.  The stacking order
    (which stage sits on which) is encoded by the parent attribute of each
    Stage: a stage with parent=None sits on the fixed lab frame; all other
    stages sit on the stage named by their parent.

    The class supports arbitrary geometries including:
      - You (1999) psic  4S+2D  six-circle
      - Lohmeier & Vlieg (1993) sixc  six-circle
      - Busing & Levy (1967) fourc  four-circle Eulerian
      - Kappa four-circle

    The basis vectors for the lab frame are specified at construction time,
    allowing the same physical geometry to be described in different
    coordinate conventions.

    Parameters
    ----------
    name : str
        Name of the diffractometer geometry (e.g. 'psic', 'sixc', 'fourc').
    stages : list of Stage
        All stages in the geometry.  Order within the list does not determine
        the stacking; the parent attribute of each Stage does.
    basis : dict, optional
        Mapping from physical direction names to Cartesian basis vectors.
        Default is the You (1999) convention:
            {'vertical': XHAT, 'longitudinal': YHAT, 'lateral': ZHAT}
    description : str, optional
        Free-text description of the geometry.

    Attributes
    ----------
    sample_stages : list of Stage
        Stages with role='sample', in stacking order (floor first).
    detector_stages : list of Stage
        Stages with role='detector', in stacking order (floor first).
    """

    DEFAULT_BASIS = {
        "vertical": XHAT,
        "longitudinal": YHAT,
        "lateral": ZHAT,
    }

    def __init__(self, name, stages, basis=None, description=""):
        self.name = name
        self.description = description
        self.basis = basis if basis is not None else dict(self.DEFAULT_BASIS)

        # Validate basis vectors
        self._check_basis()

        # Store stages in a dict keyed by name for fast lookup
        self._stages = {s.name: s for s in stages}
        if len(self._stages) != len(stages):
            raise ValueError("Stage names must be unique.")

        # Validate parent references
        for s in stages:
            if s.parent is not None and s.parent not in self._stages:
                raise ValueError(
                    f"Stage {s.name!r} has parent {s.parent!r} "
                    f"which is not in the stage list."
                )

        # Check for cycles in the parent graph
        self._check_no_cycles()

        # Build ordered lists per role
        self.sample_stages = self._ordered_stages("sample")
        self.detector_stages = self._ordered_stages("detector")

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _check_basis(self):
        """
        Validate that the basis contains exactly three mutually orthogonal,
        non-zero, 3-dimensional vectors.

        The basis is a labelled dict; its keys have no defined ordering, so
        right-handedness cannot be checked here (the cross product of two
        vectors depends on which is 'first').  Right-handedness is the
        caller's responsibility and should be verified by the geometry
        factory functions or the user.

        Raises
        ------
        ValueError
            If the basis does not have exactly three entries, if any vector
            is not 3-dimensional or is zero, or if any two vectors are not
            mutually orthogonal (normalised dot product exceeds tolerance).
        """
        vecs = list(self.basis.values())
        names = list(self.basis.keys())

        if len(vecs) != 3:
            raise ValueError(f"Basis must contain exactly 3 vectors; got {len(vecs)}.")

        for name, v in zip(names, vecs, strict=False):
            v = np.asarray(v, dtype=float)
            if v.shape != (3,):
                raise ValueError(
                    f"Basis vector {name!r} must be 3-dimensional; got shape {v.shape}."
                )
            if np.linalg.norm(v) == 0.0:
                raise ValueError(f"Basis vector {name!r} must be non-zero.")

        atol = 1e-10
        for i, (n1, v1) in enumerate(zip(names, vecs, strict=False)):
            for n2, v2 in zip(names[i + 1 :], vecs[i + 1 :], strict=False):
                v1n = np.asarray(v1, dtype=float) / np.linalg.norm(v1)
                v2n = np.asarray(v2, dtype=float) / np.linalg.norm(v2)
                dot = np.dot(v1n, v2n)
                if abs(dot) > atol:
                    raise ValueError(
                        f"Basis vectors {n1!r} and {n2!r} are not orthogonal "
                        f"(normalised dot product = {dot:.6g}, tolerance {atol})."
                    )

    def _check_no_cycles(self):
        """Raise ValueError if the parent graph contains a cycle."""
        for start in self._stages:
            visited = set()
            node = self._stages[start].parent
            while node is not None:
                if node in visited:
                    raise ValueError(
                        f"Cycle detected in parent chain starting from {start!r}."
                    )
                visited.add(node)
                node = self._stages[node].parent

    def _ordered_stages(self, role):
        """
        Return stages of the given role in stacking order (floor-most first).

        Uses a topological sort on the parent graph restricted to stages of
        the requested role.
        """
        role_stages = [s for s in self._stages.values() if s.role == role]

        # Topological sort: repeatedly pick stages whose parent is not in the
        # remaining set (i.e. parent is None or parent belongs to a different role)
        remaining = list(role_stages)
        ordered = []
        role_names = {s.name for s in role_stages}

        max_iter = len(remaining) + 1
        while remaining:
            max_iter -= 1
            if max_iter < 0:
                raise RuntimeError(
                    "Could not determine stacking order; check parent chain."
                )
            for s in remaining:
                if s.parent is None or s.parent not in role_names:
                    ordered.append(s)
                    remaining.remove(s)
                    role_names.discard(s.name)
                    break

        return ordered

    # ------------------------------------------------------------------
    # Public interface
    # ------------------------------------------------------------------

    def stage(self, name):
        """Return the Stage with the given name."""
        return self._stages[name]

    def set_angle(self, name, angle_deg):
        """Set the angle of a named stage in degrees."""
        self._stages[name].angle = angle_deg

    def sample_rotation_matrix(self):
        """
        Compute the total sample rotation matrix Z = R_floor * ... * R_top.

        Stages are applied in stacking order: the floor-most stage is applied
        first (leftmost in the product), the innermost stage last.

        Returns
        -------
        Z : numpy.ndarray, shape (3, 3)
        """
        Z = np.eye(3)
        for s in self.sample_stages:
            Z = s.rotation_matrix() @ Z
        return Z

    def detector_rotation_matrix(self):
        """
        Compute the total detector rotation matrix in stacking order.

        Returns
        -------
        D : numpy.ndarray, shape (3, 3)
        """
        D = np.eye(3)
        for s in self.detector_stages:
            D = s.rotation_matrix() @ D
        return D

    def summary(self):
        """Print a human-readable summary of the geometry."""
        lines = [
            f"Diffractometer geometry: {self.name}",
            f"  {self.description}",
            "",
            "  Basis vectors:",
        ]
        for direction, vec in self.basis.items():
            lines.append(f"    {direction:12s} -> {vec}")
        lines.append("")
        lines.append("  Sample stages (floor first):")
        for s in self.sample_stages:
            lines.append(
                f"    {s.name:8s}  axis={s.axis}  parent={s.parent}  angle={s.angle} deg"
            )
        lines.append("")
        lines.append("  Detector stages (floor first):")
        for s in self.detector_stages:
            lines.append(
                f"    {s.name:8s}  axis={s.axis}  parent={s.parent}  angle={s.angle} deg"
            )
        print("\n".join(lines))

    def __repr__(self):
        return (
            f"AdHocDiffractometer(name={self.name!r}, "
            f"stages={list(self._stages.keys())})"
        )


# ---------------------------------------------------------------------------
# Predefined geometry factories
# ---------------------------------------------------------------------------


def geometry_psic():
    """
    You (1999) '4S+2D' six-circle diffractometer (psic geometry).

    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.

    Sample stack (floor first):
        mu  (S2-1): vertical,     +xHat, right-handed
        eta (S2-2): lateral,      -zHat, left-handed
        chi (S2-3): longitudinal, +yHat, right-handed
        phi (S2-4): lateral,      -zHat, left-handed

    Detector stack (floor first):
        nu    (S1-1): vertical, +xHat, right-handed
        delta (S1-2): lateral,  -zHat, left-handed

    Reference: H. You, J. Appl. Cryst. 32, 614-623 (1999).
               DOI: 10.1107/S0021889899001223
    """
    stages = [
        Stage("mu", +XHAT, parent=None, role="sample"),
        Stage("eta", -ZHAT, parent="mu", role="sample"),
        Stage("chi", +YHAT, parent="eta", role="sample"),
        Stage("phi", -ZHAT, parent="chi", role="sample"),
        Stage("nu", +XHAT, parent=None, role="detector"),
        Stage("delta", -ZHAT, parent="nu", role="detector"),
    ]
    return AdHocDiffractometer(
        name="psic",
        stages=stages,
        description="You (1999) 4S+2D six-circle diffractometer",
    )


def geometry_fourc():
    """
    Busing & Levy (1967) four-circle Eulerian diffractometer.

    Basis: xHat=longitudinal (beam), yHat=lateral, zHat=vertical.
    (Busing & Levy convention: z vertical, y along beam, x = scattering vector at zero.)

    Sample stack (floor first):
        omega (2theta/2): vertical, left-handed about zHat
        chi:              lateral,  right-handed about yHat
        phi:              vertical, left-handed about zHat

    Detector:
        two_theta: vertical, left-handed about zHat

    Reference: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
    """
    # Busing & Levy right-handed basis: x=lateral, y=longitudinal (beam), z=vertical
    # Right-handed: lateral x longitudinal = vertical => xHat x yHat = zHat
    basis = {
        "lateral": np.array([1.0, 0.0, 0.0]),  # x: scattering vector at zero angles
        "longitudinal": np.array([0.0, 1.0, 0.0]),  # y: along the beam
        "vertical": np.array([0.0, 0.0, 1.0]),  # z: vertical
    }
    ZHAT_BL = np.array([0.0, 0.0, 1.0])
    YHAT_BL = np.array([0.0, 1.0, 0.0])
    stages = [
        # Sample stack: omega sits on the floor, chi and phi stack above it.
        # two_theta is independent of the sample stack: it also sits on the
        # floor and shares the same vertical axis as omega but is mechanically
        # decoupled (Busing & Levy Fig. 1; same relationship as S1-1/S2-1 in
        # our equipment description).
        Stage("omega", -ZHAT_BL, parent=None, role="sample"),
        Stage("chi", +YHAT_BL, parent="omega", role="sample"),
        Stage("phi", -ZHAT_BL, parent="chi", role="sample"),
        Stage("two_theta", -ZHAT_BL, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="fourc",
        stages=stages,
        basis=basis,
        description="Busing & Levy (1967) four-circle Eulerian diffractometer",
    )


def geometry_sixc():
    """
    Lohmeier & Vlieg (1993) six-circle surface diffractometer (sixc geometry).

    Basis: xHat=vertical, yHat=longitudinal (beam), zHat=lateral.
    (Lohmeier & Vlieg: y along beam, z horizontal -- same as You (1999) but
    with z horizontal rather than vertical; their z is our lateral.)

    The key structural difference from psic is that both sample and detector
    stacks share a common base stage (alpha, the rotary table), making this
    a coupled (S3D2)1 geometry rather than the decoupled S4D2 of psic.

    Stack (floor first):
        alpha (shared): vertical, right-handed  [rotary table]
          --> omega (sample): lateral, right-handed
                --> chi:  lateral,  right-handed
                      --> phi: lateral, right-handed
          --> delta (detector): lateral, left-handed
                --> gamma: vertical, right-handed

    Reference: M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
    """
    stages = [
        # Shared base (rotary table) -- treated as sample for rotation product
        Stage("alpha", +XHAT, parent=None, role="sample"),
        # Sample stack on top of alpha
        Stage("omega", +YHAT, parent="alpha", role="sample"),
        Stage("chi", +YHAT, parent="omega", role="sample"),
        Stage("phi", +YHAT, parent="chi", role="sample"),
        # Detector stack also rooted at alpha
        Stage("delta", -ZHAT, parent="alpha", role="detector"),
        Stage("gamma", +XHAT, parent="delta", role="detector"),
    ]
    return AdHocDiffractometer(
        name="sixc",
        stages=stages,
        description=(
            "Lohmeier & Vlieg (1993) six-circle surface diffractometer. "
            "Sample and detector share the alpha (rotary table) base stage."
        ),
    )


def lattice_vectors(a, b, c, alpha, beta, gamma):
    """
    Compute the three Cartesian direct lattice vectors from crystal lattice parameters.

    The convention used here places a1 along xHat, a2 in the xHat-yHat plane,
    and a3 determined by the right-hand rule.

    Parameters
    ----------
    a, b, c : float
        Lattice parameters in Angstroms.
    alpha, beta, gamma : float
        Lattice angles in degrees.
        alpha = angle between b and c axes
        beta  = angle between a and c axes
        gamma = angle between a and b axes

    Returns
    -------
    a1, a2, a3 : numpy.ndarray, shape (3,)
        Cartesian direct lattice vectors in Angstroms.
    """
    alpha_r = np.deg2rad(alpha)
    beta_r = np.deg2rad(beta)
    gamma_r = np.deg2rad(gamma)

    # a1 along xHat
    a1 = np.array([a, 0.0, 0.0])

    # a2 in the xHat-yHat plane
    a2 = np.array([b * np.cos(gamma_r), b * np.sin(gamma_r), 0.0])

    # a3 determined by alpha, beta, gamma
    a3x = c * np.cos(beta_r)
    a3y = c * (np.cos(alpha_r) - np.cos(beta_r) * np.cos(gamma_r)) / np.sin(gamma_r)
    a3z = np.sqrt(max(c**2 - a3x**2 - a3y**2, 0.0))
    a3 = np.array([a3x, a3y, a3z])

    return a1, a2, a3


def reciprocal_vectors(a1, a2, a3):
    """
    Compute the three reciprocal lattice vectors from Cartesian direct lattice vectors.

    Uses the standard crystallographic definition (Lecture-2-Reciprocal-lattice-notes):

        b1 = 2*pi * (a2 x a3) / (a1 . (a2 x a3))
        b2 = 2*pi * (a3 x a1) / (a1 . (a2 x a3))
        b3 = 2*pi * (a1 x a2) / (a1 . (a2 x a3))

    where the denominator is the volume of the direct unit cell: Vc = |a1 . (a2 x a3)|

    Parameters
    ----------
    a1, a2, a3 : numpy.ndarray, shape (3,)
        Cartesian direct lattice vectors in Angstroms.

    Returns
    -------
    b1, b2, b3 : numpy.ndarray, shape (3,)
        Reciprocal lattice vectors in inverse Angstroms (with 2*pi factor).
    """
    Vc = np.dot(a1, np.cross(a2, a3))  # unit cell volume

    b1 = 2 * np.pi * np.cross(a2, a3) / Vc
    b2 = 2 * np.pi * np.cross(a3, a1) / Vc
    b3 = 2 * np.pi * np.cross(a1, a2) / Vc

    return b1, b2, b3


def b_matrix(b1, b2, b3):
    """
    Compute the B matrix from the reciprocal lattice vectors.

    The B matrix transforms Miller indices h = (h, k, l) to Cartesian
    coordinates in the crystal frame (Busing & Levy, 1967, eq. 3):

        hc = B . h

    Following the I16 diffractometer convention:

        (b1, b2, b3) = 2*pi * B.T

    so the columns of B.T are the reciprocal lattice vectors divided by 2*pi,
    or equivalently the rows of B are b1, b2, b3 divided by 2*pi.

    B is not in general orthonormal (the crystal need not be cubic).

    Parameters
    ----------
    b1, b2, b3 : numpy.ndarray, shape (3,)
        Reciprocal lattice vectors in inverse Angstroms (with 2*pi factor).

    Returns
    -------
    B : numpy.ndarray, shape (3, 3)
        B matrix in inverse Angstroms (no 2*pi factor).
    """
    # columns of B.T are b1, b2, b3 / (2*pi)
    # so B.T = np.column_stack([b1, b2, b3]) / (2*pi)
    # and B = (np.column_stack([b1, b2, b3]) / (2*pi)).T
    B = np.column_stack([b1, b2, b3]).T / (2 * np.pi)
    return B
