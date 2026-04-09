# Diffractometer with Specified Reference Frame

Given the equipment described in problem1.md, I now wish to use it with a different choice of basis vectors.  In this description:

- xHat: vertical
- yHat: longitudinal
- zHat: lateral

The sample will be a single crystal with the usual lattice constants: a, b, c, alpha, beta, gamma

The incident monochromatic radiation will have wavelength lambda

## Directions:

- Compare the geometry with this [article](./1999-JAppl-Cryst-32-614-623-H-You-psic-4S+2D/hn0093.pdf): 1999 J Appl Cryst 32 614-623.

## Analysis

### Basis Vector Assignment

The basis vectors specified here differ from those chosen in problem1.md.  The
new assignment is:

- xHat: vertical
- yHat: longitudinal (along the incoming beam direction)
- zHat: lateral

This is a valid right-handed system: xHat x yHat = zHat (vertical x longitudinal
= lateral).  This is the same coordinate convention used by You (1999).  From
section 2 of that paper: "the x axis is defined along the vertical mu and nu
axes and the y axis is defined along the incoming beam direction."

### Correspondence with You (1999)

The You (1999) paper describes a 4S+2D six-circle diffractometer: four
sample-orienting stages and two detector stages.  This matches our equipment
exactly (S2-1 through S2-4 for the sample stack, S1-1 and S1-2 for the detector
stack).  DOI: 10.1107/S0021889899001223

Using R to denote the reference rotation matrix with 1 on the diagonal identifying
the invariant axis:

- R xHat: 1 at position [1,1], rotation about the vertical axis
- R yHat: 1 at position [2,2], rotation about the longitudinal axis
- R zHat: 1 at position [3,3], rotation about the lateral axis

Right-handed rotation uses +sin below the diagonal; left-handed (equivalently,
rotation about the negated axis) uses -sin below the diagonal, i.e. +sin above.

The rotation matrices in You (1999) (equations 5, 7, 8) correspond to our stages
as follows:

| You angle | You matrix | Our stage | Physical axis | Axis vector |
|-----------|------------|-----------|---------------|-------------|
| mu        | U          | S2-1      | vertical      | +xHat       |
| eta       | X          | S2-2      | lateral       | -zHat       |
| chi       | H          | S2-3      | longitudinal  | +yHat       |
| phi       | M          | S2-4      | lateral       | -zHat       |
| nu        | P          | S1-1      | vertical      | +xHat       |
| delta     | D          | S1-2      | lateral       | -zHat       |

Notes:
- mu and nu share the same vertical axis (+xHat) and the same right-handed sense
  of rotation.  Their rotation axes are colinear; the stages are mechanically
  independent.
- eta, phi, and delta share the same lateral axis (-zHat) and the same
  left-handed sense of rotation.
- chi is the only stage with a longitudinal axis (+yHat), with right-handed
  rotation.

### Left-Handed Angles and Sign Convention

A left-handed angle about an axis is equivalent to a right-handed angle of
opposite sign about the same axis, or equivalently, a right-handed angle of the
same sign about the negated axis:

    R_left-handed(+nHat, theta) = R_right-handed(-nHat, theta)
                                 = R_right-handed(+nHat, -theta)

For stages where You (1999) uses a left-handed convention (eta, phi, delta), the
signed axis vector is given as -zHat rather than +zHat.  The physical rotation
axes are the same; only the sign convention for positive rotation differs.

### Crystallography Connection and Common Naming

The introduction of lattice constants a, b, c, alpha, beta, gamma and wavelength
lambda sets up the UB matrix formalism of Busing & Levy (1967), Acta Cryst. 22, 
457-464.

#### The B matrix

Busing & Levy define the B matrix (their equation 3) to transform Miller indices
h = (h,k,l) to Cartesian coordinates in the crystal frame:

    hc = B h

B is constructed from the reciprocal lattice parameters derived from a, b, c,
alpha, beta, gamma.  B is not in general orthonormal (the crystal need not be
cubic).

#### The U matrix

Busing & Levy define U (their equation 4) as the orthogonal matrix relating the
phi-axis system (attached to the innermost sample stage) to the crystal Cartesian
frame:

    h_phi = U hc

U is orthonormal and corrects for the misalignment between the crystal axes and
the diffractometer axes when all motor angles are zero.  Busing & Levy call U the
"orientation matrix", as do most other sources.  However, Walko (2016) explicitly
notes that UB is also sometimes called the "orientation matrix", making the term
ambiguous.  To avoid this ambiguity entirely, we do not use "orientation matrix"
at all.  We adopt the following unambiguous names:

| Symbol | Name we use | Meaning                                                     |
|--------|-------------|-------------------------------------------------------------|
| B      | B matrix    | maps Miller indices to crystal Cartesian coords;            |
|        |             | encodes a, b, c, alpha, beta, gamma                         |
| U      | U matrix    | orthonormal; relates crystal Cartesian frame to the         |
|        |             | phi-axis (innermost sample stage) frame                     |
| UB     | UB matrix   | maps Miller indices directly to the phi-axis frame;         |
|        |             | can be determined as a unit from three reflections          |

Busing & Levy treat UB as a single practical entity (their equations 29-31):

    UB = Hc H^{-1}

where Hc and H are matrices of observed and indexed reflection vectors
respectively.  This allows UB to be determined even when lattice parameters are
unknown.

#### Full diffraction equation

The full diffraction equation from You (1999) (equations 10-11) relates Miller
indices h = (h, k, l)^T to the sample stack rotation and the detector position
in the laboratory frame:

    h^M = M H X U_mu . UB . h

where:
- UB   is the UB matrix (U times B, as defined above)
- U_mu, X, H, M are the motor rotation matrices for mu, eta, chi, phi
- h^M  is the diffraction vector in the laboratory frame

The detector position is determined separately by the detector stack:

    kf = k P D kf0

where D and P are the motor rotation matrices for delta and nu respectively,
k = 2*pi/lambda is the wave number, and kf0 is the forward beam direction.

Note: You (1999) uses the symbol U for both the mu motor rotation matrix and
for the U matrix.  We use U_mu for the mu motor rotation to avoid ambiguity
with U (the U matrix).
