# Case Study — Coordinate Convention and UB Matrix

This page analyzes the [case study diffractometer](problem1.md) to show
how a **basis vector assignment** (the mapping of physical directions to
Cartesian axes) leads to the B, U, and UB matrices used throughout
`ad_hoc_diffractometer`.

The worked example below uses the You (1999) convention, but the
procedure is identical for any right-handed orthogonal basis — only the
numerical values of the axis vectors change.  See
{doc}`howto/basis_vectors` for a tutorial on choosing a basis.

**Reference:** H. You, *J. Appl. Cryst.* **32**, 614–623 (1999).
DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)

---

## Basis vector assignment

Any right-handed orthogonal mapping of the three physical directions to
Cartesian unit vectors is valid.  The You (1999) convention assigns:

- **xHat**: vertical (along the mu and nu rotation axes)
- **yHat**: longitudinal (along the incoming beam direction)
- **zHat**: transverse

This is a valid right-handed system: xHat × yHat = zHat.
From You (1999) §2: *"the x axis is defined along the vertical mu and nu
axes and the y axis is defined along the incoming beam direction."*

## Correspondence with You (1999)

The You (1999) paper describes a 4S+2D six-circle geometry: four
sample-orienting stages (mu, eta, chi, and phi) and two detector stages (nu,
delta).  This matches the case study equipment exactly.

Using **R** to denote a rotation matrix whose invariant axis is identified
by a 1 on the diagonal:

- **R xHat**: 1 at position [1,1] — rotation about the vertical axis
- **R yHat**: 1 at position [2,2] — rotation about the longitudinal axis
- **R zHat**: 1 at position [3,3] — rotation about the transverse axis

Right-handed rotation places +sin below the diagonal; left-handed rotation
(equivalently, right-handed rotation about the negated axis) places +sin
above the diagonal.

| You angle | You matrix | Stage | Physical axis | Axis vector |
|-----------|------------|-------|---------------|-------------|
| mu        | U          | S2-1  | vertical      | +xHat       |
| eta       | X          | S2-2  | transverse       | −zHat       |
| chi       | H          | S2-3  | longitudinal  | +yHat       |
| phi       | M          | S2-4  | transverse       | −zHat       |
| nu        | P          | S1-1  | vertical      | +xHat       |
| delta     | D          | S1-2  | transverse       | −zHat       |

Notes:

- mu and nu share a colinear vertical axis (+xHat) with the same
  right-handed sense of rotation; the stages are mechanically independent.
- eta, phi, and delta share the transverse axis (−zHat) with left-handed
  rotation.
- chi is the only stage with a longitudinal axis (+yHat), right-handed.

## Sign convention

A left-handed rotation about an axis is equivalent to a right-handed
rotation about the negated axis:

```
R_left-handed(+nHat, θ) = R_right-handed(−nHat, θ)
                         = R_right-handed(+nHat, −θ)
```

For stages where You (1999) uses a left-handed convention (eta, phi,
delta), the signed axis vector is −zHat rather than +zHat.  The physical
rotation axes are the same; only the sign convention for positive rotation
differs.

## The B, U, and UB matrices

The introduction of lattice constants a, b, c, α, β, γ and a fixed
wavelength λ sets up the UB matrix formalism of Busing & Levy (1967).

### B matrix

The B matrix (Busing & Levy 1967, eq. 3) transforms Miller indices
**h** = (h, k, l)ᵀ to the scattering vector in Cartesian crystal-frame
coordinates:

```
Q_c = B h
```

B is constructed from the reciprocal lattice parameters derived from
a, b, c, α, β, γ.  B is not in general orthonormal.
See [Direct Lattice](direct-lattice.md) for the explicit construction.

### U matrix

The U matrix (Busing & Levy 1967, eq. 4) is the orthogonal matrix relating
the phi-axis frame (attached to the innermost sample stage) to the crystal
Cartesian frame:

```
h_phi = U Q_c = U B h
```

U corrects for the misalignment between the crystal axes and the
diffractometer axes when all motor angles are zero.

To avoid the ambiguity noted by Walko (2016) — where both U and UB are
sometimes called the "orientation matrix" — this package uses the
following unambiguous names:

| Symbol | Name | Meaning |
|--------|------|---------|
| B | B matrix | Maps Miller indices to crystal Cartesian coords; encodes a, b, c, α, β, γ |
| U | U matrix | Orthonormal; relates crystal Cartesian frame to the phi-axis frame |
| UB | UB matrix | Maps Miller indices directly to the phi-axis frame; determinable from reflections alone |

### UB as a practical entity

Busing & Levy treat UB as a single practical entity (eqs. 29–31):

```
UB = Hc H⁻¹
```

where **Hc** and **H** are matrices of observed and indexed reflection
vectors respectively.  This allows UB to be determined even when lattice
parameters are unknown.

## Full diffraction equation

The full diffraction equation (You 1999, eqs. 10–11) relates Miller indices
**h** to the sample rotation matrices and the detector position:

```
h^M = M H X U_mu · UB · h
```

where U_mu, X, H, M are the motor rotation matrices for mu, eta, chi, and phi
respectively, and **h**^M is the diffraction vector in the laboratory frame.

The detector position is determined by:

```
kf = k P D kf0
```

where D and P are the rotation matrices for delta and nu, k = 2π/λ is the
wave number, and **kf0** is the forward beam direction.

```{note}
You (1999) uses the symbol U for both the mu motor rotation matrix and the
U (orientation) matrix.  This package uses U_mu for the mu motor rotation
to avoid ambiguity.
```

## References

- W.R. Busing & H.A. Levy, *Acta Cryst.* **22**, 457–464 (1967).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
- H. You, *J. Appl. Cryst.* **32**, 614–623 (1999).
  DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)
- D.A. Walko, *Reference Module in Materials Science and Materials
  Engineering*, Elsevier (2016).
