(problem2)=
# Case Study: Choice of Basis and the UB Matrix

This page works through the two questions posed at the end of the
[case study diffractometer](problem1.md) — the six-circle equipment
described there.  The objective is to demonstrate that the **choice of
basis is arbitrary**: different basis assignments produce numerically
different axis vectors and U/UB matrices that are related by a fixed
rotation, but the underlying physics — the conversion of motor angles
to (h, k, l) and back — is invariant.  The package's flexibility around
{data}`~ad_hoc_diffractometer.factories.BASIS_YOU`,
{data}`~ad_hoc_diffractometer.factories.BASIS_BL`, and any other
right-handed basis exists because of this invariance, not in spite of
it.

For background on the B, U, UB matrices themselves and the You (1999)
full diffraction equation cited below, see
{doc}`concepts`.

**References:**
- W.R. Busing & H.A. Levy, *Acta Cryst.* **22**, 457–464 (1967).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
- H. You, *J. Appl. Cryst.* **32**, 614–623 (1999).
  DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)

---

## Question 1 — assignment using the You (1999) basis

The You (1999) convention assigns the three observable physical
directions to Cartesian unit vectors as:

- **xHat**: vertical (along the mu and nu rotation axes)
- **yHat**: longitudinal (along the incoming beam direction)
- **zHat**: transverse

This is a valid right-handed system: xHat × yHat = zHat.
From You (1999) §2: *"the x axis is defined along the vertical mu and
nu axes and the y axis is defined along the incoming beam direction."*

The You (1999) paper describes the same 4S+2D six-circle geometry as
the case-study equipment: four sample-orienting stages (mu, eta, chi,
phi) and two detector stages (nu, delta).  Mapping the case-study
stages onto the You convention:

| You angle | You matrix | Case-study stage | Physical axis  | Axis vector |
|-----------|------------|------------------|----------------|-------------|
| mu        | U          | sample stage 1   | vertical       | +xHat       |
| eta       | X          | sample stage 2   | transverse     | −zHat       |
| chi       | H          | sample stage 3   | longitudinal   | +yHat       |
| phi       | M          | sample stage 4   | transverse     | −zHat       |
| nu        | P          | detector stage 1 | vertical       | +xHat       |
| delta     | D          | detector stage 2 | transverse     | −zHat       |

Notes on this assignment:

- mu and nu share a colinear vertical axis (+xHat) with the same
  right-handed sense of rotation; the stages are mechanically
  independent.
- eta, phi, and delta share the transverse axis (−zHat) — that is,
  left-handed rotation about +zHat.  The signed axis vector encodes
  the handedness; see {doc}`concepts` for the sign-convention
  identity.
- chi is the only stage with a longitudinal axis (+yHat),
  right-handed.

### Computing U in the You basis

Let $\mathsf{U}_\text{You}$ denote the orthogonal matrix that relates
the crystal Cartesian frame to the phi-axis frame in this basis (the
"U matrix" of Busing & Levy 1967, eq. 4).  The procedure to determine
$\mathsf{U}_\text{You}$ from the geometry above is:

1. Construct the B matrix from the sample lattice parameters
   $(a, b, c, \alpha, \beta, \gamma)$ as described in
   {doc}`direct-lattice` and Busing & Levy (1967, eq. 3).  B encodes
   the reciprocal-lattice metric only and is **independent of the
   basis assignment** for the diffractometer axes.
2. Measure two non-collinear reflections at known motor angles —
   each provides a Q-vector in the phi-axis frame computed from the
   product of the four sample-stage rotation matrices in the You
   basis (mu about +xHat, eta about −zHat, chi about +yHat, phi
   about −zHat) acting on the lab-frame scattering vector.
3. Solve for the orthogonal $\mathsf{U}_\text{You}$ that aligns the
   crystal Cartesian frame onto the phi-axis frame using Busing &
   Levy (1967, eqs. 27–28), or equivalently use the package helper
   {func}`~ad_hoc_diffractometer.orientation.ub_from_two_reflections`.

The resulting $\mathsf{UB}_\text{You} = \mathsf{U}_\text{You}\,
\mathsf{B}$ maps Miller indices directly to the phi-axis frame in the
You basis.

---

## Question 2 — assignment using a different basis

Any right-handed orthogonal mapping of the three physical directions
to Cartesian axes is also valid.  The Busing & Levy (1967) convention
is the natural alternative:

- **xHat**: transverse
- **yHat**: longitudinal
- **zHat**: vertical

This is again a valid right-handed system: xHat × yHat = zHat.
This is the convention used in {data}`~ad_hoc_diffractometer.factories.BASIS_BL`
and in the SPEC control system.

The **physical equipment is unchanged** — the same six rotation
stages, the same axes of rotation, the same handedness.  Only the
labels change.  The same case-study stages now have these axis
vectors:

| Stage            | Physical axis  | You-basis vector | BL-basis vector |
|------------------|----------------|------------------|-----------------|
| sample stage 1   | vertical       | +xHat (You)      | +zHat (BL)      |
| sample stage 2   | transverse     | −zHat (You)      | −xHat (BL)      |
| sample stage 3   | longitudinal   | +yHat (You)      | +yHat (BL)      |
| sample stage 4   | transverse     | −zHat (You)      | −xHat (BL)      |
| detector stage 1 | vertical       | +xHat (You)      | +zHat (BL)      |
| detector stage 2 | transverse     | −zHat (You)      | −xHat (BL)      |

The longitudinal direction is +yHat in both bases (so chi's axis
vector reads identically as +yHat); only the assignments of the
vertical and transverse directions to xHat and zHat are swapped.

### Computing U in the BL basis

The procedure is structurally identical to Question 1 — same B
matrix, same two reflections, same package helper — but every sample-
stage rotation is now constructed about the BL-basis axis vector
listed in the table above.  Call the resulting orthogonal matrix
$\mathsf{U}_\text{BL}$ and the product
$\mathsf{UB}_\text{BL} = \mathsf{U}_\text{BL}\,\mathsf{B}$.

Numerically $\mathsf{U}_\text{BL} \neq \mathsf{U}_\text{You}$, and
$\mathsf{UB}_\text{BL} \neq \mathsf{UB}_\text{You}$.  The two
matrices differ by exactly the rotation that carries one basis
into the other — see the next section.

---

## The invariance result

Define $\mathsf{R}_\text{You}^\text{BL}$ as the orthogonal matrix
that re-expresses any vector from the You Cartesian frame to the
BL Cartesian frame.  Because both bases are right-handed and assign
the same three physical directions to *some* permutation of
$\pm\hat{x}$, $\pm\hat{y}$, $\pm\hat{z}$,
$\mathsf{R}_\text{You}^\text{BL}$ is a fixed signed-permutation
matrix that depends only on the two basis dictionaries — not on the
sample, the lattice, the wavelength, or any motor angle.

The U and UB matrices in the two bases are related by a single
similarity-style transformation:

$$
\mathsf{U}_\text{BL}
\;=\;
\mathsf{R}_\text{You}^\text{BL}\,
\mathsf{U}_\text{You}\,
\bigl(\mathsf{R}_\text{You}^\text{BL}\bigr)^{\mathsf{T}}.
$$

The B matrix is unchanged because it depends only on the reciprocal
lattice, not on the diffractometer-axis basis, so

$$
\mathsf{UB}_\text{BL}
\;=\;
\mathsf{R}_\text{You}^\text{BL}\,
\mathsf{UB}_\text{You}\,
\bigl(\mathsf{R}_\text{You}^\text{BL}\bigr)^{\mathsf{T}}.
$$

The physical Bragg condition — the equation
$\mathsf{h}^M = M\,H\,X\,U_\mu \cdot \mathsf{UB} \cdot \mathsf{h}$
of You (1999, eqs. 10–11) — is therefore invariant under the basis
change: every motor rotation matrix on the left-hand side picks up
the same conjugation by $\mathsf{R}_\text{You}^\text{BL}$, and the
conjugations cancel.  In words:

> **The motor angles produced by `forward(h, k, l)` and the (h, k, l)
> recovered by `inverse(angles)` are the same regardless of which
> right-handed basis was used to construct the geometry.**

The choice of basis is therefore a representational convenience.  It
fixes which Cartesian symbols appear in the printed axis vectors and
in the printed U/UB matrices, but it does not affect the diffraction
calculation.  This is why
{class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
exposes a `basis` constructor argument and ships several pre-made
choices ({data}`~ad_hoc_diffractometer.factories.BASIS_YOU`,
{data}`~ad_hoc_diffractometer.factories.BASIS_BL`): users can pick
the convention they are most comfortable comparing against published
literature, without changing what the package computes.

---

## See also

- {doc}`problem1` — the case-study equipment described in detail
- {doc}`concepts` — B, U, UB matrices; sign convention; full
  diffraction equation
- {doc}`howto/basis_vectors` — practical tutorial on choosing a basis
- {data}`~ad_hoc_diffractometer.factories.BASIS_YOU`,
  {data}`~ad_hoc_diffractometer.factories.BASIS_BL`
- {func}`~ad_hoc_diffractometer.orientation.ub_from_two_reflections`
