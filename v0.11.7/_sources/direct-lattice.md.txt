# Direct Lattice in Crystallography

This page summarises the vector mathematics underlying the direct lattice
representation used throughout `ad_hoc_diffractometer`.  See the
{mod}`~ad_hoc_diffractometer.lattice` module for the implementation.

## Lattice vectors

The unit cell of a crystal is described by three edge vectors
**A**, **B**, **C** expressed in a Cartesian frame.  With **A** along the
$x$-axis and **B** in the $xy$-plane, the standard crystallographic choice is:

| Vector | Cartesian components |
|--------|----------------------|
| **A** | $(a,\ 0,\ 0)$ |
| **B** | $(b\cos\gamma,\ b\sin\gamma,\ 0)$ |
| **C** | $(c\cos\beta,\ c\,\tfrac{\cos\alpha - \cos\beta\cos\gamma}{\sin\gamma},\ c\,v)$ |

where

$$v = \frac{\sqrt{1 - \cos^2\!\alpha - \cos^2\!\beta - \cos^2\!\gamma
      + 2\cos\alpha\cos\beta\cos\gamma}}{\sin\gamma}$$

and $a, b, c$ are the unit-cell edge lengths; $\alpha, \beta, \gamma$ are
the inter-edge angles ($\alpha$ between **B** and **C**, $\beta$ between
**A** and **C**, $\gamma$ between **A** and **B**).

## Unit cell volume

The unit cell volume is the scalar triple product of the three lattice vectors:

$$V = \mathbf{A} \cdot (\mathbf{B} \times \mathbf{C})$$

The cross product is evaluated first (standard vector-mathematics precedence),
giving a vector perpendicular to **B** and **C**; the dot product with **A**
then yields the signed volume of the parallelepiped spanned by the three
vectors.  With the lattice vectors defined above, this reduces to:

$$V = a\,b\,c\,v\,\sin\gamma$$

## Reciprocal lattice

The reciprocal lattice vectors $\mathbf{b}_1, \mathbf{b}_2, \mathbf{b}_3$
satisfy the orthogonality condition

$$\mathbf{b}_i \cdot \mathbf{a}_j = 2\pi\,\delta_{ij}$$

where $\mathbf{a}_1 = \mathbf{A}$, $\mathbf{a}_2 = \mathbf{B}$,
$\mathbf{a}_3 = \mathbf{C}$.  Each reciprocal vector includes the $2\pi$
factor (the convention used in `ad_hoc_diffractometer` and by Busing & Levy 1967).

## The B matrix

The **B matrix** encodes the reciprocal lattice vectors as its columns:

$$\mathbf{B} = [\mathbf{b}_1\ \mathbf{b}_2\ \mathbf{b}_3]$$

It maps Miller indices $\mathbf{h} = (h, k, l)^T$ to the scattering vector
in Cartesian crystal-frame coordinates (Busing & Levy 1967, eq. 3):

$$\mathbf{Q}_c = \mathbf{B}\,\mathbf{h}$$

The magnitude $|\mathbf{Q}_c| = |\mathbf{B}\,\mathbf{h}| = 2\pi / d_{hkl}$,
where $d_{hkl}$ is the interplanar spacing.

### Explicit construction

The B matrix is built from the reciprocal-lattice magnitudes
$a^\star$, $b^\star$, $c^\star$ and the reciprocal-cell angles
$\alpha^\star$, $\beta^\star$, $\gamma^\star$, using the explicit
upper-triangular form (eq. 3 of Busing & Levy 1967):

```
B[0,0] = a*               B[0,1] = b* cos γ*        B[0,2] = c* cos β*
B[1,0] = 0                B[1,1] = b* sin γ*        B[1,2] = -c* sin β* cos α
B[2,0] = 0                B[2,1] = 0                B[2,2] = 2π / c
```

This places the **reciprocal-lattice $\mathbf{a}^\star$** along
crystal-Cartesian $\hat{x}$, **$\mathbf{b}^\star$** in the
crystal-$xy$-plane (with positive $y$ component), and
**$\mathbf{c}^\star$** completing the right-handed orthonormal triple.
The B matrix is therefore upper-triangular by construction for every
crystal system.

The reciprocal-lattice magnitudes carry the $2\pi$ factor (BL1967
convention), giving $|\mathbf{B}\,\mathbf{h}| = 2\pi / d_{hkl}$ and
$\mathbf{b}_i \cdot \mathbf{a}_j = 2\pi\,\delta_{ij}$.

### Cross-solver equivalence

This construction is **numerically identical** (to machine precision)
to the B matrix used by

- **SPEC** (psic and other four-circle modes),
- **hkl_soleil** (the libhkl library) with its default
  `HKL_TAU = 2π` compile-time setting,
- **diffcalc-core** (`Crystal._set_reciprocal_cell`).

The equivalence holds for every crystal system, including
non-orthogonal cells (hexagonal, trigonal, monoclinic, triclinic).
For orthogonal cells (cubic, tetragonal, orthorhombic), the
direct-lattice vectors $\mathbf{a}$, $\mathbf{b}$, $\mathbf{c}$ are
parallel to the reciprocal vectors $\mathbf{a}^\star$,
$\mathbf{b}^\star$, $\mathbf{c}^\star$ respectively, so the BL1967
crystal-Cartesian frame coincides with the
"direct-lattice-$\mathbf{a}$-along-$\hat{x}$" frame.  For
non-orthogonal cells the two frames differ by a rotation in the
crystal-Cartesian space; both encode the same physical crystal, but
only the BL1967 frame is the one SPEC and hkl_soleil cross-validate
against.

The no-$2\pi$ alternative ($\mathbf{b}_i \cdot \mathbf{a}_j =
\delta_{ij}$, $|\mathbf{B}\,\mathbf{h}| = 1/d_{hkl}$), used by
FullProf, CrysFML, parts of CCP4, and `libhkl` compiled with
`HKL_TAU = 1`, is not exposed by this package.

## Reference

- W.R. Busing & H.A. Levy, *Acta Cryst.* **22**, 457–464 (1967).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
