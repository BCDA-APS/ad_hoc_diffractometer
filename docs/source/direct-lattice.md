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

## Reference

- W.R. Busing & H.A. Levy, *Acta Cryst.* **22**, 457–464 (1967).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
