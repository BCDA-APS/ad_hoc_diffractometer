(howto-basis-vectors)=
# Choose and Understand Basis Vectors

A **basis** maps the three observable physical directions of the
diffractometer — vertical, longitudinal, and transverse — to Cartesian
unit vectors (x, y, z).  The physics does not depend on which mapping
you choose; the basis is purely a **notational convention** that
determines the numerical representation of axis vectors, rotation
matrices, and the UB matrix.

This page explains what a basis is, why it matters, how to choose one,
and how to define a custom basis.

---

## The three physical directions

Every diffractometer has three mutually perpendicular directions that
can be identified by looking at the instrument:

| Direction | Definition |
|---|---|
| **vertical** | opposite to gravitational acceleration (upward) |
| **longitudinal** | in the horizontal plane, conventionally along the incident beam toward the equipment |
| **transverse** | orthogonal to both; positive sense completes a right-handed system (vertical x longitudinal) |

These directions are physical observables — they do not depend on any
coordinate convention.  See {doc}`../problem1` for the full derivation.

---

## What a basis does

A basis assigns each physical direction to a Cartesian axis.  For
example, the You (1999) convention assigns:

```
vertical     -> +x   (XHAT)
longitudinal -> +y   (YHAT)
transverse   -> +z   (ZHAT)
```

The Busing & Levy (1967) convention makes a different choice:

```
vertical     -> +z   (ZHAT)
longitudinal -> +y   (YHAT)
transverse   -> +x   (XHAT)
```

The NeXus [convention](https://manual.nexusformat.org/design.html#the-nexus-coordinate-system) makes a different choice:

```
vertical     -> +y   (YHAT)
longitudinal -> +z   (ZHAT)
transverse   -> +x   (XHAT)
```

All are valid right-handed systems.  The same physical rotation (e.g.
"right-handed about the vertical axis") has different matrix
representations in each basis, but the diffraction calculation gives the
same motor angles regardless.

---

## Requirements for a valid basis

The package validates the basis at geometry construction time.  A valid
basis must satisfy:

1. **Exactly three vectors** — one for each physical direction.
2. **Each vector is 3-dimensional** and non-zero.
3. **Mutually orthogonal** — every pair of vectors has zero dot product
   (within a tolerance of 1e-10).

```{note}
The package does **not** check right-handedness.  It is the caller's
responsibility to ensure that `vertical x longitudinal = transverse`
(i.e. the three vectors form a right-handed system).  All pre-built
geometries satisfy this constraint; custom geometries should verify it
explicitly.
```

---

## Pre-built conventions

The package exports two basis constants used by the preset geometries:

| Basis | vertical | longitudinal | transverse | Presets |
|---|---|---|---|---|
| {data}`~ad_hoc_diffractometer.factories.BASIS_YOU` | +x | +y | +z | psic, sixc, kappa6c, zaxis, s2d2, fivec |
| {data}`~ad_hoc_diffractometer.factories.BASIS_BL` | +z | +y | +x | fourcv, fourch, kappa4cv, kappa4ch |

Other conventions exist in the broader community (e.g. NeXus, Hkl) —
see {ref}`concepts` for a tabulated comparison.

---

## Defining a custom basis

You may pass any valid basis dict to
{class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`.  The keys
must be `"vertical"`, `"longitudinal"`, and `"transverse"`; the values
are 3-element array-like objects:

```python
import numpy as np
import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.stage import Stage

# NeXus convention: vertical=+y, longitudinal=+z, transverse=+x
nexus_basis = {
    "vertical":     np.array([0.0, 1.0, 0.0]),
    "longitudinal": np.array([0.0, 0.0, 1.0]),
    "transverse":   np.array([1.0, 0.0, 0.0]),
}

# Use this basis with a custom geometry
g = ahd.AdHocDiffractometer(
    name="my_geometry",
    stages=[
        Stage("omega", np.array([0.0, 1.0, 0.0]), role="sample"),
        Stage("ttheta", np.array([0.0, 1.0, 0.0]), role="detector"),
    ],
    basis=nexus_basis,
    description="Minimal two-circle with NeXus basis",
)
```

The axis vectors on each {class}`~ad_hoc_diffractometer.stage.Stage`
are expressed in the same Cartesian frame as the basis — so when you
change the basis, the Stage axis vectors must change accordingly.

---

## What changes when you change the basis

Changing the basis changes the **numerical representation** of every
vector quantity in the geometry:

- Stage axis vectors (e.g. `+XHAT` vs `+ZHAT` for the vertical axis)
- Rotation matrices
- The B, U, and UB matrices
- Q-vectors computed by {func}`~ad_hoc_diffractometer.orientation.angles_to_phi_vector`

What does **not** change:

- Motor angles returned by `forward()` and accepted by `inverse()`
- Bragg angles and d-spacings
- Physical rotation directions (right-handed vs left-handed)
- The number and type of solutions

The basis is a bookkeeping convention.  Two geometries that describe the
same physical diffractometer with different bases will produce identical
motor angles for the same reflection.

---

## Common mistakes

**Swapping vertical and transverse.**
If you accidentally assign `vertical -> +z` and `transverse -> +z` (or
repeat any axis), the constructor will raise a `ValueError` because the
vectors are not mutually orthogonal.

**Left-handed system.**
If `vertical x longitudinal` points in the opposite direction from your
`transverse` vector, the system is left-handed.  The package does not
detect this automatically — it will produce incorrect results silently.
Always verify: `np.cross(vertical, longitudinal)` should equal
`transverse` (not `-transverse`).

**Mixing conventions.**
If you define stages using axis vectors from one convention but pass a
basis dict from another, the geometry description will be internally
inconsistent.  The axis vectors and basis must use the same Cartesian
frame.

---

## See also

- {ref}`concepts` — coordinate conventions, axis signs, and four
  tabulated basis mappings
- {doc}`../problem1` — the case study that defines the physical
  reference frame
- {doc}`../problem2` — worked example deriving B/U/UB from a basis
  assignment
- {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer` — the
  `basis` constructor argument
- {data}`~ad_hoc_diffractometer.factories.BASIS_YOU`,
  {data}`~ad_hoc_diffractometer.factories.BASIS_BL` — pre-built basis
  constants
