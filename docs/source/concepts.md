(concepts)=
# Concepts

Key ideas behind `ad_hoc_diffractometer`.  Each section gives a brief
overview and links to richer detail in the how-to guides and background
pages.

---

## Coordinate convention

Diffractometer stages are described in terms of three **observable physical
directions** that can be identified directly in the laboratory:

| Physical direction | Lab meaning |
|---|---|
| **vertical** | opposite to gravitational acceleration |
| **longitudinal** | a chosen direction in the plane perpendicular to vertical, conventionally aligned with the nominal incident beam; a property of the instrument installation |
| **lateral** | orthogonal to both; positive sense completes a right-handed system (vertical × longitudinal) |

The package uses a right-handed Cartesian frame internally.  Different authors
assigned different Cartesian letters (x, y, z) to these physical directions —
historically a source of confusion when diffractometer geometries are compared.
The package supports both major conventions via the `basis` argument to each
factory.

::::{tab-set}

:::{tab-item} You1999 (default)
Used by: `psic`, `sixc`, `kappa6c`, `zaxis`, `s2d2`, `fivec`

| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +x | `XHAT` |
| longitudinal | +y | `YHAT` |
| lateral | +z | `ZHAT` |

Pass `basis=BASIS_YOU` (the default for these geometries).
:::

:::{tab-item} BL1967
Used by: `fourcv`, `fourch`, `kappa4cv`, `kappa4ch`

Convention of Busing & Levy.

Also used by:
- [SPEC](https://certif.com)

| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +z | `ZHAT` |
| longitudinal | +y | `YHAT` |
| lateral | +x | `XHAT` |

Pass `basis=BASIS_BL` (the default for these geometries).
:::

:::{tab-item} NeXus
Used by: [NeXus](https://manual.nexusformat.org/design.html#the-nexus-coordinate-system)

Also used by:
- [hklpy2](https://blueskyproject.io/hklpy2/)

| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +y | `YHAT` |
| longitudinal | +z | `ZHAT` |
| lateral | +x | `XHAT` |
:::

:::{tab-item} Hkl
Used by: [Hkl](https://people.debian.org/~picca/hkl/hkl.html#org4569ec8)

| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +z | `ZHAT` |
| longitudinal | +x | `XHAT` |
| lateral | +y | `YHAT` |
:::

::::

The `BASIS_YOU` and `BASIS_BL` constants are exported from the package.

---

## Axis sign convention

Each stage's rotation axis is a **signed unit vector**: `+nHat` means
right-handed rotation, `-nHat` means left-handed (equivalent to
right-handed about the negated axis).  Physical direction names
(`"vertical"`, `"lateral"`, `"longitudinal"`) are resolved against
the geometry's basis dict.

See {func}`~ad_hoc_diffractometer.axes.parse_axis`.

---

## Stage stacking

Stages are stacked: each stage sits on its parent and its rotation modifies
the orientation of everything above it.  The `parent` attribute names the
stage directly below (`None` for floor-mounted stages).  The combined sample
rotation matrix is the ordered product from floor to innermost stage.

See {class}`~ad_hoc_diffractometer.stage.Stage`.

---

## Monochromatic radiation

The package assumes **monochromatic radiation** throughout — all
calculations are performed at a fixed wavelength.  Energy and wavelength
are related by $hc = 12.3984\,\text{keV·Å}$ exactly (2019 SI redefinition).

```python
g.wavelength = 1.5406  # Å  (Cu Kα)
```

See {doc}`howto/wavelength` and {mod}`~ad_hoc_diffractometer.radiation`.

---

## The B, U, and UB matrices

Three matrices connect Miller indices to motor angles:

| Symbol | Name | Role |
|---|---|---|
| **B** | B matrix | Encodes the reciprocal lattice; maps hkl → crystal Cartesian frame |
| **U** | U matrix | Orthonormal; encodes crystal mounting on the diffractometer |
| **UB** | UB matrix | Maps hkl → phi-axis frame; determined from orienting reflections |

The B matrix is constructed from unit-cell parameters $(a, b, c, \alpha,
\beta, \gamma)$.  U is determined by measuring two or more Bragg reflections.
UB = U × B maps Miller indices directly to the lab frame.

See {doc}`howto/orient`, {doc}`howto/lattice`, {doc}`problem2`, and
{class}`~ad_hoc_diffractometer.lattice.Lattice`.

---

## Diffraction modes

A **diffraction mode** describes how `forward()` will compute the motor
angles and which ones remain constant.  Modes fix or couple specific stages
(e.g. bisecting: ω = 2θ/2).  Available modes depend on the geometry.

```python
g.mode_name = "bisecting"
```

See {doc}`howto/modes` and {mod}`~ad_hoc_diffractometer.mode`.

---

## Forward and inverse computations

- **Forward** ({meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.forward`):
  given (h, k, l), find the motor angles satisfying the Bragg condition.
  Returns a **list** of 0 to ~12 solutions depending on geometry and mode.
- **Inverse** ({meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.inverse`):
  given motor angles, find the unique (h, k, l) in the Bragg condition.
  Requires a UB matrix.

See {doc}`howto/forward`.

---

## The ψ angle

Two definitions of ψ appear in the literature:

- **You (1999)**: azimuthal angle of a reference vector about **Q** —
  constant for a given (hkl, UB); a crystal-orientation diagnostic.
  See {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.psi`.
- **Busing & Levy (1967)**: angle of sample rotation about **Q** relative
  to a reference orientation — the quantity physically varied in a ψ scan.
  See {func}`~ad_hoc_diffractometer.psi_trajectory`.

See {doc}`howto/trajectory`.
