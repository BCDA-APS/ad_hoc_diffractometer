(concepts)=
# Concepts

This page explains the key ideas behind `ad_hoc_diffractometer`.
Understanding these concepts helps you use the package effectively and
interpret its outputs correctly.

---

## Coordinate convention

The package uses a right-handed Cartesian frame in which physical directions
are mapped to named basis vectors.  The default convention follows You (1999):

| Basis vector | Physical direction | Lab meaning |
|---|---|---|
| `XHAT` (+x) | vertical | out of the floor |
| `YHAT` (+y) | longitudinal | along the beam, toward equipment |
| `ZHAT` (+z) | lateral | to our left when facing the equipment |

This is the convention used by the psic, sixc, kappa, zaxis, s2d2, and fivec
factories.  The Busing & Levy (1967) convention (used by fourcv, fourch,
kappa4cv, kappa4ch) swaps x and z:

| Basis vector | Physical direction |
|---|---|
| `XHAT` (+x) | lateral |
| `YHAT` (+y) | longitudinal |
| `ZHAT` (+z) | vertical |

Both conventions are supported through the `basis` argument to each factory
function.  The `BASIS_YOU` and `BASIS_BL` constants provide the two standard
dictionaries.

---

## Axis sign convention

Each stage has a rotation axis described as a signed unit vector.  The sign
encodes the handedness of positive rotation:

- `+nHat` — right-handed rotation about `nHat` (thumb along `+nHat`, fingers
  curl in the positive direction)
- `-nHat` — left-handed rotation about `nHat` (equivalent to right-handed
  with a negated angle)

For example, `eta` in psic has axis `−z` (lateral), meaning positive eta
rotates left-handed about the lateral axis — the same physical sense as
You (1999) equation 7.

Axis labels accepted by `parse_axis()`:

| String | Meaning |
|---|---|
| `"+x"`, `"+y"`, `"+z"` | right-handed about Cartesian axis |
| `"-x"`, `"-y"`, `"-z"` | left-handed about Cartesian axis |
| `"vertical"`, `"+vertical"`, `"-vertical"` | resolved against basis dict |
| `"longitudinal"`, `"lateral"` (± prefix) | resolved against basis dict |

---

## Stage stacking

Stages are stacked: each stage sits on its parent and its rotation modifies
the orientation of everything above it.  The `parent` attribute of each
{class}`~ad_hoc_diffractometer.stage.Stage` names the stage directly below
it (`None` for floor-mounted stages).

The sample rotation matrix at a given set of motor angles is the ordered
product of all sample-stack rotation matrices, from floor (outermost) to
innermost:

$$\mathbf{Z} = R_{\text{floor}} \cdot R_{\text{next}} \cdots R_{\text{innermost}}$$

The detector rotation matrix is computed the same way for detector-stack stages.

---

## Monochromatic radiation

The package assumes **monochromatic radiation** throughout.  All diffraction
calculations — Bragg angles, Q-vector magnitudes, forward and inverse
problems — are performed at a single fixed wavelength set on the geometry:

```python
g.wavelength = 1.5406  # Å  (Cu Kα)
```

Energy and wavelength are related by:

$$E \,[\text{keV}] = \frac{hc}{\lambda \,[\text{Å}]} \approx \frac{12.3984}{\lambda}$$

where $hc = 12.3984\,\text{keV·Å}$ exactly (2019 SI redefinition).
See {func}`~ad_hoc_diffractometer.wavelength_to_energy` and
{func}`~ad_hoc_diffractometer.energy_to_wavelength`.

---

## The B, U, and UB matrices

These three matrices connect Miller indices to motor angles.

### B matrix

The **B matrix** encodes the reciprocal lattice.  It transforms Miller
indices $\mathbf{h} = (h, k, l)^\top$ to the scattering vector
$\mathbf{Q}_c$ in Cartesian crystal-frame coordinates
(Busing & Levy 1967, eq. 3):

$$\mathbf{Q}_c = \mathbf{B}\,\mathbf{h}$$

$|\mathbf{Q}_c| = 2\pi / d_{hkl}$.  B is constructed from the unit-cell
parameters $a, b, c, \alpha, \beta, \gamma$ and is not in general
orthonormal.  See {class}`~ad_hoc_diffractometer.lattice.Lattice`.

### U matrix

The **U matrix** is orthonormal.  It relates the Cartesian crystal frame
to the phi-axis frame (the frame of the innermost sample stage):

$$\mathbf{Q}_\phi = \mathbf{U}\,\mathbf{Q}_c = \mathbf{U}\,\mathbf{B}\,\mathbf{h}$$

U encodes how the crystal is mounted: when all motor angles are zero, U
corrects for any misalignment between the crystal axes and the diffractometer
axes.  U is determined from orienting reflections.

### UB matrix

The **UB matrix** is the practical product $\mathbf{U}\mathbf{B}$.  It maps
Miller indices directly to the phi-axis frame and can be determined from
reflections alone, without prior knowledge of the lattice parameters
(Busing & Levy 1967, eqs. 29–31):

$$\mathbf{UB} = \mathbf{H}_c\,\mathbf{H}^{-1}$$

where $\mathbf{H}_c$ and $\mathbf{H}$ are matrices of observed and indexed
reflection vectors.

| Symbol | Name | Role |
|---|---|---|
| B | B matrix | Encodes lattice; maps hkl → crystal Cartesian |
| U | U matrix | Encodes crystal mounting; orthonormal |
| UB | UB matrix | Maps hkl → phi-axis frame; determined from reflections |

---

## Diffraction modes

A **diffraction mode** constrains which motor angles are free during a
forward calculation.  Each mode fixes one or more stages (frozen angles),
or defines a relationship between stages (e.g. bisecting: ω = 2θ/2).

Available modes depend on the geometry.  For four-circle geometries:

| Mode | Description |
|---|---|
| `bisecting` | ω = ttheta/2; chi and phi free |
| `fixed_chi` | chi held at a fixed value |
| `fixed_phi` | phi held at a fixed value |

Set and query the active mode via:

```python
g.mode_name = "bisecting"
print(g.mode_name)
```

When no mode is set (`mode_name` is `None`), all stages are free and the
forward solver explores the full solution space.

---

## The forward problem

The **forward problem** finds the motor angles that satisfy the Bragg
condition for a given reflection $(h, k, l)$:

$$\mathbf{Z} \cdot \mathbf{q}_\phi = \hat{\mathbf{Q}}_\text{lab}$$

where $\mathbf{Z}$ is the sample rotation matrix, $\mathbf{q}_\phi =
\mathbf{UB}\,\mathbf{h} / |\mathbf{UB}\,\mathbf{h}|$ is the unit
scattering vector in the phi frame, and $\hat{\mathbf{Q}}_\text{lab}$
points along the bisector of the incident and diffracted beams.

Multiple solutions may exist (different chi branches, different phi
settings).  The active diffraction mode filters these to the
physically useful subset.

```python
solutions = g.forward(1, 0, 0)
for s in solutions:
    print(s)
```

---

## The ψ angle

Two definitions of the azimuthal angle ψ appear in the literature.

**You (1999) ψ** (`geometry.psi()`): the azimuthal angle of the reference
vector **n** about **Q**, measured from the projection of the beam direction
onto the Q-perpendicular plane.  For a given (hkl, UB, **n**) this value is
*constant* across all motor-angle solutions satisfying the Bragg condition —
it is a crystal-orientation diagnostic, not a motor-angle observable.

**Busing & Levy (1967) ψ** (`psi_trajectory()`): the angle through which
the sample has been rotated about the scattering vector **Q** relative to a
chosen reference orientation.  This is the quantity physically varied during
a ψ scan on a real diffractometer.

See {func}`~ad_hoc_diffractometer.psi_trajectory` and
{meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.psi`.
