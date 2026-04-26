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
| **transverse** | orthogonal to both; positive sense completes a right-handed system (vertical × longitudinal) |

The package uses a right-handed Cartesian frame internally.  Different authors
assigned different Cartesian letters (x, y, z) to these physical directions —
historically a source of confusion when diffractometer geometries are compared.
The package accepts any right-handed orthogonal basis via the `basis` argument
to {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`.

::::{tab-set}

:::{tab-item} You1999 (default)
| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +x | `XHAT` |
| longitudinal | +y | `YHAT` |
| transverse | +z | `ZHAT` |

Used by: `psic`, `sixc`, `kappa6c`, `zaxis`, `s2d2`, `fivec`

Pass `basis=BASIS_YOU` (the default for these geometries).
:::

:::{tab-item} BL1967
| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +z | `ZHAT` |
| longitudinal | +y | `YHAT` |
| transverse | +x | `XHAT` |

Convention of Busing & Levy.

Used by: `fourcv`, `fourch`, `kappa4cv`, `kappa4ch`

Also used by:
- [SPEC](https://certif.com)

Pass `basis=BASIS_BL` (the default for these geometries).
:::

:::{tab-item} NeXus
| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +y | `YHAT` |
| longitudinal | +z | `ZHAT` |
| transverse | +x | `XHAT` |

Used by: [NeXus](https://manual.nexusformat.org/design.html#the-nexus-coordinate-system)

Also used by:
- [hklpy2](https://blueskyproject.io/hklpy2/)
:::

:::{tab-item} Hkl
| Physical direction | Cartesian | Constant |
|---|---|---|
| vertical | +z | `ZHAT` |
| longitudinal | +x | `XHAT` |
| transverse | +y | `YHAT` |

Used by: [Hkl](https://people.debian.org/~picca/hkl/hkl.html#org4569ec8)
:::

::::

The `BASIS_YOU` and `BASIS_BL` constants are exported from the package.

---

## Axis sign convention

Each stage's rotation axis is a **signed unit vector**: `+nHat` means
right-handed rotation, `-nHat` means left-handed (equivalent to
right-handed about the negated axis).  Physical direction names
(`"vertical"`, `"transverse"`, `"longitudinal"`) are resolved against
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

A **diffraction mode** is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
that describes how `forward()` resolves the free degrees of freedom: which
stages are fixed, which are coupled, and which are solved freely.
Available modes depend on the geometry.

```python
# Four-circle geometries use "bisecting"
g = ahd.presets.fourcv()
g.mode_name = "bisecting"

# Six-circle psic uses named variants
g = ahd.presets.psic()
g.mode_name = "bisecting_vertical"   # vertical scattering plane
g.mode_name = "bisecting_horizontal" # horizontal scattering plane
```

Modes can also be added at run time:

```python
from ad_hoc_diffractometer import ConstraintSet, SampleConstraint

g.modes["my_chi45"] = ConstraintSet([SampleConstraint("chi", 45.0)])
g.mode_name = "my_chi45"
```

See {doc}`howto/modes`, {doc}`howto/constraints`, and
{mod}`~ad_hoc_diffractometer.mode`.

---

## Diffraction constraints

Specifying (h, k, l) provides exactly **3 equations** on the motor angles.
A geometry with N real axes therefore has **N − 3 free parameters** that must
each be resolved by a constraint.  Every mode is a
{class}`~ad_hoc_diffractometer.mode.ConstraintSet` — an ordered list of
constraints equal in length to N − 3.

Three constraint categories exist:

**Sample constraints** fix one sample motor angle at a declared value, or
express the bisecting relational condition:

```python
from ad_hoc_diffractometer import SampleConstraint, BisectConstraint

SampleConstraint("chi", 90.0)       # chi fixed at 90°
BisectConstraint("eta", "delta")    # eta = delta / 2  (psic bisecting)
```

**Detector constraints** fix one detector stage at a declared value, or
constrain the azimuthal angle of Q (the ``"qaz"`` pseudo-angle from You 1999
eq. 18: ``tan(qaz) = tan(delta) / sin(nu)``):

```python
from ad_hoc_diffractometer import DetectorConstraint

DetectorConstraint("nu", 0.0)       # nu fixed at 0°
DetectorConstraint("qaz", 90.0)     # Q in the vertical plane
```

**Reference constraints** express a condition between Q and an external
reference vector n̂ (surface normal, polarisation axis, etc.):

```python
from ad_hoc_diffractometer import ReferenceConstraint

ReferenceConstraint("alpha_i", 5.0)  # incidence angle fixed
ReferenceConstraint("a_eq_b", True)  # alpha_i = beta_out (symmetric)
```

Taxonomy rules: at most one {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`,
at most one {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`.

Two checks distinguish solver availability from prerequisite satisfaction:

- `constraint.is_implemented(geometry)` — returns `True` when a forward
  solver exists for this constraint on this geometry.
- `rc.has_reference_vector(geometry)` — returns `True` when the required
  n̂ vector is set on the geometry (a prerequisite for reference constraints,
  independent of solver availability).

See {doc}`howto/constraints` and {mod}`~ad_hoc_diffractometer.mode`.

---

## Kappa virtual angles

Kappa geometries (``kappa4cv``, ``kappa4ch``, ``kappa6c``) have **real motor
angles** (komega, kappa, kphi) and **virtual Eulerian pseudoangles** (omega,
chi, phi) that are more intuitive to specify.

The conversion follows Walko (2016) eq. [16] with a fixed tilt angle
α₀ (default 50°):

```python
from ad_hoc_diffractometer.kappa import kappa_to_eulerian, eulerian_to_kappa

# Real kappa angles → virtual Eulerian angles
omega, chi, phi = kappa_to_eulerian(komega, kappa, kphi, alpha_deg=50.0)

# Virtual Eulerian angles → real kappa angles (two branches)
komega, kappa, kphi = eulerian_to_kappa(omega, chi, phi, alpha_deg=50.0)
```

Branch selection: the positive branch (kappa ≥ 0) corresponds to positive
chi; the negative branch (kappa < 0) to negative chi.

Kappa modes accept virtual angle names directly in
{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:

```python
from ad_hoc_diffractometer import ConstraintSet, SampleConstraint, BisectConstraint

g = ahd.presets.kappa4cv()
# "chi" is a virtual angle — the kappa inversion solver handles it
g.modes["fixed_chi"] = ConstraintSet([SampleConstraint("chi", 90.0)])
```

See {func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian`,
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa`, and the
{doc}`howto/constraints` guide.

---

## Surface geometry and reference vector

Some diffraction modes and pseudo-angle functions require an external
reference vector supplied as **Miller indices (h, k, l)**.  Two separate
vectors may be set:

- **`surface_normal`** — direction perpendicular to the sample surface;
  used by incidence/exit angle functions and surface diffraction modes.
- **`azimuthal_reference`** — direction defining ψ = 0; used by
  `psi_angle` and `fixed_psi_*` modes.

```python
g.surface_normal = (0, 0, 1)    # (001)-cut sample
g.azimuthal_reference = (1, 0, 0)
```

Vectors are stored as Miller indices and converted to the lab frame
internally via the UB matrix.

See {doc}`howto/surface` and {mod}`~ad_hoc_diffractometer.reference`.

---

## Custom exceptions

Two exceptions signal specific failure modes of the forward solver:

{exc}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
: Raised when |Q| > 4π/λ — the requested reflection cannot be reached at
  the current wavelength regardless of motor angles.  Carries attributes
  `q_mag`, `q_max`, and `wavelength`.

{exc}`~ad_hoc_diffractometer.mode.ConstraintViolation`
: Raised when a solver returns a solution that violates a declared
  constraint beyond the display-precision tolerance (indicates a solver
  error or an unimplemented virtual-angle constraint).  Carries attributes
  `solution_index`, `constraint_repr`, `residual`, and `tolerance`.

```python
from ad_hoc_diffractometer import EwaldSphereViolation, ConstraintViolation

try:
    solutions = g.forward(10, 10, 10)   # likely unreachable
except EwaldSphereViolation as e:
    print(f"|Q| = {e.q_mag:.3f} Å⁻¹ exceeds Ewald sphere (max {e.q_max:.3f} Å⁻¹)")
```

See {doc}`howto/forward` and {mod}`~ad_hoc_diffractometer.mode`.

---

## Forward and inverse computations

- **Forward** ({meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward`):
  given (h, k, l), find the motor angles satisfying the Bragg condition.
  Returns a **list** of 0 to ~12 solutions depending on geometry and mode.
- **Inverse** ({meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.inverse`):
  given motor angles, find the unique (h, k, l) in the Bragg condition.
  Requires a UB matrix.

See {doc}`howto/forward`.

---

## The ψ angle

Two definitions of ψ appear in the literature:

- **You (1999)**: azimuthal angle of a reference vector about **Q** —
  constant for a given (hkl, UB); a crystal-orientation diagnostic.
  See {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.psi`.
- **Busing & Levy (1967)**: angle of sample rotation about **Q** relative
  to a reference orientation — the quantity physically varied in a ψ scan.
  See {func}`~ad_hoc_diffractometer.psi_trajectory`.

See {doc}`howto/trajectory`.

---

## Serialization

The complete diffractometer state — geometry, wavelength, lattice, reflections,
UB matrix, and all parameters — can be saved and restored via
`to_dict()` / `from_dict()` on {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`.
The dict contains only JSON-compatible types; save to JSON (stdlib) or YAML
(`pyyaml`) without loss.

See {doc}`howto/serialize`.
