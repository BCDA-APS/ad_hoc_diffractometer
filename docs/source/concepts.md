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

Each stage's rotation axis is a **signed unit vector**: `+n_axis` means
right-handed rotation, `-n_axis` means left-handed (equivalent to
right-handed about the negated axis).  Physical direction names
(`"vertical"`, `"transverse"`, `"longitudinal"`) are resolved against
the geometry's basis dict.

```{note}
``n_axis`` here is the **per-stage rotation axis** vector — internal
to the stage definition.  It is **not** the same as ``n̂`` (written as
the ``n_hat`` key in mode ``extras``, or as
``g.surface_normal`` / ``g.azimuth`` on the geometry),
which is the *user-facing reference vector* required by surface and
azimuthal modes.  See the {ref}`glossary <glossary>` entries for
"n̂ (reference vector)" and "Stage rotation axis" for the full
disambiguation, and {doc}`howto/surface` for how to provide the
reference vector at run time.
```

The two equivalences in full:

$$
R_\text{left-handed}(+\hat{n},\,\theta)
\;=\;
R_\text{right-handed}(-\hat{n},\,\theta)
\;=\;
R_\text{right-handed}(+\hat{n},\,-\theta).
$$

For stages where a published convention uses a left-handed sense
(for example eta, phi, and delta in You 1999 about the transverse
axis), the package stores the signed axis vector as $-\hat{n}$
rather than $+\hat{n}$.  The physical rotation axes are the same;
only the sign convention for positive rotation differs.

See {func}`~ad_hoc_diffractometer.axes.parse_axis`.

---

## Stage stacking

Stages are stacked: each stage sits on its parent and its rotation modifies
the orientation of everything above it.  The `parent` attribute names the
stage directly below (`None` for floor-mounted stages).  The combined sample
rotation matrix is the ordered product from floor to innermost stage.

See {class}`~ad_hoc_diffractometer.stage.Stage`.

---

## Declarative geometries

A diffractometer geometry can be described in two ways:

- **Declaratively**, as a YAML file that lists the basis, the stages, and
  the diffraction modes.  This is the recommended approach for any
  geometry that fits the schema.  The 10 demo geometries shipped under
  `src/ad_hoc_diffractometer/geometries/` are written in this form;
  they are *demonstrations* of the schema rather than authoritative
  descriptions of any production diffractometer.
- **Programmatically**, as a Python factory function that returns an
  {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`.
  Use this when stage axes depend on a runtime parameter the schema
  does not express, or when the stage layout itself is computed.

Every declarative file begins with the schema marker
``ad_hoc_diffractometer_geometry: {schema_revision: 1}``, which
identifies the document and selects the schema revision.  The loader
({func}`~ad_hoc_diffractometer.geometry_loader.load_geometry_file` and
{func}`~ad_hoc_diffractometer.geometry_loader.register_geometry_file`)
parses the YAML, validates it, and constructs the geometry.  See
{ref}`declarative-geometry-schema` for the complete reference and
{ref}`howto-custom-geometry` for a task-oriented walkthrough.

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
\beta, \gamma)$.  U is determined by measuring two or more Bragg
reflections.  UB = U × B maps Miller indices directly to the phi-axis
frame.

### B matrix

The B matrix (Busing & Levy 1967, eq. 3) transforms Miller indices
$\mathsf{h} = (h, k, l)^{\mathsf{T}}$ to the scattering vector in
Cartesian crystal-frame coordinates:

$$
\mathsf{Q}_c \;=\; \mathsf{B}\,\mathsf{h}.
$$

B is constructed from the reciprocal lattice parameters derived from
$(a, b, c, \alpha, \beta, \gamma)$ and is not in general orthonormal.
See {doc}`direct-lattice` for the explicit construction.  Note that
B depends only on the reciprocal lattice and is **independent of the
basis assignment** for the diffractometer axes.

### U matrix

The U matrix (Busing & Levy 1967, eq. 4) is the orthogonal matrix
relating the phi-axis frame (attached to the innermost sample stage)
to the crystal Cartesian frame:

$$
\mathsf{h}_\varphi \;=\; \mathsf{U}\,\mathsf{Q}_c \;=\; \mathsf{U}\,\mathsf{B}\,\mathsf{h}.
$$

U corrects for the misalignment between the crystal axes and the
diffractometer axes when all motor angles are zero.  Unlike B, U
**does** depend on the basis assignment — choosing a different
right-handed basis for the diffractometer axes rotates the U matrix
by the fixed signed-permutation matrix that relates the two bases.
The physical Bragg condition is invariant under that change; see
{doc}`problem2` for the worked case study.

### Nomenclature

To avoid the ambiguity noted by Walko (2016) — where both U and UB
are sometimes called the "orientation matrix" — this package uses
the following unambiguous names:

| Symbol | Name | Meaning |
|---|---|---|
| **B**  | B matrix  | Maps Miller indices to crystal Cartesian coords; encodes $a, b, c, \alpha, \beta, \gamma$ |
| **U**  | U matrix  | Orthonormal; relates the crystal Cartesian frame to the phi-axis frame |
| **UB** | UB matrix | Maps Miller indices directly to the phi-axis frame; determinable from reflections alone |

### UB as a practical entity

Busing & Levy treat UB as a single practical entity (eqs. 29–31):

$$
\mathsf{UB} \;=\; \mathsf{H}_c\,\mathsf{H}^{-1}
$$

where $\mathsf{H}_c$ and $\mathsf{H}$ are matrices of observed and
indexed reflection vectors respectively.  This allows UB to be
determined even when lattice parameters are unknown.

### Full diffraction equation

The full diffraction equation (You 1999, eqs. 10–11) relates Miller
indices $\mathsf{h}$ to the sample rotation matrices and the detector
position:

$$
\mathsf{h}^M \;=\; M\,H\,X\,U_\mu \cdot \mathsf{UB} \cdot \mathsf{h},
$$

where $U_\mu$, $X$, $H$, $M$ are the motor rotation matrices for
mu, eta, chi, and phi respectively, and $\mathsf{h}^M$ is the
diffraction vector in the laboratory frame.

The detector position is determined by:

$$
\mathsf{k}_f \;=\; k\,P\,D\,\mathsf{k}_{f0},
$$

where $D$ and $P$ are the rotation matrices for delta and nu,
$k = 2\pi/\lambda$ is the wave number, and $\mathsf{k}_{f0}$ is the
forward beam direction.

```{note}
You (1999) uses the symbol U for both the mu motor rotation matrix
and the U (orientation) matrix.  This package writes the mu motor
rotation as $U_\mu$ to avoid the ambiguity.
```

See {doc}`howto/orient`, {doc}`howto/lattice`,
{class}`~ad_hoc_diffractometer.lattice.Lattice`, and the
{doc}`problem2` case study showing that two different right-handed
basis assignments applied to the same equipment produce U/UB
matrices related by a fixed rotation while leaving the physics
invariant.

---

## Diffraction modes

A **diffraction mode** is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
that describes how `forward()` resolves the free degrees of freedom: which
stages are fixed, which are coupled, and which are solved freely.
Available modes depend on the geometry.

```python
# Four-circle geometries use "bisecting"
g = ahd.make_geometry("fourcv")
g.mode_name = "bisecting"

# Six-circle psic uses named variants
g = ahd.make_geometry("psic")
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
reference vector n̂ (surface normal, polarization axis, etc.):

```python
from ad_hoc_diffractometer import ReferenceConstraint

ReferenceConstraint("incidence", 5.0)  # incidence angle fixed
ReferenceConstraint("specular", True)    # incidence = emergence (symmetric)
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

Kappa geometries (``kappa4cv``, ``kappa4ch``, ``kappa6c``) have **real
motor angles** (``komega``, ``kappa``, ``kphi``) and **virtual
Eulerian pseudoangles** (``omega``, ``chi``, ``phi``) that are more
intuitive to specify.

### Geometry-aware decomposition

The conversion is derived directly from the preset's actual signed
stage axes via the rotation-matrix identity

$$
R(\hat{n}_{\kappa\varphi},\,\kappa\varphi) \cdot
R(\hat{n}_{\kappa},\,\kappa) \cdot
R(\hat{n}_{\kappa\omega},\,\kappa\omega)
\;=\;
R(\hat{n}_{\kappa\varphi},\,\varphi) \cdot
R(\hat{n}_{\chi,\,\text{eq}},\,\chi) \cdot
R(\hat{n}_{\kappa\omega},\,\omega).
$$

Each kappa preset declares the four signed axis vectors
``(n_komega, n_kappa, n_kphi, n_chi_eq)`` in a
{class}`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention`
instance attached to ``geometry.kappa_pseudo_angle_convention``.
The conversion functions
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes` and
{func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian_axes` solve the
identity above analytically — no Newton iteration is required:

```python
from ad_hoc_diffractometer.kappa import (
    eulerian_to_kappa_axes, kappa_to_eulerian_axes,
)

g = ahd.make_geometry("kappa4cv")
convention = g.kappa_pseudo_angle_convention

# Virtual Eulerian angles → real kappa motor angles (two branches)
komega, kappa, kphi = eulerian_to_kappa_axes(
    omega, chi, phi, convention, branch=+1
)

# Real kappa angles → virtual Eulerian pseudoangles
omega, chi, phi = kappa_to_eulerian_axes(komega, kappa, kphi, convention)
```

Branch selection: ``branch=+1`` (default) returns the kappa solution
with the smaller ``|κ|`` (the natural identity branch); ``branch=-1``
returns the chi-mirrored solution.

The kappa rotation axis itself is computed per preset to match the
published reference figures (issue #252).  In **physical-direction
language** the kappa axis is **inclined by α from the omega axis**,
lying in the plane that contains both omega and the
equivalent-Eulerian chi axis, and tilted from omega toward that chi
direction (Walko 2016 §4.1; ITC Vol. C §2.2.6.2: *"the κ axis is
inclined at 50° to the ω axis"*).  Per preset:

| Preset | omega line | chi-eq line | Axis vector |
|---|---|---|---|
| ``kappa4cv`` | transverse | +vertical | $\hat{n}_{\kappa} = \cos\alpha \cdot \hat{T} + \sin\alpha \cdot \hat{V}$ |
| ``kappa4ch`` | vertical | +longitudinal | $\hat{n}_{\kappa} = \cos\alpha \cdot \hat{V} + \sin\alpha \cdot \hat{L}$ |
| ``kappa6c``  | transverse | +vertical | $\hat{n}_{\kappa} = \cos\alpha \cdot \hat{T} + \sin\alpha \cdot \hat{V}$ (same as ``kappa4cv``) |

For ``kappa4cv`` and ``kappa6c`` this matches Walko (2016) Figure 3
and Thorkildsen *et al.* (2006) Table 1; for ``kappa4ch`` it
matches Wyckoff (1985) Figure 2(b) on p. 334.

Note that the kappa axis vector is not in general equal to
``cos α · n_komega + sin α · n_chi_eq``: the omega axis carries a
*signed handedness* (e.g. ``komega = −TRANSVERSE`` for left-handed
omega about the transverse line), while the kappa axis is computed
from the *unsigned* basis-direction lines so the kappa arm
physically extends into the upper half-space toward the published
direction.  The four-axis kappa pseudoangle convention stored in
{class}`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention`
records the actual signed stage axes; the closed-form solver in
``kappa.py`` works directly from those signed axes and does *not*
rely on the algebraic identity above.

### Divergence from Walko (2016) eq. [16]

The original Walko closed form

$$
\sin(\chi/2) = \sin(\kappa/2) \cdot \sin(\alpha_0), \qquad
\text{offset} = \arccos\bigl(\cos(\kappa/2)/\cos(\chi/2)\bigr)
$$

is correct **only for the axis convention assumed in Walko's
derivation** — omega about the transverse axis, chi about the
longitudinal axis, phi about the transverse axis, all with a
specific handedness.  No preset shipped with this package matches
that convention exactly: ``kappa4cv`` (BL) places ``komega`` along
``−TRANSVERSE``; ``kappa4ch`` (BL) along ``−VERTICAL``; ``kappa6c``
(You) along ``−TRANSVERSE`` with a horizontal ``mu`` base.  The
textbook formula therefore does **not** preserve the scattering
vector for any non-zero ``chi`` in any of these presets, which
manifested as silent ``"No solutions"`` returns from the kappa
virtual-angle solver (issue #241).

The textbook helpers
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa` and
{func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian` are retained
as reference implementations of the published closed form (with
deprecation warnings in their docstrings) but are **not** used
inside the solver.

### Handedness convention (Walko vs. ITC §2.2.6.2)

The presets shipped with this package follow Walko's left-handed
sign convention for ``omega/kphi/2theta`` (encoded as
``−TRANSVERSE`` or ``−VERTICAL`` in the ``Stage`` axis vectors).
ITC Vol. C §2.2.6.2 (2006) instead specifies the standard signs
of ``omega/chi/phi`` as right-handed (only ``2θ`` is left-handed
in Hamilton's choice).  The two conventions are equivalent up to
motor-angle sign flips and yield the same physical orientations.
Users who prefer the ITC convention can construct their own
geometry by negating the relevant ``Stage`` axis vectors — see
{func}`~ad_hoc_diffractometer.factories.register_geometry` for the
preset-construction API.

### Modes accept virtual angle names

Kappa modes accept the virtual angle names directly in
{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:

```python
from ad_hoc_diffractometer import ConstraintSet, SampleConstraint

g = ahd.make_geometry("kappa4cv")
# "chi" is a virtual angle — the kappa inversion solver handles it
g.modes["fixed_chi"] = ConstraintSet([SampleConstraint("chi", 90.0)])
```

See {func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes`,
{func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian_axes`,
{class}`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention`,
and the {doc}`howto/constraints` guide.

---

## Surface geometry and reference vector

Some diffraction modes and pseudo-angle functions require an external
reference vector supplied as **Miller indices (h, k, l)**.  Two separate
vectors may be set:

- **`surface_normal`** — direction perpendicular to the sample surface;
  used by incidence/exit angle functions and surface diffraction modes.
- **`azimuth`** — direction defining ψ = 0; used by
  `psi_angle` and `fixed_psi_*` modes.

```python
g.surface_normal = (0, 0, 1)    # (001)-cut sample
g.azimuth = (1, 0, 0)
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
