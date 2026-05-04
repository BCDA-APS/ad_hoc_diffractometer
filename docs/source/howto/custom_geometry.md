(howto-custom-geometry)=
# Build a Custom Diffractometer Geometry

This guide walks through defining a diffractometer geometry that is not
one of the demo geometries.  Each design decision is explained as
it arises, with notes on what happens if you get it wrong.

The worked example builds a **three-circle powder diffractometer** with
a theta-2theta arm and a sample spinner (phi).  This is simple enough to
follow but exercises every decision: basis, axis signs, parent chains,
roles, and modes.

```{tip}
If your instrument matches one of the pre-built geometries
({func}`~ad_hoc_diffractometer.presets.psic`,
{func}`~ad_hoc_diffractometer.presets.fourcv`, etc.), use the factory
function instead.  This guide is for instruments that differ from the
:ref:`demo geometries <geometries>`.
```

---

## The instrument

Imagine a laboratory powder diffractometer.  The designer has chosen a
**horizontal scattering plane** geometry (the reasons — whatever they may be —
are out of scope here).

That single design decision determines the rotation axis of the scattering arm.
The scattering plane is the plane containing the incident beam k_i, the
scattered beam k_f, and the scattering vector Q = k_f − k_i.  Because k_f sweeps
in the plane perpendicular to the 2θ rotation axis, **the scattering plane is
always perpendicular to the 2θ axis**.  A horizontal scattering plane therefore
*requires* 2θ — and the matching sample stage θ that tracks half of it — to
rotate about the **vertical** axis.  (Conversely, the
{func}`~ad_hoc_diffractometer.presets.fourcv` geometry places 2θ about a
horizontal/transverse axis to obtain a vertical scattering plane; see
{doc}`/concepts` for the wider discussion of the ``v`` / ``h`` suffix
convention.)

With that constraint fixed, the three rotary stages of this instrument
are:

1. **theta** — rotates the sample about the vertical axis
   (right-handed: positive angles rotate the sample face toward the
   beam, conventional laboratory sense).
2. **phi** — a sample spinner mounted on top of theta, rotating about
   the axis perpendicular to the sample face (which, when theta = 0,
   points along the beam).  Right-handed.
3. **ttheta** — the detector arm, rotating about the same vertical axis
   as theta but mechanically independent.  Right-handed.

---

## Step 1 — Choose a basis

The basis maps the three physical directions to Cartesian axes.  Any
right-handed orthogonal mapping works (see {doc}`basis_vectors`).

For a laboratory instrument the Busing & Levy convention is natural:

```python
import numpy as np
import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.stage import Stage

basis = BASIS_BL
VERTICAL     = basis["vertical"]      # +z  [0, 0, 1]
LONGITUDINAL = basis["longitudinal"]  # +y  [0, 1, 0]
TRANSVERSE   = basis["transverse"]    # +x  [1, 0, 0]
```

```{note}
The basis choice does not affect the physics — only the numerical
representation of axis vectors and matrices.  Choosing a different
basis produces the same motor angles from `forward()`.  See
{doc}`basis_vectors` for details.
```

---

## Step 2 — Determine axis vectors from physical observation

Stand at the instrument and observe each stage:

1. **What physical direction is the rotation axis?**
   Look at the rotation axis of the stage.  Theta and ttheta both
   rotate about the vertical axis.  Phi rotates about the longitudinal
   axis (along the beam) when theta = 0.

2. **What is the handedness?**
   Curl the fingers of your right hand around the axis so they point in
   the direction of positive rotation.  If your thumb points along the
   positive physical direction, the rotation is right-handed (`+`).  If
   it points opposite, the rotation is left-handed (`-`).

For this instrument:

| Stage | Physical axis | Handedness | Axis vector |
|-------|---------------|------------|-------------|
| theta | vertical | right-handed | `+VERTICAL` |
| phi | longitudinal | right-handed | `+LONGITUDINAL` |
| ttheta | vertical | right-handed | `+VERTICAL` |

```{important}
Getting the sign wrong inverts the rotation sense.  The symptom is that
`forward()` returns angles with the wrong sign — e.g. theta = -20
instead of +20 for a reflection you know requires positive rotation.
If you see this, negate the axis vector for that stage.
```

---

## Step 3 — Define parent relationships

The `parent` argument on each
{class}`~ad_hoc_diffractometer.stage.Stage` says which stage this one
physically sits on:

- `parent=None` — the stage is mounted on the fixed laboratory frame
  (the floor, the optical table, etc.).
- `parent="theta"` — the stage is mounted on top of the theta stage;
  when theta rotates, this stage rotates with it.

For this instrument:

- **theta** sits on the floor: `parent=None`
- **phi** sits on top of theta: `parent="theta"`
- **ttheta** is the detector arm, mechanically independent of the
  sample stages: `parent=None`

```{important}
Parent relationships encode the mechanical coupling between stages.
If you set `parent=None` on a stage that is actually mounted on
another stage, the forward solver will treat them as independent and
compute incorrect angles.

Conversely, if two stacks share a common base stage (as in the
{func}`~ad_hoc_diffractometer.presets.zaxis` geometry, where the
sample and detector both sit on the alpha stage), then both stacks
must list that shared stage as their parent.
```

---

## Step 4 — Assign roles

Each stage has a `role` string:

- `"sample"` — stages that orient the sample.  Their rotation matrices
  combine to form the sample rotation matrix Z.
- `"detector"` — stages that position the detector.  Their rotation
  matrices combine to form the detector rotation matrix D.

Other roles (e.g. `"analyzer"`, `"polarizer"`) are accepted but are not
used by the forward/inverse solvers.

For this instrument:

| Stage | Role |
|-------|------|
| theta | `"sample"` |
| phi | `"sample"` |
| ttheta | `"detector"` |

```{important}
Assigning the wrong role causes the stage to be placed in the wrong
rotation matrix (Z vs D).  The symptom is that `forward()` produces
solutions that do not satisfy the Bragg condition — `inverse()` of the
returned angles does not give back the original (h, k, l).
```

---

## Step 5 — Create the stages

Putting the decisions together:

```python
stages = [
    Stage("theta",  +VERTICAL,     parent=None,     role="sample"),
    Stage("phi",    +LONGITUDINAL, parent="theta",  role="sample"),
    Stage("ttheta", +VERTICAL,     parent=None,     role="detector"),
]
```

The order of stages in the list does not matter — the package builds the
stacking order from the `parent` references.  By convention, list sample
stages first (floor-most first), then detector stages.

---

## Step 6 — Define diffraction modes

A {class}`~ad_hoc_diffractometer.mode.ConstraintSet` resolves the free
degrees of freedom left after the Bragg condition is satisfied.  An
N-stage geometry has N - 3 free DOF:

- 3-circle: 3 - 3 = 0 free DOF (no constraints needed)
- 4-circle: 4 - 3 = 1 constraint
- 6-circle: 6 - 3 = 3 constraints

This 3-circle geometry has **zero free DOF**, so each mode needs zero
constraints.  However, modes are still useful for documenting the
intended operating condition and for the `computed` field (which lists
the stages the solver writes to):

```python
modes = {
    "theta_ttheta": ahd.ConstraintSet(
        [],
        computed=["theta", "phi", "ttheta"],
    ),
}
```

For geometries with free DOF, you must add constraints.  See
{doc}`constraints` for the full constraint framework.  Common patterns:

- **Bisecting**: use
  {class}`~ad_hoc_diffractometer.mode.BisectConstraint` to tie a sample
  stage to half the detector angle (e.g. `theta = ttheta / 2`).
- **Fixed angle**: use
  {class}`~ad_hoc_diffractometer.mode.SampleConstraint` to hold a stage
  at a declared value (e.g. `phi = 0`).
- **Frozen detector**: use
  {class}`~ad_hoc_diffractometer.mode.DetectorConstraint` to hold a
  detector stage fixed (e.g. `nu = 0` in psic).

```{tip}
The number of constraints must equal N - 3.  Too few and the solver
is under-determined; too many and it is over-determined.  The package
does not enforce this at construction time — it checks at solve time
via
{meth}`~ad_hoc_diffractometer.mode.ConstraintSet.is_fully_constrained`.
```

---

## Step 7 — Assemble the geometry

Pass everything to
{class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`:

```python
g = ahd.AdHocDiffractometer(
    name="powder3c",
    stages=stages,
    basis=basis,
    description="Three-circle powder diffractometer, horizontal scattering plane",
    modes=modes,
    default_mode="theta_ttheta",
)
```

The constructor validates:

- Stage names are unique
- Parent references point to existing stages
- No cycles in the parent graph
- The basis contains exactly 3 mutually orthogonal, non-zero vectors

If any check fails, a `ValueError` is raised with a descriptive message.

---

## Step 8 — Verify the geometry

```python
g.wavelength = 1.5406  # Cu K-alpha, Angstrom
g.sample.lattice = ahd.Lattice(a=5.431)  # cubic silicon
ahd.ub_identity(g.sample)

# Inspect
g.summary()

# Forward calculation: find motor angles for (1, 1, 1)
solutions = g.forward(1, 1, 1)
for sol in solutions:
    print(sol)

# Round-trip check: inverse should recover (1, 1, 1)
for sol in solutions:
    hkl = g.inverse(sol)
    print(f"  -> hkl = ({hkl[0]:.4f}, {hkl[1]:.4f}, {hkl[2]:.4f})")
```

If the round-trip does not recover the original (h, k, l) to within
floating-point precision, check:

1. Are the axis signs correct? (Step 2)
2. Are the parent relationships correct? (Step 3)
3. Are the roles correct? (Step 4)

---

## Optional: register as a factory function

To make the geometry available via
{func}`~ad_hoc_diffractometer.factories.list_geometries` and
{func}`~ad_hoc_diffractometer.factories.make_geometry`, wrap it in a
function decorated with
{func}`~ad_hoc_diffractometer.factories.register_geometry`:

```python
from ad_hoc_diffractometer import register_geometry
from ad_hoc_diffractometer.factories import BASIS_BL

@register_geometry
def powder3c(basis=BASIS_BL):
    VERTICAL     = basis["vertical"]
    LONGITUDINAL = basis["longitudinal"]

    stages = [
        Stage("theta",  +VERTICAL,     parent=None,    role="sample"),
        Stage("phi",    +LONGITUDINAL, parent="theta", role="sample"),
        Stage("ttheta", +VERTICAL,     parent=None,    role="detector"),
    ]
    modes = {
        "theta_ttheta": ahd.ConstraintSet(
            [],
            computed=["theta", "phi", "ttheta"],
        ),
    }
    return ahd.AdHocDiffractometer(
        name="powder3c",
        stages=stages,
        basis=basis,
        description="Three-circle powder diffractometer",
        modes=modes,
        default_mode="theta_ttheta",
    )
```

The function name becomes the geometry name in the registry.  Names must
be unique — duplicates raise `ValueError`.

---

## Advanced: shared base stages

Some geometries have sample and detector stages that share a common base
(e.g. the {func}`~ad_hoc_diffractometer.presets.zaxis` geometry, where
both the sample Z stage and the detector delta stage sit on the alpha
stage).  In this pattern:

```python
stages = [
    Stage("alpha", +VERTICAL,     parent=None,    role="sample"),
    Stage("Z",     +LONGITUDINAL, parent="alpha", role="sample"),
    Stage("delta", -TRANSVERSE,   parent="alpha", role="detector"),
    Stage("gamma", +VERTICAL,     parent="delta", role="detector"),
]
```

The alpha stage is listed with `role="sample"` because it contributes to
the sample rotation matrix Z.  The detector stages (delta, gamma)
reference `parent="alpha"` to express the mechanical coupling: when
alpha rotates, the detector arm moves with it.

This pattern is used for surface diffraction geometries where the
incidence angle is controlled by a shared tilt stage.

---

## Decision checklist

Use this checklist when defining a new geometry:

- [ ] **Basis** — which physical direction maps to which Cartesian axis?
      (See {doc}`basis_vectors`)
- [ ] **Axis vector for each stage** — what physical direction is the
      rotation axis, and is the rotation right-handed (+) or
      left-handed (-)?
- [ ] **Parent chain** — which stage sits on which?  Are sample and
      detector stacks independent or do they share a base stage?
- [ ] **Roles** — is each stage `"sample"` or `"detector"`?
- [ ] **Modes** — how many constraints are needed (N - 3)?  Which
      constraint types resolve the free DOF?
- [ ] **Verification** — does `inverse(forward(h, k, l))` recover
      (h, k, l)?

---

## See also

- {doc}`../quick_start` — tutorial building a fourcv step by step
- {doc}`basis_vectors` — choosing and understanding basis vectors
- {doc}`constraints` — the full constraint framework
- {doc}`modes` — switching between diffraction modes
- {doc}`../problem1` — case study defining the physical reference frame
- {doc}`../problem2` — case study showing that the choice of basis is
  arbitrary (different basis assignments produce U/UB matrices that
  differ by a fixed rotation, leaving the physics invariant)
- {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.stage.Stage`
- {func}`~ad_hoc_diffractometer.factories.register_geometry`
