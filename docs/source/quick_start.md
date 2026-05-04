(quick-start)=
# Quick Start

This guide walks through the steps to build a four-circle diffractometer
by hand — without using a demo geometry.  Building it step by step makes
the design choices explicit and shows how every diffractometer in the
package is constructed.  This builds a diffractometer for the vertical
scattering plane.

The concise demo geometry equivalent appears at the end.

```{tip}
If you already have the package installed and just want a five-line
example that uses one of the demo geometries, see the
{ref}`Demonstration Geometries Quick start <quick-start-demo>` page
instead.
```

---

## 1. Choose a coordinate basis

Any right-handed orthogonal mapping of the three physical directions (vertical,
longitudinal, transverse) to Cartesian unit vectors (x, y, z) is a valid basis.
Two conventions (shown below) are used in the demo examples. You may supply
**any basis** (dictionary) that satisfies those constraints:

| Basis | vertical | longitudinal | transverse | Used by demos |
|---|---|---|---|---|
| {data}`~ad_hoc_diffractometer.factories.BASIS_YOU` | +x | +y | +z | {func}`~ad_hoc_diffractometer.presets.psic`, {func}`~ad_hoc_diffractometer.presets.sixc`, {func}`~ad_hoc_diffractometer.presets.kappa6c`, {func}`~ad_hoc_diffractometer.presets.zaxis`, {func}`~ad_hoc_diffractometer.presets.s2d2`, {func}`~ad_hoc_diffractometer.presets.fivec` |
| {data}`~ad_hoc_diffractometer.factories.BASIS_BL` | +z | +y | +x | {func}`~ad_hoc_diffractometer.presets.fourcv`, {func}`~ad_hoc_diffractometer.presets.fourch`, {func}`~ad_hoc_diffractometer.presets.kappa4cv`, {func}`~ad_hoc_diffractometer.presets.kappa4ch` |

See {doc}`howto/basis_vectors` for a full tutorial on choosing and
understanding basis vectors.

For a standard four-circle diffractometer (Busing & Levy 1967) choose
{data}`~ad_hoc_diffractometer.factories.BASIS_BL`:

```python
import numpy as np
import ad_hoc_diffractometer as ahd

BASIS = ahd.BASIS_BL
TRANSVERSE      = BASIS["transverse"]       # +x
LONGITUDINAL = BASIS["longitudinal"]  # +y
VERTICAL     = BASIS["vertical"]      # +z
```

---

## 2. Define the stage stack

A four-circle diffractometer has three sample stages and one detector stage.
Each {class}`~ad_hoc_diffractometer.stage.Stage` is described by:

- **name** — motor name used in angle dictionaries and `forward()` results
- **axis** — rotation axis vector; a leading `−` means left-handed rotation
- **parent** — the stage this one sits on (`None` = directly on the floor)
- **role** — `"sample"` or `"detector"`

The fourcv (vertical scattering plane, synchrotron) stack is:

```python
stages = [
    # Sample stack — base stage first
    ahd.Stage("omega",  -TRANSVERSE,      parent=None,    role="sample"),
    ahd.Stage("chi",    +LONGITUDINAL, parent="omega", role="sample"),
    ahd.Stage("phi",    -TRANSVERSE,      parent="chi",   role="sample"),
    # Detector — independent of the sample stack
    ahd.Stage("ttheta", -TRANSVERSE,      parent=None,    role="detector"),
]
```

`omega` and `ttheta` both rotate about the transverse axis, so their scattering
plane is **vertical** (the synchrotron convention).  For the laboratory
(horizontal scattering plane) convention swap every `TRANSVERSE` for `VERTICAL`
and every `LONGITUDINAL` for `TRANSVERSE` — that is the {func}`~ad_hoc_diffractometer.presets.fourch` geometry.

---

## 3. Define diffraction modes

A {class}`~ad_hoc_diffractometer.mode.ConstraintSet` specifies which degrees of
freedom are constrained during a forward (hkl → motor angles) calculation.

Three constraint types cover most cases:

- {class}`~ad_hoc_diffractometer.mode.BisectConstraint` — ties a sample stage to
  half the detector angle, placing the sample symmetrically in the beam.
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint` — holds one sample stage
  at a fixed angle.
- {class}`~ad_hoc_diffractometer.mode.DetectorConstraint` — holds one detector
  stage at a fixed angle.

```python
modes = {
    "bisecting": ahd.ConstraintSet(
        [ahd.BisectConstraint("omega", "ttheta")],
        computed=["omega", "chi", "phi", "ttheta"],
    ),
    "fixed_chi": ahd.ConstraintSet(
        [ahd.SampleConstraint("chi", 90.0)],
        computed=["omega", "phi", "ttheta"],
    ),
    "fixed_phi": ahd.ConstraintSet(
        [ahd.SampleConstraint("phi", 0.0)],
        computed=["omega", "chi", "ttheta"],
    ),
}
```

---

## 4. Assemble the geometry

Pass the stage list, basis, and modes to
{class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`:

```python
g = ahd.AdHocDiffractometer(
    name="my_fourcv",
    stages=stages,
    basis=BASIS,
    description="Four-circle Eulerian, vertical scattering plane",
    modes=modes,
    default_mode="bisecting",
)
```

---

## 5. Set the wavelength and sample lattice

```python
g.wavelength = 1.5406          # Å (Cu Kα)
g.sample.lattice = ahd.Lattice(a=5.431)  # cubic silicon
```

---

## 6. Inspect the geometry

```python
print(g.summary())
```

Example output:

```
Geometry: my_fourcv
  Four-circle Eulerian, vertical scattering plane
  Wavelength: 1.5406 Å   Energy: 8.0478 keV
  Mode: bisecting

  Sample stages:
    omega   axis=-transverse    angle=  0.000°  limits=(-180.0, 180.0)
    chi     axis=+longitudinal  angle=  0.000°  limits=(-180.0, 180.0)
    phi     axis=-transverse    angle=  0.000°  limits=(-180.0, 180.0)

  Detector stages:
    ttheta  axis=-transverse    angle=  0.000°  limits=(-180.0, 180.0)
```

---

## 7. Solve the forward problem

Given a reflection (hkl), find the motor angles that satisfy Bragg's law.
First orient the crystal (see {doc}`howto/orient`), then call `forward()`:

```python
# Minimal orientation: identity U matrix (crystal axes || diffractometer axes)
ahd.ub_identity(g.sample)

# Solve for the (0, 0, 4) reflection
solutions = g.forward(0, 0, 4)
for sol in solutions:
    print(sol)
```

---

## Concise form — the factory function

The code above is exactly what the built-in `fourcv()` factory does.
If you do not need to customise anything, use it directly:

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.fourcv()       # Busing & Levy (1967) four-circle, vertical plane
g.wavelength = 1.5406          # Å
g.sample.lattice = ahd.Lattice(a=5.431)
print(g.summary())
```

See {doc}`geometries/fourcv` for the full geometry reference, or
{doc}`geometries/fourch` for the horizontal-plane (laboratory) variant.

---

## See also

- {ref}`Demonstration Geometries Quick start <quick-start-demo>` — the
  five-line demo geometry example
- {doc}`howto/custom_geometry` — how to build a diffractometer that is
  not one of the demos
- {doc}`howto/basis_vectors` — choosing and understanding basis vectors
- {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.stage.Stage`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.lattice.Lattice`
- {doc}`howto/lattice`
- {doc}`howto/orient`
- {doc}`howto/forward`
- {doc}`geometries/fourcv`
- {doc}`geometries/fourch`
