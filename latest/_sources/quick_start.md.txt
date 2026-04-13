(quick-start)=
# Quick Start

This guide walks through the steps to build a four-circle diffractometer
geometry by hand — without using the pre-built factory function.  Building
it step by step makes the design choices explicit and shows how every
geometry in the package is constructed.

The concise factory-function equivalent appears at the end.

---

## 1. Choose a coordinate basis

The package supports two standard coordinate conventions.  Each maps the
physical directions of the instrument to Cartesian unit vectors:

| Basis | vertical | longitudinal | lateral | Used by |
|---|---|---|---|---|
| `BASIS_YOU` | +x | +y | +z | psic, sixc, kappa6c, zaxis, s2d2, fivec |
| `BASIS_BL` | +z | +y | +x | fourcv, fourch, kappa4cv, kappa4ch |

For a standard four-circle diffractometer (Busing & Levy 1967) choose
`BASIS_BL`:

```python
import numpy as np
import ad_hoc_diffractometer as ahd

BASIS = ahd.BASIS_BL
LATERAL      = BASIS["lateral"]       # +x
LONGITUDINAL = BASIS["longitudinal"]  # +y
VERTICAL     = BASIS["vertical"]      # +z
```

---

## 2. Define the stage stack

A four-circle diffractometer has three sample stages and one detector stage.
Each {class}`~ad_hoc_diffractometer.Stage` is described by:

- **name** — motor name used in angle dictionaries and `forward()` results
- **axis** — rotation axis vector; a leading `−` means left-handed rotation
- **parent** — the stage this one sits on (`None` = directly on the floor)
- **role** — `"sample"` or `"detector"`

The fourcv (vertical scattering plane, synchrotron) stack is:

```python
stages = [
    # Sample stack — base stage first
    ahd.Stage("omega",  -LATERAL,      parent=None,    role="sample"),
    ahd.Stage("chi",    +LONGITUDINAL, parent="omega", role="sample"),
    ahd.Stage("phi",    -LATERAL,      parent="chi",   role="sample"),
    # Detector — independent of the sample stack
    ahd.Stage("ttheta", -LATERAL,      parent=None,    role="detector"),
]
```

`omega` and `ttheta` both rotate about the lateral axis, so their scattering
plane is **vertical** (the synchrotron convention).  For the laboratory
(horizontal scattering plane) convention swap every `LATERAL` for `VERTICAL`
and every `LONGITUDINAL` for `LATERAL` — that is the `fourch` geometry.

---

## 3. Define diffraction modes

A {class}`~ad_hoc_diffractometer.DiffractionMode` specifies which degrees of
freedom are constrained during a forward (hkl → motor angles) calculation.

Two built-in mode types cover most cases:

- {class}`~ad_hoc_diffractometer.BisectingMode` — ties a sample stage to
  half the detector angle, placing the sample symmetrically in the beam.
- {class}`~ad_hoc_diffractometer.FixedAngleMode` — holds one stage at a
  preset angle.

```python
modes = {
    "bisecting": ahd.BisectingMode(
        sample_stage="omega",
        detector_stage="ttheta",
    ),
    "fixed_chi": ahd.FixedAngleMode(stage="chi", value=90.0),
    "fixed_phi": ahd.FixedAngleMode(stage="phi", value=0.0),
}
```

---

## 4. Assemble the geometry

Pass the stage list, basis, and modes to
{class}`~ad_hoc_diffractometer.AdHocDiffractometer`:

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
    omega   axis=-lateral    angle=  0.000°  limits=(-180.0, 180.0)
    chi     axis=+longitudinal  angle=  0.000°  limits=(-180.0, 180.0)
    phi     axis=-lateral    angle=  0.000°  limits=(-180.0, 180.0)

  Detector stages:
    ttheta  axis=-lateral    angle=  0.000°  limits=(-180.0, 180.0)
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

g = ahd.fourcv()               # Busing & Levy (1967) four-circle, vertical plane
g.wavelength = 1.5406          # Å
g.sample.lattice = ahd.Lattice(a=5.431)
print(g.summary())
```

See {doc}`geometries/fourcv` for the full geometry reference, or
{doc}`geometries/fourch` for the horizontal-plane (laboratory) variant.

---

## See also

- {class}`~ad_hoc_diffractometer.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.Stage`
- {class}`~ad_hoc_diffractometer.BisectingMode`
- {class}`~ad_hoc_diffractometer.FixedAngleMode`
- {class}`~ad_hoc_diffractometer.Lattice`
- {doc}`howto/lattice`
- {doc}`howto/orient`
- {doc}`howto/forward`
- {doc}`geometries/fourcv`
- {doc}`geometries/fourch`
