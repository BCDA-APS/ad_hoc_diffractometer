(howto-modes)=
# Switch Diffraction Modes

A **diffraction mode** describes how `forward()` will compute the motor
angles and which ones remain constant.  See [Concepts](../concepts.md)
for background.

## List available modes

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourcv()
print(list(g.modes.keys()))
# ['bisecting', 'fixed_chi', 'fixed_phi']
```

## Get and set the active mode

```python
# No mode active by default
print(g.mode_name)   # None

# Set a mode
g.mode_name = "bisecting"
print(g.mode_name)   # 'bisecting'
```

## Mode reference

### Bisecting mode

The sample stage angle is constrained to half the detector stage angle,
placing the sample symmetrically between the incident and diffracted beams.
Available on: fourcv (omega = ttheta/2), fourch (omega = ttheta/2),
psic (eta = delta/2), and kappa geometries (komega = ttheta/2 or
komega = delta/2).

```python
g.mode_name = "bisecting"
```

### Fixed chi / fixed phi / fixed kphi / fixed mu

A fixed-angle mode holds one stage at its **current angle** during
`forward()`.  The caller presets the stage angle with `g.set_angle()`
before activating the mode:

```python
# Freeze chi at 90° (surface diffraction geometry)
g.set_angle("chi", 90.0)
g.mode_name = "fixed_chi"
solutions = g.forward(h, k, l)   # chi held at 90° throughout

# Freeze chi at a different angle — same mode, new preset
g.set_angle("chi", 45.0)
solutions = g.forward(h, k, l)   # chi held at 45°
```

```python
# Freeze phi at 0°
g.set_angle("phi", 0.0)
g.mode_name = "fixed_phi"
```

```python
# Freeze mu at 0° on psic or kappa6c
g = ahd.psic()
g.set_angle("mu", 0.0)
g.mode_name = "fixed_mu"
```

### Custom fixed-angle mode

To add a fixed-angle mode that is not pre-built into the geometry:

```python
from ad_hoc_diffractometer import FixedAngleMode

g.modes["my_chi"] = FixedAngleMode(stage="chi", value=90.0)
g.set_angle("chi", 90.0)   # preset the stage angle
g.mode_name = "my_chi"
```

The `value` argument to `FixedAngleMode` sets the **initial** stage angle
as a convenience default.  The solver always reads the stage's current angle
at call time, so a subsequent `g.set_angle()` overrides it.

## Custom modes

Subclass `DiffractionMode` to implement any constraint:

```python
from ad_hoc_diffractometer import DiffractionMode

class MyMode(DiffractionMode):
    @property
    def constrained_stages(self):
        return ["phi"]

    def solve(self, geometry, h, k, l):
        ...  # return list of angle dicts
```

## Clear the active mode

```python
g.mode_name = None   # all stages free
```

## See also

- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`
- {class}`~ad_hoc_diffractometer.mode.ModeDict`
- [Concepts — Diffraction modes](../concepts.md)
