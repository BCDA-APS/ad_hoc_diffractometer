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

Constrains omega = ttheta/2 (sample bisects the scattered beam).
Chi and phi are free.  Available on: fourcv, fourch, psic, kappa geometries.

```python
g.mode_name = "bisecting"
```

### Fixed chi

Holds chi at a specified angle.  Useful for surface diffraction or
when the goniometer range is limited.

```python
from ad_hoc_diffractometer import FixedAngleMode
g.modes["my_chi"] = FixedAngleMode(stage="chi", value=90.0)
g.mode_name = "my_chi"
```

### Fixed phi

Holds phi at a specified angle.

```python
g.mode_name = "fixed_phi"
# Default value is 0.0; override at construction or via:
g.modes["fixed_phi"] = FixedAngleMode(stage="phi", value=45.0)
```

### Fixed mu (psic / kappa6c)

Holds the outer mu stage fixed.

```python
g = ahd.psic()
g.mode_name = "fixed_mu"
```

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
