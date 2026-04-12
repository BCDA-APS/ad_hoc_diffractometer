(geometry-psic)=
# psic — Eulerian Six-Circle, 4S+2D (You 1999)

You (1999) 4S+2D six-circle diffractometer. Four sample stages (mu, eta, chi, and phi) and two detector stages (nu, delta). Lateral detector, vertical scattering plane. Standard synchrotron six-circle.

**Coordinate basis:** You (1999) (`BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.psic()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.psic` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L439) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``mu`` | +vertical (+x) | right-handed |
| ``eta`` | −lateral (−z) | left-handed |
| ``chi`` | +longitudinal (+y) | right-handed |
| ``phi`` | −lateral (−z) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``nu`` | +vertical (+x) | right-handed |
| ``delta`` | −lateral (−z) | left-handed |

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`. Preset any fixed-stage angles with `g.set_angle(name, value)` before calling `forward()`. See {doc}`../howto/modes` for usage details and {class}`~ad_hoc_diffractometer.mode.DiffractionMode` for the base class.

### `bisecting`

**Class:** {class}`~ad_hoc_diffractometer.mode.BisectingMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L227))

`eta` is constrained to `delta` / 2, placing the sample symmetrically between the incident and diffracted beams (the bisecting condition).

Additional stages frozen at their current angles: `mu` = current angle, `nu` = current angle.

### `fixed_chi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `chi` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("chi", value)` before calling `forward()`.

### `fixed_phi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `phi` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("phi", value)` before calling `forward()`.

### `fixed_mu`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `mu` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("mu", value)` before calling `forward()`.


## API reference

- {func}`~ad_hoc_diffractometer.psic`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- You, *J. Appl. Cryst.* **32**, 614–623 (1999). DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
