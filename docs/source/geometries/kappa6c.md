(geometry-kappa6c)=
# kappa6c — Kappa Six-Circle

Six-circle kappa diffractometer with psic-style outer axes (mu, nu). The inner sample axes (komega, kappa, and kphi) replace the Eulerian chi circle. Lateral detector, vertical scattering plane.

**Coordinate basis:** You (1999) (`BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.kappa6c()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.kappa6c` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L794) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``mu`` | +vertical (+x) | right-handed |
| ``komega`` | −lateral (−z) | left-handed |
| ``kappa`` | tilted axis, α=50° | right-handed |
| ``kphi`` | −lateral (−z) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``nu`` | +vertical (+x) | right-handed |
| ``delta`` | −lateral (−z) | left-handed |

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`. Preset any fixed-stage angles with `g.set_angle(name, value)` before calling `forward()`. See {doc}`../howto/modes` for usage details and {class}`~ad_hoc_diffractometer.mode.DiffractionMode` for the base class.

### `bisecting`

**Class:** {class}`~ad_hoc_diffractometer.mode.BisectingMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L227))

`komega` is constrained to `delta` / 2, placing the sample symmetrically between the incident and diffracted beams (the bisecting condition).

Additional stages frozen at their current angles: `mu` = current angle, `nu` = current angle.

### `fixed_kphi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `kphi` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("kphi", value)` before calling `forward()`.

### `fixed_mu`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `mu` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("mu", value)` before calling `forward()`.


## API reference

- {func}`~ad_hoc_diffractometer.kappa6c`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
