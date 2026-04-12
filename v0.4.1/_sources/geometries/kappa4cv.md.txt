(geometry-kappa4cv)=
# kappa4cv — Kappa Four-Circle (Synchrotron)

Four-circle kappa diffractometer, vertical scattering plane. The chi circle is replaced by a kappa axis tilted at α = 50° from the vertical toward the lateral axis.

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy (`BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.kappa4cv()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.kappa4cv` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L668) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``komega`` | −lateral (−x BL) | left-handed |
| ``kappa`` | tilted axis, α=50° | right-handed |
| ``kphi`` | −lateral (−x BL) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``ttheta`` | −lateral (−x BL) | left-handed |

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`. Preset any fixed-stage angles with `g.set_angle(name, value)` before calling `forward()`. See {doc}`../howto/modes` for usage details and {class}`~ad_hoc_diffractometer.mode.DiffractionMode` for the base class.

### `bisecting`

**Class:** {class}`~ad_hoc_diffractometer.mode.BisectingMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L227))

`komega` is constrained to `ttheta` / 2, placing the sample symmetrically between the incident and diffracted beams (the bisecting condition).

### `fixed_kphi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `kphi` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("kphi", value)` before calling `forward()`.


## API reference

- {func}`~ad_hoc_diffractometer.kappa4cv`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
