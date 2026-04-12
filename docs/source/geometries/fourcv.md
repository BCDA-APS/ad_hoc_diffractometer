(geometry-fourcv)=
# fourcv — Eulerian Four-Circle (Synchrotron)

Busing & Levy (1967) four-circle Eulerian diffractometer, vertical scattering plane. ω and 2θ rotate about the lateral axis. Standard synchrotron convention.

**Walko (2016) designation:** S3D1

**Coordinate basis:** Busing & Levy (`BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourcv()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.fourcv` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L502) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``omega`` | −lateral (−x BL) | left-handed |
| ``chi`` | +longitudinal (+y BL) | right-handed |
| ``phi`` | −lateral (−x BL) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``ttheta`` | −lateral (−x BL) | left-handed |

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`. Preset any fixed-stage angles with `g.set_angle(name, value)` before calling `forward()`. See {doc}`../howto/modes` for usage details and {class}`~ad_hoc_diffractometer.mode.DiffractionMode` for the base class.

### `bisecting`

**Class:** {class}`~ad_hoc_diffractometer.mode.BisectingMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L227))

`omega` is constrained to `ttheta` / 2, placing the sample symmetrically between the incident and diffracted beams (the bisecting condition).

### `fixed_chi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `chi` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("chi", value)` before calling `forward()`.

### `fixed_phi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `phi` fixed at its **current angle** during `forward()`. Preset the angle with `g.set_angle("phi", value)` before calling `forward()`.


## API reference

- {func}`~ad_hoc_diffractometer.fourcv`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- Busing & Levy, *Acta Cryst.* **22**, 457–464 (1967). DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
