(geometry-fourch)=
# fourch — Eulerian Four-Circle (Laboratory)

Busing & Levy (1967) four-circle Eulerian diffractometer, horizontal scattering plane. ω and 2θ rotate about the vertical axis. Standard laboratory convention.

**Walko (2016) designation:** S3D1

**Coordinate basis:** Busing & Levy (`BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourch()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.fourch` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L556) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``omega`` | −vertical (−z BL) | left-handed |
| ``chi`` | +longitudinal (+y BL) | right-handed |
| ``phi`` | −vertical (−z BL) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``ttheta`` | −vertical (−z BL) | left-handed |

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`. See {doc}`../howto/modes` for usage details and {class}`~ad_hoc_diffractometer.mode.DiffractionMode` for the base class.

### `bisecting`

**Class:** {class}`~ad_hoc_diffractometer.mode.BisectingMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L215))

The bisecting condition: the sample stage angle equals half the detector angle, placing the sample symmetrically between the incident and diffracted beams.

### `fixed_chi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `chi` fixed at 90.0° during `forward()`. All other free stages are computed by the solver.

### `fixed_phi`

**Class:** {class}`~ad_hoc_diffractometer.mode.FixedAngleMode` ([source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/mode.py#L156))

Holds `phi` fixed at 0.0° during `forward()`. All other free stages are computed by the solver.


## API reference

- {func}`~ad_hoc_diffractometer.fourch`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- Busing & Levy, *Acta Cryst.* **22**, 457–464 (1967). DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
