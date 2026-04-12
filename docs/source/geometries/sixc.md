(geometry-sixc)=
# sixc — Eulerian Six-Circle, Surface (Lohmeier & Vlieg 1993)

Six-circle surface diffractometer. Sample and detector share a common alpha (rotary table) base stage. Designed for surface diffraction.

**Walko (2016) designation:** (S3D2)1

**Coordinate basis:** You (1999) (`BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.sixc()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.sixc` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L609) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``alpha`` | +vertical (+x) | right-handed, shared base |
| ``omega`` | −lateral (−z) | left-handed |
| ``chi`` | +longitudinal (+y) | right-handed |
| ``phi`` | −lateral (−z) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``delta`` | −lateral (−z) | left-handed |
| ``gamma`` | +vertical (+x) | right-handed |

**Shared stage:** alpha (base stage shared between sample and detector stacks)

## Diffraction modes

*(No modes defined for this geometry.  All stages are free during*
*`forward()` — the solver explores the full solution space.)*

## API reference

- {func}`~ad_hoc_diffractometer.sixc`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- Lohmeier & Vlieg, *J. Appl. Cryst.* **26**, 706–716 (1993). DOI: [10.1107/S0021889893006198](https://doi.org/10.1107/S0021889893006198)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
