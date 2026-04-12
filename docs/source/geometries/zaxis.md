(geometry-zaxis)=
# zaxis — Z-Axis Four-Circle (Surface)

Z-axis four-circle diffractometer for surface diffraction. The sample surface normal is parallel to the Z-axis. Sample and detector share an alpha base stage.

**Walko (2016) designation:** (S1D2)1

**Coordinate basis:** You (1999) (`BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.zaxis()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.zaxis` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L872) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``alpha`` | +vertical (+x) | right-handed, shared base |
| ``Z`` | +longitudinal (+y) | right-handed |

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

- {func}`~ad_hoc_diffractometer.zaxis`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- Bloch, *J. Appl. Cryst.* **18**, 33–36 (1985). DOI: [10.1107/S0021889885009858](https://doi.org/10.1107/S0021889885009858)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
