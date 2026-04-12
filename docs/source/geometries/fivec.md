(geometry-fivec)=
# fivec — Eulerian Five-Circle (Vlieg et al. 1987)

Five-circle diffractometer: a standard fourcv (Eulerian four-circle) mounted on a vertical mu base stage. Sample and detector are coupled through mu.

**Walko (2016) designation:** (S3D1)1

**Coordinate basis:** You (1999) (`BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.fivec()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.fivec` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L974) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``mu`` | +vertical (+x) | right-handed, shared base |
| ``omega`` | −lateral (−z) | left-handed |
| ``chi`` | +longitudinal (+y) | right-handed |
| ``phi`` | −lateral (−z) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``ttheta`` | −lateral (−z) | left-handed |

**Shared stage:** mu (base stage shared between sample and detector stacks)

## Diffraction modes

*(No modes defined for this geometry.  All stages are free during*
*`forward()` — the solver explores the full solution space.)*

## API reference

- {func}`~ad_hoc_diffractometer.fivec`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.DiffractionMode`
- {class}`~ad_hoc_diffractometer.mode.BisectingMode`
- {class}`~ad_hoc_diffractometer.mode.FixedAngleMode`

## References

- Vlieg et al., *J. Appl. Cryst.* **20**, 330–337 (1987). DOI: [10.1107/S0021889887087266](https://doi.org/10.1107/S0021889887087266)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
