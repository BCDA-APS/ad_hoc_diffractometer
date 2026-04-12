(geometry-fourcv)=
# fourcv — Four-Circle Eulerian (Synchrotron)

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

Available modes: `bisecting`, `fixed_chi`, `fixed_phi`

## API reference

- {func}`~ad_hoc_diffractometer.fourcv`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- Busing & Levy, *Acta Cryst.* **22**, 457–464 (1967). DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
