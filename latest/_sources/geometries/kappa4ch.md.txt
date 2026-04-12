(geometry-kappa4ch)=
# kappa4ch — Kappa Four-Circle (Laboratory)

Four-circle kappa diffractometer, horizontal scattering plane. Kappa axis tilted at α = 50° from vertical. Laboratory convention.

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy (`BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.kappa4ch()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``komega`` | −vertical (−z BL) | left-handed |
| ``kappa`` | tilted axis, α=50° | right-handed |
| ``kphi`` | −vertical (−z BL) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``ttheta`` | −vertical (−z BL) | left-handed |

## Diffraction modes

Available modes: `bisecting`, `fixed_kphi`

## API reference

- {func}`~ad_hoc_diffractometer.kappa4ch`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
