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

Available modes: `bisecting`, `fixed_kphi`

## API reference

- {func}`~ad_hoc_diffractometer.kappa4cv`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
