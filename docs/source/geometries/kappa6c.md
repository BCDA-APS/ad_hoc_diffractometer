(geometry-kappa6c)=
# kappa6c — Six-Circle Kappa (Synchrotron)

Six-circle kappa diffractometer with psic-style outer axes (mu, nu). The inner sample axes (komega, kappa, and kphi) replace the Eulerian chi circle. Lateral detector, vertical scattering plane.

**Coordinate basis:** You (1999) (`BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.kappa6c()
g.wavelength = 1.0  # Å
print(g.summary())
```

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

Available modes: `bisecting`, `fixed_kphi`, `fixed_mu`

## API reference

- {func}`~ad_hoc_diffractometer.kappa6c`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
