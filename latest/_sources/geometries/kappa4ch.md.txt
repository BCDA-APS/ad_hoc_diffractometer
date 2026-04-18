(geometry-kappa4ch)=
# kappa4ch — Kappa Four-Circle (Laboratory)

Four-circle kappa diffractometer, horizontal scattering plane. Kappa axis tilted at α = 50° from vertical. Laboratory convention.

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy ({data}`~ad_hoc_diffractometer.factories.BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.kappa4ch()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.factories.kappa4ch` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L964) for the complete stage
and mode configuration.

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

**Virtual Eulerian angles** (computed from real kappa angles via Walko 2016 eq. [16]):
omega, chi, phi.  Used as constraint names for stub modes; converted back to
komega, kappa, kphi by the kappa inversion solver (Issue I / #153).

## Diffraction modes

Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 1 constraint
(N − 3 = 1 for N = 4 DOF).
Identical mode set to {doc}`kappa4cv`.
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

### `bisecting` *(default)*

{class}`~ad_hoc_diffractometer.mode.BisectConstraint`:
`komega = ttheta / 2` (approximates the bisecting condition).

> **Note:** The correct bisecting condition is virtual `omega_euler = 0`.
> Corrected in Issue I / #153.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | — |

### `fixed_kphi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
`kphi` held at declared value (default 0°) — real stage, no kappa inversion needed.

| | |
|---|---|
| **Computed** | komega, kappa, ttheta |
| **Constant during** `forward()` | kphi |

### `constant_omega` *(stub)*

Fix virtual Eulerian omega (default 0°). Requires Issue I / #153.

### `constant_chi` *(stub)*

Fix virtual Eulerian chi (default 90°). Requires Issue I / #153.

### `constant_phi` *(stub)*

Fix virtual Eulerian phi (default 0°). Requires Issue I / #153.

### `psi_constant` *(stub)*

Fix azimuthal angle ψ of n̂ about Q.
Requires kappa inversion (Issue I) and reference infrastructure (Issue J / #157).

| | |
|---|---|
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

## API reference

- {func}`~ad_hoc_diffractometer.factories.kappa4ch`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {class}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {class}`~ad_hoc_diffractometer.mode.ConstraintViolation`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), eq. [16].
