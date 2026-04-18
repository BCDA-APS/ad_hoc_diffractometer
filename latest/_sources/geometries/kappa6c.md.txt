(geometry-kappa6c)=
# kappa6c — Kappa Six-Circle

Six-circle kappa diffractometer with psic-style outer axes (mu, nu). The inner sample axes (komega, kappa, and kphi) replace the Eulerian chi circle. Lateral detector, vertical scattering plane.

**Coordinate basis:** You (1999) ({data}`~ad_hoc_diffractometer.factories.BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.kappa6c()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.presets.kappa6c` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L1049) for the complete stage
and mode configuration.

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

**Virtual Eulerian angles** (computed from komega, kappa, kphi via Walko 2016 eq. [16]):
omega, chi, phi.  Used in stub mode descriptions; kappa inversion required (Issue I / #153).

**Bisect pairs:**

- Vertical: komega (lateral) ↔ delta (lateral) → `komega = delta/2`
- Horizontal: mu (vertical) ↔ nu (vertical) → `mu = nu/2`

## Diffraction modes

Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 3 constraints
(N − 3 = 3 for N = 6 DOF).
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

### `bisecting_vertical` *(default)*

{class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`:
`komega = delta/2`, `mu = 0`, `nu = 0`.
Vertical scattering plane (psic-style).

| | |
|---|---|
| **Computed** | komega, kappa, kphi, delta |
| **Constant during** `forward()` | mu = 0, nu = 0 |

### `bisecting_horizontal`

{class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`:
`mu = nu/2`, `komega = 0`, `delta = 0`.
Horizontal scattering plane.

| | |
|---|---|
| **Computed** | mu, kappa, kphi, nu |
| **Constant during** `forward()` | komega = 0, delta = 0 |

### `fixed_kphi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
`kphi` held at declared value (default 0°), `mu = 0`, `nu = 0`.

| | |
|---|---|
| **Computed** | komega, kappa, delta |
| **Constant during** `forward()` | kphi, mu = 0, nu = 0 |

### `fixed_mu`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`:
`mu` held at declared value (default 0°), `komega = delta/2`, `nu = 0`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, delta |
| **Constant during** `forward()` | mu, nu = 0 |

### `fixed_nu`

{class}`~ad_hoc_diffractometer.mode.DetectorConstraint` + {class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
`nu` held at declared value (default 0°), `komega = delta/2`, `mu = 0`.
Analogous to psic `fixed_nu`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, delta |
| **Constant during** `forward()` | nu, mu = 0 |

### `fixed_delta`

{class}`~ad_hoc_diffractometer.mode.DetectorConstraint` + {class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
`delta` held at declared value (default 0°), `mu = nu/2`, `komega = 0`.
Horizontal plane with delta frozen.

| | |
|---|---|
| **Computed** | mu, kappa, kphi, nu |
| **Constant during** `forward()` | delta, komega = 0 |

### `lifting_detector_mu`

Out-of-plane mode: mu and komega frozen, nu and delta solved via the qaz
constraint (``tan(qaz) = tan(delta) / sin(nu)``, You 1999 eq. 18).
``qaz = 90°`` constrains the scattering to the vertical plane.

| | |
|---|---|
| **Computed** | mu, nu, delta |
| **Constant during** `forward()` | mu = 0, komega = 0 |

### `lifting_detector_kphi`

Out-of-plane mode: kphi and mu frozen, nu and delta solved via the qaz
constraint (``tan(qaz) = tan(delta) / sin(nu)``, You 1999 eq. 18).
``qaz = 90°`` constrains the scattering to the vertical plane.

| | |
|---|---|
| **Computed** | kphi, nu, delta |
| **Constant during** `forward()` | kphi = 0, mu = 0 |

### `psi_constant_vertical`

Vertical bisecting with azimuthal angle ψ validation.
Set ``g.azimuthal_reference = (h, k, l)`` before calling ``forward()``.
The solver returns bisecting solutions only when the natural ψ for the
requested (h,k,l) matches the stored target.  See {doc}`../howto/surface`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, delta |
| **Constant during** `forward()` | mu = 0, nu = 0 |
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

### `psi_constant_horizontal`

Horizontal bisecting with azimuthal angle ψ validation.
Symmetric with `psi_constant_vertical` in the horizontal plane.
Set ``g.azimuthal_reference = (h, k, l)`` before calling ``forward()``.

| | |
|---|---|
| **Computed** | mu, kappa, kphi, nu |
| **Constant during** `forward()` | komega = 0, delta = 0 |
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

### `double_diffraction_vertical`

Full 4D simultaneous solver in the vertical scattering plane: finds motor
angles where both the primary (h₁,k₁,l₁) and secondary (h₂,k₂,l₂)
reflections satisfy the Ewald sphere condition.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, delta |
| **Constant during** `forward()` | mu = 0, nu = 0 |
| **Extras (input)** | h₂, k₂, l₂ (secondary reflection Miller indices) |

### `double_diffraction_horizontal`

Full 4D simultaneous solver in the horizontal scattering plane.

| | |
|---|---|
| **Computed** | mu, kappa, kphi, nu |
| **Constant during** `forward()` | komega = 0, delta = 0 |
| **Extras (input)** | h₂, k₂, l₂ (secondary reflection Miller indices) |

## API reference

- {func}`~ad_hoc_diffractometer.presets.kappa6c`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {class}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {class}`~ad_hoc_diffractometer.mode.ConstraintViolation`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- You, *J. Appl. Cryst.* **32**, 614–623 (1999). DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), eq. [16].
