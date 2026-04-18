(geometry-psic)=
# psic — Eulerian Six-Circle, 4S+2D (You 1999)

You (1999) 4S+2D six-circle diffractometer. Four sample stages (mu, eta, chi, and phi) and two detector stages (nu, delta). Lateral detector, vertical scattering plane. Standard synchrotron six-circle.

**Coordinate basis:** You (1999) ({data}`~ad_hoc_diffractometer.factories.BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.psic()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.factories.psic` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L454) for the complete stage
and mode configuration.

## Stage layout

**Sample stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``mu`` | +vertical (+x) | right-handed |
| ``eta`` | −lateral (−z) | left-handed |
| ``chi`` | +longitudinal (+y) | right-handed |
| ``phi`` | −lateral (−z) | left-handed |

**Detector stages (base first):**

| Stage | Axis | Handedness |
|---|---|---|
| ``nu`` | +vertical (+x) | right-handed |
| ``delta`` | −lateral (−z) | left-handed |

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`.
Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 3 constraints
(N − 3 = 3 for N = 6 DOF).
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

**Bisect pairs:**

- Vertical plane: eta (lateral) ↔ delta (lateral) → `eta = delta/2`
- Horizontal plane: mu (vertical) ↔ nu (vertical) → `mu = nu/2`

### `bisecting_vertical` *(default)*

{class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`:
`eta = delta/2`, `mu = 0`, `nu = 0`.
Vertical scattering plane bisecting condition (You 1999, §5.3).

| | |
|---|---|
| **Computed** | eta, chi, phi, delta |
| **Constant during** `forward()` | mu = 0, nu = 0 |

### `fixed_chi`

`chi` held at declared value (default 90°), `eta = delta/2`, `nu = 0`.
The caller chooses the chi value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` — see {doc}`../howto/constraints`.

| | |
|---|---|
| **Computed** | eta, phi, delta |
| **Constant during** `forward()` | chi, mu = 0, nu = 0 |

### `fixed_phi`

`phi` held at declared value (default 0°), `eta = delta/2`, `nu = 0`.

| | |
|---|---|
| **Computed** | eta, chi, delta |
| **Constant during** `forward()` | phi, mu = 0, nu = 0 |

### `fixed_mu`

`mu` held at declared value (default 0°), `eta = delta/2`, `nu = 0`.

| | |
|---|---|
| **Computed** | eta, chi, phi, delta |
| **Constant during** `forward()` | mu, nu = 0 |

### `bisecting_horizontal`

{class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`:
`mu = nu/2`, `eta = 0`, `delta = 0`.
Horizontal scattering plane bisecting condition (You 1999, §5.1).

| | |
|---|---|
| **Computed** | mu, chi, phi, nu |
| **Constant during** `forward()` | eta = 0, delta = 0 |

### `fixed_nu`

`nu` held at declared value (default 0°), `eta = delta/2`, `mu = 0`.

| | |
|---|---|
| **Computed** | eta, chi, phi, delta |
| **Constant during** `forward()` | nu, mu = 0 |

### `double_diffraction_vertical`

{class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`
with secondary-reflection extras.
Bisecting solver runs; simultaneous diffraction logic not yet implemented.

| | |
|---|---|
| **Computed** | eta, chi, phi, delta |
| **Constant during** `forward()` | mu = 0, nu = 0 |
| **Extras (input)** | h₂, k₂, l₂ (secondary reflection Miller indices) |

### `lifting_detector_mu` *(stub)*

Out-of-plane mode: mu and eta frozen, nu and delta free.
Requires the qaz pseudo-angle solver — not yet implemented.

| | |
|---|---|
| **Computed** | mu, nu, delta |
| **Constant during** `forward()` | mu, eta |

### `lifting_detector_phi` *(stub)*

Out-of-plane mode: phi and mu frozen, nu and delta free.
Requires the qaz pseudo-angle solver — not yet implemented.

| | |
|---|---|
| **Computed** | phi, nu, delta |
| **Constant during** `forward()` | phi, mu |

### `psi_constant_vertical` *(stub)*

Vertical bisecting with azimuthal angle ψ of n̂ about Q fixed.
Requires reference vector infrastructure (Issue J / #157).

| | |
|---|---|
| **Computed** | eta, chi, phi, delta |
| **Constant during** `forward()` | mu = 0, nu = 0 |
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

### `psi_constant_horizontal` *(stub)*

Horizontal bisecting with azimuthal angle ψ fixed.
Requires reference vector infrastructure (Issue J / #157).

| | |
|---|---|
| **Computed** | mu, chi, phi, nu |
| **Constant during** `forward()` | eta = 0, delta = 0 |
| **Extras (input)** | n̂, ψ |
| **Extras (output)** | psi |

## API reference

- {func}`~ad_hoc_diffractometer.factories.psic`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {class}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {class}`~ad_hoc_diffractometer.mode.ConstraintViolation`

## References

- You, *J. Appl. Cryst.* **32**, 614–623 (1999). DOI: [10.1107/S0021889899001223](https://doi.org/10.1107/S0021889899001223)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
