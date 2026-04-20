(geometry-sixc)=
# sixc — Eulerian Six-Circle, Surface (Lohmeier & Vlieg 1993)

Six-circle surface diffractometer. Sample and detector share a common alpha (rotary table) base stage. Supports both bulk crystallography (four-circle mode) and surface diffraction.

**Walko (2016) designation:** (S3D2)1

**Coordinate basis:** You (1999) ({data}`~ad_hoc_diffractometer.factories.BASIS_YOU`): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.sixc()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.presets.sixc` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L689) for the complete stage
and mode configuration.

## Stage layout

### Stage coupling

```{graphviz}
digraph sixc {
    rankdir=BT;
    label="sixc";
    labelloc=t;
    fontsize=14;
    node [shape=box, style=filled, fontsize=11];

    alpha [label="alpha\naxis: +vertical\nRH", fillcolor="#a8d8ea"];
    omega [label="omega\naxis: -lateral\nLH", fillcolor="#a8d8ea"];
    chi [label="chi\naxis: +longitudinal\nRH", fillcolor="#a8d8ea"];
    phi [label="phi\naxis: -lateral\nLH", fillcolor="#a8d8ea"];
    delta [label="delta\naxis: -lateral\nLH", fillcolor="#f8a5a5"];
    gamma [label="gamma\naxis: +vertical\nRH", fillcolor="#f8a5a5"];

    omega -> alpha;
    chi -> omega;
    phi -> chi;
    delta -> alpha;
    gamma -> delta;

    // Legend
    subgraph cluster_legend {
        label="Legend";
        fontsize=8;
        style=dashed;
        color=gray;
        sample_legend [label="sample", fillcolor="#a8d8ea", shape=box, style=filled, fontsize=7];
        detector_legend [label="detector", fillcolor="#f8a5a5", shape=box, style=filled, fontsize=7];
        sample_legend -> detector_legend [style=invis];
    }
}
```

### Axis overview

![sixc stage axes](../_static/geometries/sixc/sixc_all.svg)

### Per-stage axis diagrams

::::{tab-set}
:::{tab-item} alpha
![alpha axis](../_static/geometries/sixc/alpha.svg)
:::
:::{tab-item} omega
![omega axis](../_static/geometries/sixc/omega.svg)
:::
:::{tab-item} chi
![chi axis](../_static/geometries/sixc/chi.svg)
:::
:::{tab-item} phi
![phi axis](../_static/geometries/sixc/phi.svg)
:::
:::{tab-item} delta
![delta axis](../_static/geometries/sixc/delta.svg)
:::
:::{tab-item} gamma
![gamma axis](../_static/geometries/sixc/gamma.svg)
:::
::::

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

**Shared stage:** alpha (rotary table base shared between sample and detector stacks)

## Diffraction modes

Set the active mode with `g.mode_name = "<mode>"`.
Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 3 constraints
(N − 3 = 3 for N = 6 DOF).
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

### `bisecting_4c` *(default)*

{class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint` + {class}`~ad_hoc_diffractometer.mode.BisectConstraint`:
`alpha = 0`, `gamma = 0`, `omega = delta / 2`.
Reduces to standard four-circle bisecting geometry.

| | |
|---|---|
| **Computed** | omega, chi, phi, delta |
| **Constant during** `forward()` | alpha = 0, gamma = 0 |

### `fixed_gamma_5c`

{class}`~ad_hoc_diffractometer.mode.DetectorConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.BisectConstraint`:
`alpha = 0`, `omega = delta / 2`.
`gamma` is held at the value declared in the constraint (factory default: 0°).
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`; the constraint
persists until replaced — see {doc}`../howto/constraints`.

| | |
|---|---|
| **Computed** | omega, chi, phi, delta, alpha |
| **Constant during** `forward()` | gamma, alpha = 0 |

### `fixed_alpha_5c`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.BisectConstraint` + {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`:
`omega = delta / 2`, `gamma = 0`.
`alpha` is held at the value declared in the constraint (factory default: 0°).
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | omega, chi, phi, delta, gamma |
| **Constant during** `forward()` | alpha, gamma = 0 |

### `fixed_alpha_zaxis`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint` × 2 + {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`:
Z-axis mode with fixed incidence angle. Requires ``g.surface_normal = (h, k, l)`` — see {doc}`../howto/surface`.

| | |
|---|---|
| **Computed** | omega, delta, gamma |
| **Constant during** `forward()` | alpha (= β_in), chi, phi |
| **Extras (input)** | n̂ (surface normal) |
| **Extras (output)** | alpha_i (incidence angle), beta_out (exit angle) |

### `fixed_beta_zaxis`

{class}`~ad_hoc_diffractometer.mode.DetectorConstraint` + {class}`~ad_hoc_diffractometer.mode.SampleConstraint` + {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`:
Z-axis mode with fixed exit angle. Requires ``g.surface_normal = (h, k, l)`` — see {doc}`../howto/surface`.

| | |
|---|---|
| **Computed** | omega, delta, alpha |
| **Constant during** `forward()` | gamma (= β_out), chi |
| **Extras (input)** | n̂ |
| **Extras (output)** | alpha_i, beta_out |

### `alpha_eq_beta_zaxis`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint` × 2 + {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`:
Z-axis mode, symmetric reflection (α = γ, β_in = β_out). Requires ``g.surface_normal = (h, k, l)`` — see {doc}`../howto/surface`.

| | |
|---|---|
| **Computed** | omega, delta, alpha, gamma |
| **Constant during** `forward()` | chi, phi |
| **Extras (input)** | n̂ |
| **Extras (output)** | alpha_i, beta_out |

## API reference

- {func}`~ad_hoc_diffractometer.presets.sixc`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {exc}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {exc}`~ad_hoc_diffractometer.mode.ConstraintViolation`

## References

- Lohmeier & Vlieg, *J. Appl. Cryst.* **26**, 706–716 (1993). DOI: [10.1107/S0021889893006198](https://doi.org/10.1107/S0021889893006198)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
