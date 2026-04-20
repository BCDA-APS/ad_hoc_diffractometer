(geometry-kappa4cv)=
# kappa4cv — Kappa Four-Circle (Synchrotron)

Four-circle kappa diffractometer, vertical scattering plane. The chi circle is replaced by a kappa axis tilted at α = 50° from the vertical toward the lateral axis.

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy ({data}`~ad_hoc_diffractometer.factories.BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.kappa4cv()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.presets.kappa4cv` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L871) for the complete stage
and mode configuration.

## Stage layout

### Stage coupling

```{graphviz}
digraph kappa4cv {
    rankdir=BT;
    label="kappa4cv";
    labelloc=t;
    fontsize=14;
    node [shape=box, style=filled, fontsize=11];

    komega [label="komega\naxis: -lateral\nLH", fillcolor="#a8d8ea"];
    kappa [label="kappa\naxis: [0.766, 0, 0.6428]\nRH", fillcolor="#a8d8ea"];
    kphi [label="kphi\naxis: -lateral\nLH", fillcolor="#a8d8ea"];
    ttheta [label="ttheta\naxis: -lateral\nLH", fillcolor="#f8a5a5"];

    { rank=same; komega; ttheta; }

    kappa -> komega;
    kphi -> kappa;

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

![kappa4cv stage axes](../_static/geometries/kappa4cv/kappa4cv_all.svg)

### Per-stage axis diagrams

::::{tab-set}
:::{tab-item} komega
![komega axis](../_static/geometries/kappa4cv/komega.svg)
:::
:::{tab-item} kappa
![kappa axis](../_static/geometries/kappa4cv/kappa.svg)
:::
:::{tab-item} kphi
![kphi axis](../_static/geometries/kappa4cv/kphi.svg)
:::
:::{tab-item} ttheta
![ttheta axis](../_static/geometries/kappa4cv/ttheta.svg)
:::
::::

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

**Virtual Eulerian angles** (computed from real kappa angles via Walko 2016 eq. [16]):
omega, chi, phi.  Used as constraint names for stub modes; converted back to
komega, kappa, kphi by the kappa inversion solver (Issue I / #153).

## Diffraction modes

Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 1 constraint
(N − 3 = 1 for N = 4 DOF).
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

### `bisecting` *(default)*

{class}`~ad_hoc_diffractometer.mode.BisectConstraint`:
`komega = ttheta / 2` (approximates the bisecting condition).

> **Note:** The correct bisecting condition is virtual `omega_euler = 0`,
> which differs from `komega = ttheta/2`. Corrected in Issue I / #153.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | — |

### `fixed_kphi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
`kphi` held at declared value (default 0°) — real stage, no kappa inversion needed.
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, ttheta |
| **Constant during** `forward()` | kphi |

### `constant_omega`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
Fix the virtual Eulerian omega at declared value (default 0°).
Uses the kappa Newton-Raphson solver — the caller chooses the value by
constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | omega (virtual) |

### `constant_chi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
Fix the virtual Eulerian chi at declared value (default 90°).
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | chi (virtual) |

### `constant_phi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
Fix the virtual Eulerian phi at declared value (default 0°).
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | phi (virtual) |

### `psi_constant`

{class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`:
azimuthal angle ψ validation filter.
Set ``g.azimuthal_reference = (h, k, l)`` before calling ``forward()``.
Returns bisecting solutions only when the natural ψ for (h,k,l) matches
the stored target.  See {doc}`../howto/surface`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

### `double_diffraction`

Full 4D simultaneous solver: finds motor angles where both the primary
(h₁,k₁,l₁) and secondary (h₂,k₂,l₂) reflections satisfy the Ewald
sphere condition.  Set ``mode.extras['h2']``, ``['k2']``, ``['l2']``
before calling ``forward()``.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Extras (input)** | h₂, k₂, l₂ (secondary reflection Miller indices) |

## API reference

- {func}`~ad_hoc_diffractometer.presets.kappa4cv`
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
