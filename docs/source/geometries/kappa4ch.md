(geometry-kappa4ch)=
# kappa4ch — Kappa Four-Circle (Laboratory)

Four-circle kappa diffractometer, horizontal scattering plane. Kappa axis tilted at α = 50° from vertical. Laboratory convention.

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy ({data}`~ad_hoc_diffractometer.factories.BASIS_BL`): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.kappa4ch()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Pre-built geometry definition

This geometry is defined by the {func}`~ad_hoc_diffractometer.presets.kappa4ch` factory
function — see the [source](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/src/ad_hoc_diffractometer/factories.py#L964) for the complete stage
and mode configuration.

## Stage layout

### Stage coupling

```{graphviz}
digraph kappa4ch {
    rankdir=BT;
    label="kappa4ch";
    labelloc=t;
    fontsize=14;
    node [shape=box, style=filled, fontsize=11];

    komega [label="komega\naxis: -vertical\nLH", fillcolor="#a8d8ea"];
    kappa [label="kappa\naxis: [0.766, 0, 0.6428]\nRH", fillcolor="#a8d8ea"];
    kphi [label="kphi\naxis: -vertical\nLH", fillcolor="#a8d8ea"];
    ttheta [label="ttheta\naxis: -vertical\nLH", fillcolor="#f8a5a5"];

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

![kappa4ch stage axes](../_static/geometries/kappa4ch/kappa4ch_all.svg)

### Per-stage axis diagrams

::::{tab-set}
:::{tab-item} komega
![komega axis](../_static/geometries/kappa4ch/komega.svg)
:::
:::{tab-item} kappa
![kappa axis](../_static/geometries/kappa4ch/kappa.svg)
:::
:::{tab-item} kphi
![kphi axis](../_static/geometries/kappa4ch/kphi.svg)
:::
:::{tab-item} ttheta
![ttheta axis](../_static/geometries/kappa4ch/ttheta.svg)
:::
::::

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

### `constant_omega`

Fix virtual Eulerian omega at declared value (default 0°) — see {doc}`kappa4cv` for details.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | omega (virtual) |

### `constant_chi`

Fix virtual Eulerian chi at declared value (default 90°).

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | chi (virtual) |

### `constant_phi`

Fix virtual Eulerian phi at declared value (default 0°).

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
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

## API reference

- {func}`~ad_hoc_diffractometer.presets.kappa4ch`
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
