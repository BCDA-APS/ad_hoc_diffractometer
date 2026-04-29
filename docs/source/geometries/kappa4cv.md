(geometry-kappa4cv)=
# kappa4cv — Kappa Four-Circle (Synchrotron)

Four-circle kappa diffractometer, vertical scattering plane. The chi
circle is replaced by a kappa axis tilted at α = 50° **from the
outer komega axis toward the equivalent Eulerian chi axis** (issue
#241; see [How the kappa axis is defined](#kappa4cv-axis-definition)
below).

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy ({data}`~ad_hoc_diffractometer.factories.BASIS_BL`): transverse=+x, longitudinal=+y, vertical=+z.

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

<iframe
  src="../_static/geometries/kappa4cv/kappa4cv.html"
  width="100%"
  height="1230px"
  style="border:none;"
  loading="lazy">
</iframe>

```{raw} html
<details>
<summary>Static fallback (click to expand if the interactive figure above is blank)</summary>
```

![kappa4cv stage layout](../_static/geometries/kappa4cv/kappa4cv.svg)

```{raw} html
</details>
```

**Sample stages (base first):**

| Stage | Axis | Handedness | Parent |
|---|---|---|---|
| ``komega`` | −transverse (−x BL) | left-handed | base |
| ``kappa`` | −x · cos α + ŷ · sin α (α = 50°) | right-handed | ``komega`` |
| ``kphi`` | −transverse (−x BL) | left-handed | ``kappa`` |

**Detector stages (base first):**

| Stage | Axis | Handedness | Parent |
|---|---|---|---|
| ``ttheta`` | −transverse (−x BL) | left-handed | base |

(kappa4cv-axis-definition)=
### How the kappa axis is defined

The kappa rotation axis lies in the plane spanned by the outer
``komega`` axis and the *equivalent Eulerian chi axis*, tilted by
``α`` from ``komega`` toward the chi-equivalent direction:

$$
\hat{n}_{\kappa} \;=\; \cos\alpha \cdot \hat{n}_{\kappa\omega} \;+\; \sin\alpha \cdot \hat{n}_{\chi,\,\text{eq}}.
$$

For ``kappa4cv`` (BL convention) this is

$$
\hat{n}_{\kappa} \;=\; \cos 50° \cdot (-\hat{x}) \;+\; \sin 50° \cdot (+\hat{y}) \;=\; (-0.643,\, 0.766,\, 0).
$$

This **differs** from earlier versions of the package, which set the
kappa axis from a textbook formula
``vertical · cos α + transverse · sin α`` regardless of the actual
``komega`` orientation.  The earlier formula is correct only when
``komega`` is along ``+vertical``; it is *not* correct for
``kappa4cv``, ``kappa4ch``, or ``kappa6c``, all of which encode
``komega`` along a non-vertical signed axis.  The mismatch caused a
silent solver gap that returned ``"No solutions"`` for several
physically reachable reflections; see issue #241 and
{class}`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention` for
the full derivation and rationale.

**Virtual Eulerian angles** ``omega``, ``chi``, ``phi`` are mapped
to / from the real motors via the geometry-aware decomposition in
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes` and
{func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian_axes`.  The
older Walko-textbook helpers
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa` and
{func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian` are retained
as reference implementations of the published closed form but are
**not** used inside the solver.

## Diffraction modes

Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 1 constraint
(N − 3 = 1 for N = 4 DOF).
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

### `bisecting` *(default)*

{class}`~ad_hoc_diffractometer.mode.VirtualBisectConstraint`:
``omega_virtual = ttheta / 2`` enforced on the **virtual** Eulerian
omega pseudoangle.  The kappa motor triple ``(komega, kappa, kphi)``
satisfies this constraint via the geometry-aware
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes`
decomposition (issue #241).

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

### `fixed_omega`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
Fix the virtual Eulerian omega at declared value (default 0°).
Solved analytically via the equivalent-Eulerian dispatch (issue
#241) — the caller chooses the value by constructing a
{class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | omega (virtual) |

### `fixed_chi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
Fix the virtual Eulerian chi at declared value (default 90°).
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | chi (virtual) |

### `fixed_phi`

{class}`~ad_hoc_diffractometer.mode.SampleConstraint`:
Fix the virtual Eulerian phi at declared value (default 0°).
The caller chooses the value by constructing a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | phi (virtual) |

### `fixed_psi`

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
- {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {class}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {class}`~ad_hoc_diffractometer.mode.ConstraintViolation`

## References

- ITC Vol. C §2.2.6 (2006). DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577)
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), eq. [16].
