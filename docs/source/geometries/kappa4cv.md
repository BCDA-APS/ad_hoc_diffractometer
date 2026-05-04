(geometry-kappa4cv)=
# kappa4cv — Kappa Four-Circle (vertical)

Four-circle kappa diffractometer, vertical scattering plane. The chi
circle is replaced by a kappa arm tilted at α = 50° from the
**transverse** axis toward the **vertical** axis, lying entirely in
the transverse–vertical plane (see
[How the kappa axis is defined](#kappa4cv-axis-definition) below).

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy ({data}`~ad_hoc_diffractometer.factories.BASIS_BL`): transverse=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.kappa4cv()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Demo geometry definition

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
| ``komega`` | transverse (−x in BL) | left-handed | base |
| ``kappa`` | (cos α)·transverse + (sin α)·vertical = (cos α)·x̂ + (sin α)·ẑ in BL (α = 50°) | right-handed | ``komega`` |
| ``kphi`` | transverse (−x in BL) | left-handed | ``kappa`` |

**Detector stages (base first):**

| Stage | Axis | Handedness | Parent |
|---|---|---|---|
| ``ttheta`` | transverse (−x in BL) | left-handed | base |

(kappa4cv-axis-definition)=
### How the kappa axis is defined

The kappa rotation axis is **inclined by α from the omega axis**,
lying in the plane that contains both omega and the
equivalent-Eulerian chi axis, and tilted from omega toward that chi
direction (Walko 2016 §4.1; ITC Vol. C §2.2.6.2: *"the κ axis is
inclined at 50° to the ω axis"*).

For ``kappa4cv`` the omega axis lies along the **transverse** line
and the equivalent-Eulerian chi axis lies along the **vertical**
direction.  The kappa arm therefore lies in the **transverse–
vertical plane**, tilted by ``α`` from the transverse line toward
+V:

$$
\hat{n}_{\kappa} \;=\; \cos\alpha \cdot \hat{T} \;+\; \sin\alpha \cdot \hat{V}.
$$

In the BL basis (T=+x̂, V=+ẑ) at α = 50° this is

$$
\hat{n}_{\kappa} \;=\; \cos 50° \cdot \hat{x} \;+\; \sin 50° \cdot \hat{z} \;=\; (0.643,\, 0,\, 0.766).
$$

The longitudinal component is *exactly zero*: the kappa arm rises
upward and outward in the vertical plane perpendicular to the
incident beam.

**Published reference figures.**  The kappa-arm orientation
described above matches Walko (2016) Figure 3 directly and the
canonical κ-goniostat axes in Thorkildsen *et al.* (2006) Table 1
and equation (3).  Wyckoff (1985) Figure 2(b) on p. 334 shows the
analogous 90° rotation onto the horizontal scattering plane (the
``kappa4ch`` convention).

The omega axis sense (``komega = −transverse`` = left-handed
rotation about +T) follows Walko's left-handed sign convention.
ITC Vol. C §2.2.6.2 prefers a right-handed sign convention for
``omega/chi/phi``; the two are equivalent up to motor-angle sign
flips.  See the {func}`~ad_hoc_diffractometer.presets` module
docstring for further discussion.

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
decomposition.

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
Solved analytically via the equivalent-Eulerian dispatch — the
caller chooses the value by constructing a
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

- D.A. Walko, *Multicircle Diffractometry Methods*, in *Reference
  Module in Materials Science and Materials Engineering* (Elsevier,
  2016), §4.1, Figure 3, equation [16].
  DOI: [10.1016/B978-0-12-803581-8.01215-7](https://doi.org/10.1016/B978-0-12-803581-8.01215-7).
- G. Thorkildsen, H.B. Larsen & J.A. Beukes, *Angle calculations
  for a three-circle goniostat*, J. Appl. Cryst. **39**, 151–157
  (2006), Table 1 and equation (3) (canonical κ-goniostat axes).
  DOI: [10.1107/S0021889805041877](https://doi.org/10.1107/S0021889805041877).
- W.R. Busing & H.A. Levy, *Angle calculations for 3- and 4-circle
  X-ray and neutron diffractometers*, Acta Cryst. **22**, 457–464
  (1967) (BL coordinate basis).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970).
- *International Tables for Crystallography*, Vol. C, §2.2.6
  (2006), p. 36 (α = 50° convention; cites Wyckoff 1985 for the
  schematic picture).
  DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577).
