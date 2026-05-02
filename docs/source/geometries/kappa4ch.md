(geometry-kappa4ch)=
# kappa4ch — Kappa Four-Circle (Laboratory)

Four-circle kappa diffractometer, horizontal scattering plane.  The
chi circle is replaced by a kappa arm tilted at α = 50° from the
**vertical** axis toward the **longitudinal** axis (along the
incident beam), lying entirely in the vertical–longitudinal plane.
Laboratory convention.

**Walko (2016) designation:** S3D1 (kappa)

**Coordinate basis:** Busing & Levy ({data}`~ad_hoc_diffractometer.factories.BASIS_BL`): transverse=+x, longitudinal=+y, vertical=+z.

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

<iframe
  src="../_static/geometries/kappa4ch/kappa4ch.html"
  width="100%"
  height="1230px"
  style="border:none;"
  loading="lazy">
</iframe>

```{raw} html
<details>
<summary>Static fallback (click to expand if the interactive figure above is blank)</summary>
```

![kappa4ch stage layout](../_static/geometries/kappa4ch/kappa4ch.svg)

```{raw} html
</details>
```

**Sample stages (base first):**

| Stage | Axis | Handedness | Parent |
|---|---|---|---|
| ``komega`` | vertical (−z in BL) | left-handed | base |
| ``kappa`` | (cos α)·vertical + (sin α)·longitudinal = (cos α)·ẑ + (sin α)·ŷ in BL (α = 50°) | right-handed | ``komega`` |
| ``kphi`` | vertical (−z in BL) | left-handed | ``kappa`` |

**Detector stages (base first):**

| Stage | Axis | Handedness | Parent |
|---|---|---|---|
| ``ttheta`` | vertical (−z in BL) | left-handed | base |

### How the kappa axis is defined

The kappa rotation axis is **inclined by α from the omega axis**,
lying in the plane that contains both omega and the
equivalent-Eulerian chi axis, and tilted from omega toward that chi
direction (Walko 2016 §4.1; ITC Vol. C §2.2.6.2: *"the κ axis is
inclined at 50° to the ω axis"*).

For ``kappa4ch`` the omega axis lies along the **vertical** line
and the equivalent-Eulerian chi axis lies along the **longitudinal**
direction (along the X-ray beam).  The kappa arm therefore lies in
the **vertical–longitudinal plane**, tilted by ``α`` from the
vertical line toward +L:

$$
\hat{n}_{\kappa} \;=\; \cos\alpha \cdot \hat{V} \;+\; \sin\alpha \cdot \hat{L}.
$$

In the BL basis (V=+ẑ, L=+ŷ) at α = 50° this is

$$
\hat{n}_{\kappa} \;=\; \cos 50° \cdot \hat{z} \;+\; \sin 50° \cdot \hat{y} \;=\; (0,\, 0.766,\, 0.643).
$$

The transverse component is *exactly zero*: the kappa arm rises
upward and toward the X-ray source/sample direction in the plane
that contains the vertical omega axis and the beam.

**Published reference figure.**  This convention matches Wyckoff
(1985) Figure 2(b) on p. 334 — the canonical Enraf-Nonius kappa
diffractometer in the horizontal-scattering laboratory layout.
ITC Vol. C §2.2.6 cites this figure as the schematic for the
kappa goniostat.  Walko (2016) Figure 3 shows the same instrument
rotated 90° onto the vertical scattering plane (the ``kappa4cv``
convention).

The omega axis sense (``komega = −vertical`` = left-handed
rotation about +V) follows Walko's left-handed sign convention.
ITC Vol. C §2.2.6.2 prefers a right-handed sign convention for
``omega/chi/phi``; the two are equivalent up to motor-angle sign
flips.

**Virtual Eulerian angles** ``omega``, ``chi``, ``phi`` are mapped
to / from the real motors via the geometry-aware decomposition in
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes` and
{func}`~ad_hoc_diffractometer.kappa.kappa_to_eulerian_axes`.

## Diffraction modes

Each mode is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet` of 1 constraint
(N − 3 = 1 for N = 4 DOF).
Identical mode set to {doc}`kappa4cv`.
See {doc}`../howto/modes` for usage and {doc}`../howto/constraints` for
changing constraint values at run time.

### `bisecting` *(default)*

{class}`~ad_hoc_diffractometer.mode.VirtualBisectConstraint`:
``omega_virtual = ttheta / 2`` enforced on the virtual Eulerian
omega pseudoangle.  Solved via the geometry-aware
{func}`~ad_hoc_diffractometer.kappa.eulerian_to_kappa_axes`
decomposition.

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

### `fixed_omega`

Fix virtual Eulerian omega at declared value (default 0°) — see {doc}`kappa4cv` for details.

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | omega (virtual) |

### `fixed_chi`

Fix virtual Eulerian chi at declared value (default 90°).

| | |
|---|---|
| **Computed** | komega, kappa, kphi, ttheta |
| **Constant during** `forward()` | chi (virtual) |

### `fixed_phi`

Fix virtual Eulerian phi at declared value (default 0°).

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
| **Extras (input)** | n̂ (reference vector), ψ (target azimuth, degrees) |
| **Extras (output)** | psi (computed azimuth) |

## API reference

- {func}`~ad_hoc_diffractometer.presets.kappa4ch`
- {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer`
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {class}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {class}`~ad_hoc_diffractometer.mode.ConstraintViolation`

## References

- H.W. Wyckoff, *Diffractometry*, in *Methods in Enzymology* **114**,
  330–386 (1985), Figure 2(b) on p. 334 (canonical Enraf-Nonius
  kappa diffractometer).
  DOI: [10.1016/0076-6879(85)14026-7](https://doi.org/10.1016/0076-6879(85)14026-7).
- D.A. Walko, *Multicircle Diffractometry Methods*, in *Reference
  Module in Materials Science and Materials Engineering* (Elsevier,
  2016), §4.1, equation [16].
  DOI: [10.1016/B978-0-12-803581-8.01215-7](https://doi.org/10.1016/B978-0-12-803581-8.01215-7).
- W.R. Busing & H.A. Levy, *Angle calculations for 3- and 4-circle
  X-ray and neutron diffractometers*, Acta Cryst. **22**, 457–464
  (1967) (BL coordinate basis).
  DOI: [10.1107/S0365110X67000970](https://doi.org/10.1107/S0365110X67000970).
- *International Tables for Crystallography*, Vol. C, §2.2.6
  (2006), p. 36 (α = 50° convention; cites Wyckoff 1985 for the
  schematic picture).
  DOI: [10.1107/97809553602060000577](https://doi.org/10.1107/97809553602060000577).
