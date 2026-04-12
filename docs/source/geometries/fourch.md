(geometry-fourch)=
# fourch — Four-Circle Eulerian (Laboratory)

Busing & Levy (1967) four-circle Eulerian diffractometer, horizontal scattering plane. ω and 2θ rotate about the vertical axis. Standard laboratory convention.

**Walko (2016) designation:** S3D1

**Coordinate basis:** Busing & Levy (``BASIS_BL``): lateral=+x, longitudinal=+y, vertical=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourch()
g.wavelength = 1.0  # Å
print(g.summary())
```

## Stage layout

**Sample stages** (floor first):

.. list-table::
   :header-rows: 1
   :widths: 15 40 20

   * - Stage
     - Axis
     - Handedness
   * - ``omega``
     - −vertical (−z BL)
     - left-handed
   * - ``chi``
     - +longitudinal (+y BL)
     - right-handed
   * - ``phi``
     - −vertical (−z BL)
     - left-handed

**Detector stages** (floor first):

.. list-table::
   :header-rows: 1
   :widths: 15 40 20

   * - Stage
     - Axis
     - Handedness
   * - ``ttheta``
     - −vertical (−z BL)
     - left-handed

## Diffraction modes

Available modes: ``bisecting``, ``fixed_chi``, ``fixed_phi``

## API reference

- {func}`~ad_hoc_diffractometer.fourch`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- Busing & Levy, *Acta Cryst.* **22**, 457–464 (1967). DOI: `10.1107/S0365110X67000970 <https://doi.org/10.1107/S0365110X67000970>`_
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
