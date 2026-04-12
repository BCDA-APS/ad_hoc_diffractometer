(geometry-fivec)=
# fivec — Five-Circle (Vlieg et al. 1987)

Five-circle diffractometer: a standard fourcv (Eulerian four-circle) mounted on a vertical mu base stage. Sample and detector are coupled through mu, providing access to wider regions of reciprocal space.

**Walko (2016) designation:** (S3D1)1

**Coordinate basis:** You (1999) (``BASIS_YOU``): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.fivec()
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
   * - ``mu``
     - +vertical (+x)
     - right-handed, shared base
   * - ``omega``
     - −lateral (−z)
     - left-handed
   * - ``chi``
     - +longitudinal (+y)
     - right-handed
   * - ``phi``
     - −lateral (−z)
     - left-handed

**Detector stages** (floor first):

.. list-table::
   :header-rows: 1
   :widths: 15 40 20

   * - Stage
     - Axis
     - Handedness
   * - ``ttheta``
     - −lateral (−z)
     - left-handed

**Shared stage:** mu (base stage shared between sample and detector stacks)

## Diffraction modes

Available modes: *(none defined)*

## API reference

- {func}`~ad_hoc_diffractometer.fivec`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- Vlieg et al., *J. Appl. Cryst.* **20**, 330–337 (1987). DOI: `10.1107/S0021889887087266 <https://doi.org/10.1107/S0021889887087266>`_
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
