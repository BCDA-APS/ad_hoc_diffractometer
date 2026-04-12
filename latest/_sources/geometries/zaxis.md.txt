(geometry-zaxis)=
# zaxis — Z-Axis Four-Circle (Surface)

Z-axis four-circle diffractometer for surface diffraction. The sample surface normal is parallel to the Z-axis. Sample and detector share an alpha base stage. Designed for grazing-incidence work.

**Walko (2016) designation:** (S1D2)1

**Coordinate basis:** You (1999) (``BASIS_YOU``): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.zaxis()
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
   * - ``alpha``
     - +vertical (+x)
     - right-handed, shared base
   * - ``Z``
     - +longitudinal (+y)
     - right-handed

**Detector stages** (floor first):

.. list-table::
   :header-rows: 1
   :widths: 15 40 20

   * - Stage
     - Axis
     - Handedness
   * - ``delta``
     - −lateral (−z)
     - left-handed
   * - ``gamma``
     - +vertical (+x)
     - right-handed

**Shared stage:** alpha (base stage shared between sample and detector stacks)

## Diffraction modes

Available modes: *(none defined)*

## API reference

- {func}`~ad_hoc_diffractometer.zaxis`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- Bloch, *J. Appl. Cryst.* **18**, 33–36 (1985). DOI: `10.1107/S0021889885009858 <https://doi.org/10.1107/S0021889885009858>`_
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
