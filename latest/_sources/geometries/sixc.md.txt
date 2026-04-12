(geometry-sixc)=
# sixc — Six-Circle Surface (Lohmeier & Vlieg 1993)

Six-circle surface diffractometer. Sample and detector share a common alpha (rotary table) base stage. Designed for surface diffraction. Also known as the IUCr six-circle geometry.

**Walko (2016) designation:** (S3D2)1

**Coordinate basis:** You (1999) (``BASIS_YOU``): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.sixc()
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

- {func}`~ad_hoc_diffractometer.sixc`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- Lohmeier & Vlieg, *J. Appl. Cryst.* **26**, 706–716 (1993). DOI: `10.1107/S0021889893006198 <https://doi.org/10.1107/S0021889893006198>`_
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
