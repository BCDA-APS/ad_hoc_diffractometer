(geometry-psic)=
# psic — Six-Circle 4S+2D (You 1999)

You (1999) 4S+2D six-circle diffractometer. Four sample stages (mu, eta, chi, phi) and two detector stages (nu, delta). Lateral detector, vertical scattering plane. Standard synchrotron six-circle.

**Coordinate basis:** You (1999) (``BASIS_YOU``): vertical=+x, longitudinal=+y, lateral=+z.

## Quick start

```python
import ad_hoc_diffractometer as ahd

g = ahd.psic()
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
     - right-handed
   * - ``eta``
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
   * - ``nu``
     - +vertical (+x)
     - right-handed
   * - ``delta``
     - −lateral (−z)
     - left-handed

## Diffraction modes

Available modes: ``bisecting``, ``fixed_chi``, ``fixed_phi``, ``fixed_mu``

## API reference

- {func}`~ad_hoc_diffractometer.psic`
- {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer`

## References

- You, *J. Appl. Cryst.* **32**, 614–623 (1999). DOI: `10.1107/S0021889899001223 <https://doi.org/10.1107/S0021889899001223>`_
- Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016).
