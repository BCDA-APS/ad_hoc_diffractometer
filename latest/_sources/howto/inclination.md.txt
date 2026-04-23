(howto-inclination)=
# Set Diffractometer Inclination

This guide explains how to account for a diffractometer that is mounted
at a non-zero angle relative to the incident beam.

## When is inclination needed?

In the standard geometry the incident beam travels along the
**longitudinal** axis of the diffractometer coordinate frame.
Some experimental setups violate this assumption:

- The beam is deflected upwards by a grazing-incidence mirror but the
  diffractometer cannot be re-levelled.
- The instrument is intentionally tilted for a non-standard sample
  environment.
- The beam path includes a pre-monochromator deflection that is not
  corrected mechanically.
- The beam path is deflected by an optical device, such as transfocator, before
  the sample.

In each case the effective beam direction in the diffractometer frame
is no longer the longitudinal axis.  The **inclination matrix** encodes
this offset as a 3×3 rotation matrix applied to the beam direction
before any diffraction calculation.

## Set the inclination

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.fourcv()

# Tilt 2° about the transverse axis (+x in the BL1967 basis)
TRANSVERSE = g.basis["transverse"]
g.set_inclination(axis=TRANSVERSE, angle_deg=2.0)
```

{meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.set_inclination`
takes a rotation axis (any non-zero vector; it is normalised internally)
and an angle in degrees.  Internally it builds the rotation matrix using
the Rodrigues formula
({func}`~ad_hoc_diffractometer.rotation.rotation_matrix`).

## Inspect the inclination matrix

```python
import numpy as np

print(np.round(g.inclination_matrix, 4))
# [[ 1.      0.      0.    ]
#  [ 0.      0.9994 -0.0349]
#  [ 0.      0.0349  0.9994]]
```

The
{attr}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.inclination_matrix`
property returns the current 3×3 matrix.  The default is the identity
(no inclination).

## Reset to zero inclination

Pass any non-zero axis with `angle_deg=0`:

```python
g.set_inclination(axis=TRANSVERSE, angle_deg=0.0)

np.array_equal(g.inclination_matrix, np.eye(3))  # True
```

Or assign the identity matrix directly:

```python
g.inclination_matrix = np.eye(3)
```

## How inclination affects calculations

The inclination matrix **R** modifies the effective incident-beam
direction.  Instead of using the longitudinal basis vector **ŷ**
directly, the calculation uses **R**\ :sup:`T` **ŷ** as the beam
direction:

$$
\hat{y}_{\text{eff}} = R^{T} \, \hat{y}
$$

This enters the scattering-vector computation in
{func}`~ad_hoc_diffractometer.orientation.angles_to_phi_vector`:

$$
\mathbf{Q}_{\text{lab}} = \frac{2\pi}{\lambda}
\bigl( D \, \hat{y}_{\text{eff}} - \hat{y}_{\text{eff}} \bigr)
$$

where *D* is the composite detector rotation matrix.  With zero
inclination (*R* = *I*), this reduces to the standard formula.

All calculations that depend on
{func}`~ad_hoc_diffractometer.orientation.angles_to_phi_vector` —
including {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.forward`,
{meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.inverse`,
and {func}`~ad_hoc_diffractometer.scan.psi_trajectory` — automatically
respect the inclination setting.

## Example: tilted four-circle

```python
import ad_hoc_diffractometer as ahd

g = ahd.presets.fourcv()
g.wavelength = 1.5406
g.sample.lattice = ahd.Lattice(a=5.431)
ahd.ub_identity(g.sample)

# Standard (no inclination): solve for (0 0 4)
sol_std = g.forward(0, 0, 4)
print("Standard:", sol_std[0])

# Tilt 2° about the transverse axis
g.set_inclination(axis=g.basis["transverse"], angle_deg=2.0)
sol_tilt = g.forward(0, 0, 4)
print("Tilted:  ", sol_tilt[0])
```

The tilted solution differs from the standard one because the solver
must compensate for the beam no longer being perfectly aligned with the
longitudinal axis.

## Validation

The inclination matrix must be a **proper rotation**:

- Shape (3, 3)
- Orthonormal: *R*\ :sup:`T` *R* = *I*  (within 10\ :sup:`-8`)
- det(*R*) = +1  (no reflections)

Assigning an invalid matrix raises {class}`ValueError`:

```python
import numpy as np

g.inclination_matrix = np.eye(4)          # wrong shape → ValueError
g.inclination_matrix = 2 * np.eye(3)      # not orthonormal → ValueError
g.inclination_matrix = -np.eye(3)         # det = -1 → ValueError
```

## Serialisation

The inclination matrix is included in
{meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.to_dict` as a
nested list and restored by
{meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.from_dict`:

```python
import json

g.set_inclination(axis=[0, 1, 0], angle_deg=3.0)
state = g.to_dict()
print(json.dumps(state["inclination_matrix"], indent=2))

g2 = ahd.AdHocDiffractometer.from_dict(state)
np.testing.assert_allclose(g2.inclination_matrix, g.inclination_matrix)
```

See {doc}`serialize` for the full save/restore workflow.

## See also

- {attr}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.inclination_matrix`
- {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.set_inclination`
- {func}`~ad_hoc_diffractometer.rotation.rotation_matrix`
- {func}`~ad_hoc_diffractometer.orientation.angles_to_phi_vector`
- {doc}`serialize`

## References

- D.A. Walko, *Ref. Module Mater. Sci. Mater. Eng.* (2016), §4.2 —
  general-inclination geometries
