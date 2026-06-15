(howto-forward)=
# Forward and Inverse Computations

The diffractometer has two complementary computations:

- **Forward** (`g.forward(h, k, l)`) — given Miller indices, find the motor
  angles that satisfy the Bragg condition (hkl → motor angles).
- **Inverse** (`g.inverse(**angles)`) — given a set of motor angles, find
  the Miller indices of the reflection currently in the Bragg condition
  (motor angles → hkl).

This guide covers the forward computation.  The inverse is used, for example,
after manually positioning the diffractometer to identify an unknown peak.

## Prerequisites

- A geometry with a wavelength set
- A sample with a UB matrix set (see {doc}`orient`)

## Set a diffraction mode

Before calling `forward()`, set the active diffraction mode.  The mode
controls which stages are free, fixed, or coupled, and determines which
solutions are physically meaningful:

```python
import ad_hoc_diffractometer as ahd

g = ahd.make_geometry("fourcv")
g.wavelength = 1.5406   # Å
g.sample.lattice = ahd.Lattice(a=5.431)
ahd.ub_identity(g.sample)

g.mode_name = "bisecting"   # omega = ttheta/2; standard synchrotron mode
```

For fixed-angle modes, preset the stage angle **before** activating the mode:

```python
g.set_angle("chi", 90.0)   # preset chi to the desired fixed value
g.mode_name = "fixed_chi"  # chi will be held at 90° during forward()
```

Without a mode, all stages are free and the solver returns all geometrically
valid solutions.  See {doc}`modes` for the full list.

## Forward computation

```python
solutions = g.forward(1, 1, 0)
for s in solutions:
    print(s)
```

`forward()` returns a **list** of dicts, each mapping stage name → angle
(degrees).  Multiple solutions exist because the Bragg condition fixes only
the direction of the scattering vector **Q** in the laboratory frame — it
does not uniquely determine the motor angles.  For a four-circle geometry,
the same **Q** can typically be reached with chi > 0 or chi < 0 (two
branches), and for each branch there are infinitely many (omega, phi) pairs
that satisfy the condition.  The solver samples a finite set of
representative solutions; the active diffraction mode filters these to those
that satisfy any additional constraints.  It is the caller's responsibility
to select the physically appropriate solution for their experimental setup.

## Select a solution

```python
# Take the first solution
angles = g.forward(1, 1, 0)[0]
print(angles)
# {'omega': 23.65, 'chi': 35.26, 'phi': 0.0, 'ttheta': 47.30}
```

## Predict the Bragg angle only

To get d-spacing and 2θ without computing motor angles:

```python
d   = ahd.hkl_to_d(g, 1, 1, 0)
tth = ahd.hkl_to_two_theta(g, 1, 1, 0)
print(f"d = {d:.4f} Å,  2θ = {tth:.3f}°")
```

## Apply motor limits

Stage limits are enforced automatically by `forward()`.  To set limits:

```python
g.stage("chi").limits = (-10.0, 100.0)   # degrees
```

Solutions outside the limits are filtered out.

## Inverse computation

`inverse()` requires a UB matrix to be set on the sample (see {doc}`orient`).
Given a set of motor angles, it returns the Miller indices of the reflection
currently in the Bragg condition.

```{note}
`inverse()` always returns a **unique** hkl: once UB is established, a single
matrix multiplication maps motor angles unambiguously to reciprocal space.
`forward()` is the reverse: the Bragg condition constrains only the direction
of **Q** in the laboratory frame, not the individual motor angles, so the
result is a **list** that may contain anywhere from 0 (reflection
unreachable) to a geometry-dependent maximum — typically 2–12 solutions for
four- and six-circle geometries depending on the number of free stages and
the active diffraction mode.
```

```python
# After driving the diffractometer to some position manually:
hkl = g.inverse(omega=23.65, chi=35.26, phi=0.0, ttheta=47.30)
print(hkl)   # (h, k, l) as a numpy array
```

This is useful for identifying an unknown peak found during a scan, or
for verifying that the diffractometer is positioned at the expected
reflection after moving motors.

## Handle errors from forward()

`forward()` raises specific exceptions for distinct failure modes:

### EwaldSphereViolation

Raised when the requested reflection cannot be reached at the current
wavelength — `|Q| > 4π/λ` regardless of motor angles:

```python
from ad_hoc_diffractometer import EwaldSphereViolation

try:
    solutions = g.forward(10, 10, 10)
except EwaldSphereViolation as e:
    print(f"|Q| = {e.q_mag:.4f} Å⁻¹")
    print(f"Ewald sphere limit = {e.q_max:.4f} Å⁻¹  (λ = {e.wavelength} Å)")
    print("Reduce wavelength or choose a smaller reflection.")
```

Attributes: `q_mag` (requested |Q| in Å⁻¹), `q_max` (4π/λ), `wavelength`.

### ConstraintViolation

Raised when the solver returns a solution that violates a declared
constraint beyond the display-precision tolerance.  This signals either a
solver bug or an unimplemented virtual-angle constraint:

```python
from ad_hoc_diffractometer import ConstraintViolation

try:
    solutions = g.forward(1, 0, 0)
except ConstraintViolation as e:
    print(f"Solution {e.solution_index} violated {e.constraint_repr}")
    print(f"Residual = {e.residual:.2e}°  (tolerance = {e.tolerance:.2e}°)")
```

Attributes: `solution_index`, `constraint_repr`, `residual`, `tolerance`.

### NotImplementedError

Raised when the active mode is `None` or its `is_implemented(geometry)`
returns `False` (e.g. a prerequisite like ``azimuth`` or
``surface_normal`` is not set):

```python
g.mode_name = None
try:
    g.forward(1, 0, 0)
except NotImplementedError as e:
    print(e)

# fixed_psi requires azimuth to be set
g.mode_name = "fixed_psi"
g.azimuth = None
try:
    g.forward(1, 0, 0)
except NotImplementedError as e:
    print(e)  # describes which prerequisite is missing
```

## See also

- {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward`
- {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.inverse`
- {func}`~ad_hoc_diffractometer.conversions.hkl_to_d`
- {func}`~ad_hoc_diffractometer.conversions.hkl_to_two_theta`
- {exc}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {exc}`~ad_hoc_diffractometer.mode.ConstraintViolation`
- {doc}`modes`
- {doc}`constraints`
- {doc}`orient`
