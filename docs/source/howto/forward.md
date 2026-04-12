(howto-forward)=
# Solve the Forward Problem

The **forward problem** finds the motor angles that satisfy the Bragg
condition for a given reflection (hkl → motor angles).

## Prerequisites

- A geometry with a wavelength set
- A sample with a UB matrix set (see {doc}`orient`)

## Basic usage

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourcv()
g.wavelength = 1.5406   # Å
g.sample.lattice = ahd.Lattice(a=5.431)
ahd.ub_identity(g.sample)

solutions = g.forward(1, 1, 0)
for s in solutions:
    print(s)
```

`forward()` returns a list of dicts, each mapping stage name → angle (degrees).
Multiple solutions exist when different motor configurations satisfy the
Bragg condition (e.g. positive and negative chi branches).

## Select a solution

```python
# Take the first solution
angles = g.forward(1, 1, 0)[0]
print(angles)
# {'omega': 23.65, 'chi': 35.26, 'phi': 0.0, 'ttheta': 47.30}
```

## Predict the Bragg angle only

To get d-spacing and 2θ without motor angles:

```python
d   = ahd.hkl_to_d(g, 1, 1, 0)
tth = ahd.hkl_to_two_theta(g, 1, 1, 0)
print(f"d = {d:.4f} Å,  2θ = {tth:.3f}°")
```

## Apply motor limits

Stage limits are enforced automatically by `forward()`.  To set limits:

```python
g.stages["chi"].limits = (-10.0, 100.0)   # degrees
```

Solutions outside the limits are filtered out.

## Use with diffraction modes

Set an active mode to restrict which angles are free:

```python
g.mode_name = "bisecting"   # omega = ttheta/2
solutions = g.forward(0, 0, 2)
```

See {doc}`modes` for details.

## See also

- {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.forward`
- {func}`~ad_hoc_diffractometer.hkl_to_d`
- {func}`~ad_hoc_diffractometer.hkl_to_two_theta`
- {doc}`modes`
- {doc}`orient`
