(howto-orient)=
# Orient a Crystal

This guide shows how to compute the U and UB matrices from orienting
reflections.  See [Concepts](../concepts.md) for background on what
B, U, and UB mean.

## Overview

Crystal orientation requires:

1. A geometry with a wavelength set
2. A sample with a lattice defined
3. One, two, or three measured reflections with known (hkl) and motor angles

## Set up the geometry and sample

```python
import ad_hoc_diffractometer as ahd

g = ahd.make_geometry("fourcv")
g.wavelength = 1.5406   # Å  (Cu Kα)
g.sample.lattice = ahd.Lattice(a=5.431)  # cubic silicon
```

## Add orienting reflections

```python
# Add the primary reflection at its measured motor angles
g.sample.reflections.add("or1", hkl=(0, 0, 4),
                          angles={"ttheta": 69.13, "omega": 34.56,
                                  "chi": 90.0, "phi": 0.0})

# Add a second reflection 90° away in phi
g.sample.reflections.add("or2", hkl=(2, 2, 0),
                          angles={"ttheta": 47.30, "omega": 23.65,
                                  "chi": 35.26, "phi": 90.0})

# Designate which reflections to use for the UB calculation
g.sample.reflections.setor0("or1")
g.sample.reflections.setor1("or2")
```

## Compute UB from two reflections (Busing & Levy 1967)

```python
UB = ahd.ub_from_two_reflections_bl1967(g.sample)
print("UB =\n", UB)
```

This is the standard method for most beamline work.

## Compute UB from one reflection

When only one reflection is available, a provisional UB can be computed
by specifying a reference crystal direction and a reference stage:

```python
UB = ahd.ub_from_one_reflection(
    g.sample,
    reflection=g.reflections["or1"],
    reference_hkl=(0, 0, 1),   # crystal direction to align
    reference_stage="chi",      # stage axis to align it along
)
```

## Compute UB from three reflections (direct method)

For unknown lattices, UB can be determined directly from three reflections
without prior knowledge of the unit cell:

```python
g.sample.reflections.add("or3", hkl=(1, 1, 3),
                          angles={"ttheta": 56.12, "omega": 28.06,
                                  "chi": 58.0, "phi": 45.0})

UB = ahd.ub_from_three_reflections_bl1967(
    g.sample,
    g.sample.reflections["or1"],
    g.sample.reflections["or2"],
    g.sample.reflections["or3"],
)
```

## Use an identity U matrix

For a perfectly aligned crystal (crystal axes parallel to diffractometer
axes), use the identity orientation:

```python
ahd.ub_identity(g.sample)
```

## Verify the orientation — direction check

Check that `UB @ h` points in the same direction as the observed Q:

```python
import numpy as np

def direction_check(geometry, name):
    r = geometry.sample.reflections[name]
    q_phi = ahd.orientation.angles_to_phi_vector(geometry, **r.angles)
    ub_h  = geometry.sample.UB @ r.hkl
    cos_theta = np.dot(q_phi, ub_h) / (np.linalg.norm(q_phi) * np.linalg.norm(ub_h))
    angle_deg = np.degrees(np.arccos(np.clip(cos_theta, -1, 1)))
    print(f"{name}: {angle_deg:.3f}° discrepancy")

direction_check(g, "or1")
direction_check(g, "or2")
```

`0.0°` means perfect agreement.

## Full worked example

See the {doc}`fourcv_alignment_howto` notebook for a complete
step-by-step crystal alignment on a real instrument.

## See also

- {func}`~ad_hoc_diffractometer.orientation.ub_from_two_reflections_bl1967`
- {func}`~ad_hoc_diffractometer.orientation.ub_from_one_reflection`
- {func}`~ad_hoc_diffractometer.orientation.ub_from_three_reflections_bl1967`
- {func}`~ad_hoc_diffractometer.orientation.ub_identity`
- {func}`~ad_hoc_diffractometer.orientation.angles_to_phi_vector`
- {class}`~ad_hoc_diffractometer.reflection.ReflectionList`
