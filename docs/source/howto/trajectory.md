(howto-trajectory)=
# Plan a Trajectory

This guide shows how to compute motor-angle sequences along paths through
reciprocal space and for ψ scans.

## Prerequisites

- A geometry with wavelength and UB matrix set (see {doc}`orient`)

## hkl trajectory — path between two reflections

```python
import ad_hoc_diffractometer as ahd

g = ahd.fourcv()
g.wavelength = 1.5406
g.sample.lattice = ahd.Lattice(a=5.431)
ahd.ub_identity(g.sample)
g.mode_name = "bisecting"

# Define a path in reciprocal space
trajectory = [
    {"hkl_start": (1, 0, 0), "hkl_end": (2, 0, 0)},
]

points = ahd.hkl_trajectory(g, trajectory, n_points=11)
for p in points:
    print(p)
```

Each point is a dict with keys `h`, `k`, `l`, and one key per motor stage.

## ψ scan — rotate about the scattering vector

A ψ scan rotates the sample about Q while keeping the Bragg condition
satisfied.  This is used to probe azimuthal anisotropy.

```python
import numpy as np

psi_values = np.linspace(-180, 180, 37)  # degrees
points = ahd.psi_trajectory(g, 1, 1, 0, psi_values)
for p in points:
    print(p["psi"], p["omega"], p["chi"], p["phi"])
```

## Trajectory plan — with limit checking

`trajectory_plan()` flags points that violate motor limits or are
inaccessible (Q exceeds the Ewald sphere):

```python
plan = ahd.trajectory_plan(
    g,
    hkl_start=(1, 0, 0),
    hkl_end=(3, 0, 0),
    n_points=21,
)
for p in plan:
    if p.get("reachable", True):
        print("OK ", p["omega"])
    else:
        print("LIMIT VIOLATION at", p)
```

## Solution selection — avoiding branch flips

All three functions accept a `solution_key` parameter that scores
candidates against the previous point, preventing discontinuous jumps
between chi branches:

```python
points = ahd.hkl_trajectory(
    g, trajectory, n_points=51,
    solution_key=ahd.NEAREST_ANGLES,   # default; minimises motor travel
)
```

## See also

- {func}`~ad_hoc_diffractometer.hkl_trajectory`
- {func}`~ad_hoc_diffractometer.psi_trajectory`
- {func}`~ad_hoc_diffractometer.trajectory_plan`
- {data}`~ad_hoc_diffractometer.NEAREST_ANGLES`
- {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.psi`
