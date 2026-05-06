(howto-surface)=
# Surface Geometry and the Reference Vector

This guide explains how to supply the surface normal n̂ and the azimuthal
reference vector to modes that require them, and how to compute the
resulting reference pseudo-angles (α_i, β_out, ψ, naz).

## Background

Several diffraction modes — `fixed_psi_vertical`, `zaxis`, `reflectivity`,
and others — require an external reference vector to complete their constraint.
In this package the reference vector is supplied as **Miller indices (h, k, l)**,
not as a lab-frame Cartesian vector.  The package converts to the lab frame
internally using the UB matrix.

Two separate reference vectors may be set:

- **`surface_normal`** — the direction perpendicular to the sample surface,
  used by `alpha_i`, `alpha_f`, `incidence_angle`, `exit_angle`, and
  surface modes (`zaxis`, `reflectivity`, `alpha_eq_beta_zaxis`).
- **`azimuthal_reference`** — the direction used to define ψ = 0, used by
  `psi_angle` and `fixed_psi_*` modes.

They may be the same vector (e.g. the surface normal is also the azimuthal
reference) or different.

## Set the surface normal

```python
import ad_hoc_diffractometer as ahd

g = ahd.make_geometry("psic")
g.wavelength = 1.0  # Å
g.sample.lattice = ahd.Lattice(a=4.0, c=6.5)
ahd.ub_identity(g.sample)

# Surface normal for a (001)-cut sample: Miller indices (0, 0, 1)
g.surface_normal = (0, 0, 1)
print(g.surface_normal)   # (0.0, 0.0, 1.0)

# Clear the surface normal
g.surface_normal = None
```

The setter accepts any three-element sequence of numbers and raises
`ValueError` for the zero vector or wrong shape.

## Set the azimuthal reference

```python
# Azimuthal reference along the a-axis: (1, 0, 0)
g.azimuthal_reference = (1, 0, 0)
print(g.azimuthal_reference)   # (1.0, 0.0, 0.0)

# Same vector as surface normal for a (001) surface
g.azimuthal_reference = (0, 0, 1)
```

## Compute incidence and exit angles

```python
from ad_hoc_diffractometer import incidence_angle, exit_angle

g.surface_normal = (0, 0, 1)
g.mode_name = "bisecting_vertical"
solutions = g.forward(1, 0, 0)

for sol in solutions:
    ai = incidence_angle(g, angles=sol)
    af = exit_angle(g, angles=sol)
    print(f"alpha_i = {ai:.4f}°   beta_out = {af:.4f}°")
```

Both functions use the current stage angles when `angles=None`:

```python
# Set current motor positions first
for name, value in solutions[0].items():
    g.set_angle(name, value)

ai = incidence_angle(g)   # uses current stage angles
```

## Compute the azimuthal angle ψ

```python
from ad_hoc_diffractometer import psi_angle

g.azimuthal_reference = (0, 0, 1)
g.mode_name = "bisecting_vertical"
solutions = g.forward(1, 0, 0)

for sol in solutions:
    psi = psi_angle(g, angles=sol)
    print(f"psi = {psi:.4f}°")
```

ψ is the angle between the azimuthal reference vector (projected onto the
plane perpendicular to Q) and the incident beam direction in that same plane.
ψ = 0 when the reference vector lies in the scattering plane on the
incident-beam side.

## Compute naz

```python
from ad_hoc_diffractometer import naz_angle

g.surface_normal = (0, 0, 1)
naz = naz_angle(g)   # uses current stage angles
print(f"naz = {naz:.4f}°")
```

naz is the azimuthal angle of the surface normal n̂ projected onto the
horizontal plane of the lab frame.

## Symmetric reflection condition (α_i = β_out)

```python
from ad_hoc_diffractometer import incidence_angle, exit_angle

g.surface_normal = (0, 0, 1)
g.mode_name = "bisecting_vertical"
solutions = g.forward(1, 0, 0)

for sol in solutions:
    ai = incidence_angle(g, angles=sol)
    af = exit_angle(g, angles=sol)
    is_sym = abs(ai - af) < 0.01  # within 0.01°
    print(f"alpha_i={ai:.3f}°  beta_out={af:.3f}°  symmetric={is_sym}")
```

Alternatively, use the built-in `is_specular()` method on the geometry:

```python
for sol in solutions:
    for name, value in sol.items():
        g.set_angle(name, value)
    print(f"specular: {g.is_specular()}")
```

## Serialization

`surface_normal` and `azimuthal_reference` are serialized in `to_dict()`
and restored by `from_dict()`:

```python
import json
from ad_hoc_diffractometer import AdHocDiffractometer

g.surface_normal = (0, 0, 1)
g.azimuthal_reference = (1, 0, 0)

d = g.to_dict()
print(d["surface_normal"])        # [0.0, 0.0, 1.0]
print(d["azimuthal_reference"])   # [1.0, 0.0, 0.0]

g2 = AdHocDiffractometer.from_dict(d)
print(g2.surface_normal)          # (0.0, 0.0, 1.0)
```

## Reference constraint modes

Modes that use a ``ReferenceConstraint`` require the appropriate reference
vector to be set on the geometry.  ``fixed_psi_*`` modes require
``azimuthal_reference``; surface modes (``zaxis``, ``reflectivity``) require
``surface_normal``.

```python
g.azimuthal_reference = (0, 0, 1)
g.mode_name = "fixed_psi_vertical"

cs = g.modes["fixed_psi_vertical"]
rc = cs.reference_constraint

print(rc.has_reference_vector(g))   # True  — vector is set
print(rc.is_implemented(g))          # True  — solver available

# forward() returns bisecting solutions whose natural ψ matches the target
solutions = g.forward(1, 0, 0)
```

The ``fixed_psi`` solver acts as a **validation filter**: ψ is a pure
phi-frame quantity that is the same for every Bragg solution of a given
(h,k,l) and UB.  The solver computes the natural ψ and returns solutions
only if it matches the stored target — otherwise it returns an empty list.

## See also

- {doc}`constraints` — constraint framework and run-time mode customisation
- {func}`~ad_hoc_diffractometer.reference.incidence_angle`
- {func}`~ad_hoc_diffractometer.reference.exit_angle`
- {func}`~ad_hoc_diffractometer.reference.psi_angle`
- {func}`~ad_hoc_diffractometer.reference.naz_angle`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
