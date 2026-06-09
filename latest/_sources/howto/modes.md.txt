(howto-modes)=
# Switch Diffraction Modes

A **diffraction mode** is a {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
that describes which motor stages are free, which are held at declared values,
and which are related to others (e.g. bisecting).
See [Concepts](../concepts.md) for background, and
{doc}`constraints` for the full constraint framework.

## List available modes

```python
import ad_hoc_diffractometer as ahd

g = ahd.make_geometry("fourcv")
print(list(g.modes.keys()))
# ['bisecting', 'fixed_chi', 'fixed_phi', 'fixed_omega',
#  'fixed_psi', 'double_diffraction']
```

## Get and set the active mode

```python
# The default mode is set by the demo geometry
print(g.mode_name)   # 'bisecting'

# Switch to a different mode
g.mode_name = "fixed_chi"
print(g.mode_name)   # 'fixed_chi'

# Clear the active mode (all stages free — forward() will raise)
g.mode_name = None
```

## Mode reference

### Bisecting mode

The bisecting {class}`~ad_hoc_diffractometer.mode.BisectConstraint` drives the
co-axial sample stage to half the detector angle, placing the sample
symmetrically between the incident and diffracted beams.

```python
g.mode_name = "bisecting"
solutions = g.forward(1, 0, 0)
# omega = ttheta / 2 in every solution
```

### Fixed-angle modes

Fixed-angle modes hold one sample stage at its **declared constraint value**
during a single `forward()` call.  The value is part of the {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
definition and is constant only for the duration of that call.

```python
# fixed_chi: demo geometry holds chi at 90°
g.mode_name = "fixed_chi"
solutions = g.forward(1, 0, 0)
# sol["chi"] == 90.0 for every solution
```

To use a **different value**, the simplest call is
{meth}`~ad_hoc_diffractometer.mode.ConstraintSet.with_constraint_values`,
which returns a fresh `ConstraintSet` with the named constraint values
replaced:

```python
# Single value: chi = 45° on the demo fixed_chi mode
g.modes["fixed_chi"] = g.modes["fixed_chi"].with_constraint_values(chi=45.0)
g.mode_name = "fixed_chi"
sols_45 = g.forward(1, 0, 0)   # chi = 45°

# Several values at once (e.g. psic B3 mode: two sample stages plus
# the alpha_i target — three pinned values in one call):
g.modes["fixed_alpha_i_fixed_chi_fixed_phi"] = (
    g.modes["fixed_alpha_i_fixed_chi_fixed_phi"]
    .with_constraint_values(chi=15.0, phi=30.0, alpha_i=5.0)
)
```

For modes where you want to change *which* constraints appear (not just
their values), build a new {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
directly:

```python
from ad_hoc_diffractometer import ConstraintSet, SampleConstraint

g.modes["my_chi"] = ConstraintSet([SampleConstraint("chi", 60.0)])
g.mode_name = "my_chi"
```

See {doc}`constraints` for the full run-time pattern, including how to
override detector and reference-target values and the rationale for the
immutable-constraint design.

## Change a fixed-axis value

Yes, you can override the YAML default for any `fixed_*` mode at run
time.  Use
{meth}`~ad_hoc_diffractometer.mode.ConstraintSet.with_constraint_values`
for a one-call value swap (single or multiple values), or rebuild the
{class}`~ad_hoc_diffractometer.mode.ConstraintSet` from scratch when
you also need to change which constraints appear.  Worked examples for
sample, detector, and reference-target overrides are in
{ref}`howto-constraints`.

### fixed_omega

Holds `omega` at 0° (any geometry where omega is a sample stage).

```python
g.mode_name = "fixed_omega"
solutions = g.forward(1, 0, 0)
# sol["omega"] == 0.0 for every solution
```

## Inspect a mode's constraints

```python
g = ahd.make_geometry("fourcv")
cs = g.modes["fixed_chi"]

print(cs)
# ConstraintSet([SampleConstraint('chi', 90.0)])

print(cs.computed)
# ['omega', 'phi', 'ttheta']

print(cs.constant_stages)
# ['chi']

print(cs.is_fully_constrained(g))
# True  (N - 3 = 1 constraint for N = 4 DOF)

print(cs.is_implemented(g))
# True
```

## Check if a mode is implemented

Some modes require a prerequisite on the geometry.  For example,
`fixed_psi` requires ``g.azimuthal_reference`` to be set:

```python
g.mode_name = "fixed_psi"
print(g.modes["fixed_psi"].is_implemented(g))  # False (no azimuthal_reference)

g.azimuthal_reference = (0, 0, 1)
print(g.modes["fixed_psi"].is_implemented(g))  # True

# forward() raises NotImplementedError for unimplemented modes
```

## Clear the active mode

```python
g.mode_name = None   # all stages free; forward() raises NotImplementedError
```

## See also

- {doc}`constraints` — full constraint framework: DOF rule, custom modes, extras
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`
- {class}`~ad_hoc_diffractometer.mode.ModeDict`
- [Concepts — Diffraction modes](../concepts.md)
