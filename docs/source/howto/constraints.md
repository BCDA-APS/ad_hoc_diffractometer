(howto-constraints)=
# Work with Constraints and Diffraction Modes

This guide explains the constraint system in detail: how modes are defined,
how to customise constraint values, how to build entirely new modes, and how
to work with the `extras` dict for advanced modes.

## Background: degrees of freedom

Specifying a reflection (h, k, l) provides exactly **3 equations** on the
motor angles (the three components of the scattering vector Q).
A geometry with N real motor axes therefore has **N − 3 free parameters**
after the Bragg condition is satisfied.  Each free parameter must be resolved
by exactly one constraint.

| Geometry family | N | N − 3 |
|---|---|---|
| {func}`~ad_hoc_diffractometer.presets.fourcv`, {func}`~ad_hoc_diffractometer.presets.fourch`, {func}`~ad_hoc_diffractometer.presets.zaxis`, {func}`~ad_hoc_diffractometer.presets.s2d2` | 4 | 1 |
| {func}`~ad_hoc_diffractometer.presets.fivec` | 5 | 2 |
| {func}`~ad_hoc_diffractometer.presets.psic`, {func}`~ad_hoc_diffractometer.presets.kappa6c`, {func}`~ad_hoc_diffractometer.presets.sixc` | 6 | 3 |

```python
import ad_hoc_diffractometer as ahd

print(ahd.presets.fourcv().free_dof_after_bragg)   # 1
print(ahd.presets.psic().free_dof_after_bragg)     # 3
```

## Constraint categories

Every constraint belongs to one of three categories:

**Sample constraint** — fixes one sample stage at a declared value, or
declares the bisecting relational condition:

```python
from ad_hoc_diffractometer import SampleConstraint, BisectConstraint

SampleConstraint("chi", 90.0)          # chi fixed at 90°
SampleConstraint("omega", 0.0)         # omega fixed at 0°
BisectConstraint("omega", "ttheta")    # omega = ttheta / 2
BisectConstraint("eta", "delta")       # eta = delta / 2  (psic)
```

**Detector constraint** — fixes one detector stage at a declared value, or
constrains the ``"qaz"`` pseudo-angle from You (1999) eq. 18
(``tan(qaz) = tan(delta) / sin(nu)``):

```python
from ad_hoc_diffractometer import DetectorConstraint

DetectorConstraint("nu", 0.0)     # nu fixed at 0°
DetectorConstraint("gamma", 0.0)  # gamma fixed at 0°
DetectorConstraint("qaz", 90.0)   # Q confined to the vertical plane
```

``qaz = 90°`` constrains scattering to the vertical plane;
``qaz = 0°`` to the horizontal plane.
This is implemented for all geometries with two or more detector stages
(psic, kappa6c) and used by the ``lifting_detector_*`` mode family.

**Reference constraint** — expresses a condition between Q and a reference
vector n̂ (surface normal, polarisation axis, etc.).
The incidence/exit-angle constraints (``alpha_i``, ``beta_out``,
``a_eq_b``) are implemented when ``surface_normal`` is set;
``psi`` and ``naz`` are not yet implemented as forward constraints.
See {doc}`surface`:

```python
from ad_hoc_diffractometer import ReferenceConstraint

ReferenceConstraint("psi", 90.0)       # azimuthal angle of n̂ about Q
ReferenceConstraint("alpha_i", 5.0)    # incidence angle
ReferenceConstraint("a_eq_b", True)    # alpha_i = beta_out (symmetric)
```

**Rules:** at most one {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`,
at most one {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`,
remainder must be {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
or {class}`~ad_hoc_diffractometer.mode.BisectConstraint`.
Total must equal N − 3.

## Design principles

The constraint system is built on three key decisions from the #122 planning
discussion:

1. **Constraints are geometry-agnostic.**  A {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
   is defined without reference to a specific geometry; validation against
   actual DOF count and stage names happens at solve time
   (`is_fully_constrained(g)` and `is_implemented(g)`).

2. **At most one detector constraint and at most one reference constraint.**
   Fixing more than one detector angle over-constrains the scattered beam
   direction.  More than one reference constraint would also over-constrain
   the problem.  Any number of sample constraints (fixed-angle or bisect) are
   allowed, subject to the total equalling N − 3.

3. **The bisect condition is relational, not heuristic.**
   {class}`~ad_hoc_diffractometer.mode.BisectConstraint` names both stages
   explicitly (``sample_stage = detector_stage / 2``).  No geometric
   heuristics are used to infer which stage is "co-axial" with which.

These principles make it possible to define valid modes programmatically
at run time without any knowledge of geometry internals.

## Use a factory-defined mode

```python
g = ahd.presets.fourcv()
g.mode_name = "bisecting"   # 1 BisectConstraint (N-3=1)
solutions = g.forward(1, 0, 0)
```

## Set a constraint value at run time

A constraint value is **constant for the duration of a single
{func}`~ad_hoc_diffractometer.forward.compute_forward` call**.
The {class}`~ad_hoc_diffractometer.mode.ConstraintSet` object persists on the
geometry until explicitly replaced — there is no need to reassign it if the
value does not change between calls.

```python
from ad_hoc_diffractometer import ConstraintSet, SampleConstraint

g = ahd.presets.fourcv()

# Set once — all subsequent forward() calls use chi = 45°
g.modes["my_chi"] = ConstraintSet([SampleConstraint("chi", 45.0)])
g.mode_name = "my_chi"

sols_100 = g.forward(1, 0, 0)   # chi = 45°
sols_010 = g.forward(0, 1, 0)   # chi = 45° (same constraint, no reassignment needed)
sols_111 = g.forward(1, 1, 1)   # chi = 45°

# Only reassign when the value changes
g.modes["my_chi"] = ConstraintSet([SampleConstraint("chi", 60.0)])
sols_new = g.forward(1, 0, 0)   # chi = 60° from here on
```

The `"fixed_chi"` mode in the demo geometry has a default chi = 90°.
Replacing the {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
in `g.modes["fixed_chi"]` changes the value for all subsequent calls
until it is replaced again.

The `computed` field on {class}`~ad_hoc_diffractometer.mode.ConstraintSet` is
informational (documents which stages the solver computes) and does not
affect the calculation.

## Build a multi-constraint mode (six-circle)

For psic (N − 3 = 3), three constraints are needed:

```python
from ad_hoc_diffractometer import (
    ConstraintSet, BisectConstraint, SampleConstraint, DetectorConstraint
)

g = ahd.presets.psic()

# Custom bisecting mode with mu=5° (non-zero mu)
g.modes["bisecting_mu5"] = ConstraintSet(
    [
        BisectConstraint("eta", "delta"),   # sample: eta = delta/2
        SampleConstraint("mu", 5.0),        # sample: mu fixed at 5°
        DetectorConstraint("nu", 0.0),      # detector: nu fixed at 0°
    ],
    computed=["eta", "chi", "phi", "delta"],
)
g.mode_name = "bisecting_mu5"
```

## Inspect a mode

```python
cs = g.modes["bisecting"]

# All constraints in definition order
print(cs.constraints)

# Stage names held fixed (derived from constraints, not a separate dict)
print(cs.constant_stages)    # e.g. ['eta', 'mu', 'nu']

# Stage names computed by the solver
print(cs.computed)           # e.g. ['eta', 'chi', 'phi', 'delta']

# Bisect pair (sample stage, detector stage)
print(cs.bisect_stages())    # e.g. ('eta', 'delta')

# DOF check
print(cs.is_fully_constrained(g))   # True if len(constraints) == N-3

# Solver availability
print(cs.is_implemented(g))         # True if all constraints have solvers
```

## The extras dict

Some modes require additional input beyond (h, k, l), or compute additional
output quantities alongside the motor angles.
These are declared in the `extras` dict using {data}`~ad_hoc_diffractometer.mode.REQUIRED`
and {data}`~ad_hoc_diffractometer.mode.OPTIONAL` sentinels:

```python
from ad_hoc_diffractometer import REQUIRED, OPTIONAL, ConstraintSet, ReferenceConstraint

# fixed_psi mode on fourcv: requires a reference vector n̂
cs = g.modes["fixed_psi"]
print(cs.extras)
# {'n_hat': <REQUIRED>, 'psi': None}
# n_hat must be supplied; psi is an output populated by the solver
```

`REQUIRED` marks inputs that **must** be supplied before calling `forward()`.
`None` marks output slots populated by the solver.
`OPTIONAL` marks inputs with a sensible default.

## Stub modes

Modes whose constraint patterns do not yet have a solver implementation
return `False` from `is_implemented()` and raise `NotImplementedError`
when `forward()` is called.  Some modes require a prerequisite to be
set on the geometry (e.g. ``azimuthal_reference`` for psi modes,
``surface_normal`` for surface modes) — they are considered stubs until
the prerequisite is met:

```python
g = ahd.presets.fourcv()
g.mode_name = "fixed_psi"

# Without azimuthal_reference: not implemented
print(g.modes["fixed_psi"].is_implemented(g))  # False

# With azimuthal_reference: implemented
g.azimuthal_reference = (0, 0, 1)
print(g.modes["fixed_psi"].is_implemented(g))  # True
```

## Serialisation

{class}`~ad_hoc_diffractometer.mode.ConstraintSet` round-trips through `to_dict()` / `from_dict()`:

```python
import json
from ad_hoc_diffractometer import AdHocDiffractometer

g = ahd.presets.fourcv()
d = g.to_dict()
json.dumps(d)   # JSON-serialisable

g2 = AdHocDiffractometer.from_dict(d)
print(g2.modes.keys())    # same modes as g
print(g2.mode_name)       # 'bisecting'
```

## Custom constraint protocol

Any object satisfying the {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
interface can be added to a {class}`~ad_hoc_diffractometer.mode.ModeDict`.
A minimal custom constraint must implement:

- `name: str` — a unique identifier
- `category: str` — `"sample"`, `"detector"`, or `"reference"`
- `extras: dict` — input/output parameters
- `evaluate(angles, geometry) -> float` — residual (0 = satisfied)
- `is_satisfied(angles, geometry, tol) -> bool`
- `is_implemented(geometry) -> bool`
- `to_dict() / from_dict()` — serialisation

The numeric fallback solver in `forward.py` handles any custom constraint
that has `is_implemented()` returning `True`, even without an analytic solver.

## See also

- {doc}`modes` — everyday mode switching
- {class}`~ad_hoc_diffractometer.mode.ConstraintSet`
- {class}`~ad_hoc_diffractometer.mode.BisectConstraint`
- {class}`~ad_hoc_diffractometer.mode.SampleConstraint`
- {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`
- {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint`
- {data}`~ad_hoc_diffractometer.mode.REQUIRED`
- {data}`~ad_hoc_diffractometer.mode.OPTIONAL`
- {exc}`~ad_hoc_diffractometer.mode.EwaldSphereViolation`
- {exc}`~ad_hoc_diffractometer.mode.ConstraintViolation`
