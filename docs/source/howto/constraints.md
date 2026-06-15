(howto-constraints)=
# Work with Constraints and Diffraction Modes

This guide explains the constraint system in detail: how modes are defined,
how to customize constraint values, how to build entirely new modes, and how
to work with the `extras` dict for advanced modes.

## Background: degrees of freedom

Specifying a reflection (h, k, l) provides exactly **3 equations** on the
motor angles (the three components of the scattering vector Q).
A geometry with N real motor axes therefore has **N − 3 free parameters**
after the Bragg condition is satisfied.  Each free parameter must be resolved
by exactly one constraint.

| Geometry family | N | N − 3 |
|---|---|---|
| {ref}`geometry-fourcv`, {ref}`geometry-fourch`, {ref}`geometry-zaxis`, {ref}`geometry-s2d2` | 4 | 1 |
| {ref}`geometry-fivec` | 5 | 2 |
| {ref}`geometry-psic`, {ref}`geometry-kappa6c`, {ref}`geometry-sixc` | 6 | 3 |

```python
import ad_hoc_diffractometer as ahd

print(ahd.make_geometry("fourcv").free_dof_after_bragg)   # 1
print(ahd.make_geometry("psic").free_dof_after_bragg)     # 3
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
vector n̂ (surface normal, polarization axis, etc.).
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
g = ahd.make_geometry("fourcv")
g.mode_name = "bisecting"   # 1 BisectConstraint (N-3=1)
solutions = g.forward(1, 0, 0)
```

## Set a constraint value at run time

A constraint value is **constant for the duration of a single
{func}`~ad_hoc_diffractometer.forward.compute_forward` call**.
The {class}`~ad_hoc_diffractometer.mode.ConstraintSet` object persists on the
geometry until explicitly replaced — there is no need to reassign it if the
value does not change between calls.

### Why constraint values are immutable

Every constraint (`SampleConstraint`, `DetectorConstraint`,
`ReferenceConstraint`) exposes `.value` as a **read-only property** — there
is no setter and no `set_value()` method.  This is deliberate:

- Constraints are *values*, not mutable state, so `ConstraintSet.to_dict()`
  / `from_dict()` round-trip cleanly and serialized modes always describe
  the run-time state.
- The YAML file under
  `src/ad_hoc_diffractometer/geometries/` is the **single source of
  truth** for default values; run-time overrides happen by replacing the
  containing `ConstraintSet`, never by mutating the constraint object.
- Two callers holding references to the same `ConstraintSet` (e.g. the
  active mode and a snapshot taken before a scan) cannot surprise each
  other with hidden value changes.

To override a default value you therefore produce a **new**
`ConstraintSet`.  There are two ways to do that — pick whichever is more
readable for your case.

### Use `with_constraint_values` for a one-call override

`ConstraintSet.with_constraint_values(**updates)` returns a fresh
`ConstraintSet` with the named constraint values replaced.  Each keyword
argument names a constraint by its `.name` attribute (a stage name for
sample / detector constraints, a reference name like `alpha_i` /
`beta_out` / `psi` / `a_eq_b` for reference constraints).  Constraint
order, `computed`, `extras`, and `cut_points` are preserved.

```python
import ad_hoc_diffractometer as ahd

g = ahd.make_geometry("psic")

# Multiple values at once (psic B3 mode: chi, phi, and the alpha_i target):
g.modes["fixed_alpha_i_fixed_chi_fixed_phi"] = (
    g.modes["fixed_alpha_i_fixed_chi_fixed_phi"]
    .with_constraint_values(chi=15.0, phi=30.0, alpha_i=5.0)
)

# Single value (psic fixed_chi_vertical: default chi=90° → 45°):
g.modes["fixed_chi_vertical"] = (
    g.modes["fixed_chi_vertical"].with_constraint_values(chi=45.0)
)
```

Unknown keys (typos, names that do not appear in the set) raise
`KeyError` listing every unrecognised key at once so a multi-typo edit
can be fixed in one pass.  `BisectConstraint` is relational (it has no
scalar value) and is invisible to this method — any kwarg targeting a
bisect raises `KeyError`.

### Rebuild the whole `ConstraintSet`

When you want to change which constraints appear (not just their
values), construct a new `ConstraintSet` directly:

```python
from ad_hoc_diffractometer import ConstraintSet, SampleConstraint

g = ahd.make_geometry("fourcv")

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

A detector-stage value works identically — pass the new
`DetectorConstraint` in the list:

```python
from ad_hoc_diffractometer import DetectorConstraint

# psic fixed_delta-style mode: change the detector pin to a non-zero value
g.modes["my_fixed_delta"] = ConstraintSet(
    [
        BisectConstraint("eta", "delta"),
        SampleConstraint("mu", 0.0),
        DetectorConstraint("nu", 5.0),       # was 0.0; now 5.0
    ],
    computed=["eta", "chi", "phi", "delta"],
)
```

For a reference-target value (surface modes), pass a new
`ReferenceConstraint`:

```python
from ad_hoc_diffractometer import ReferenceConstraint

# psic B3 mode: pin the incidence angle alpha_i at 5° instead of the
# YAML default 0°.
g.modes["fixed_alpha_i_fixed_chi_fixed_phi"] = ConstraintSet(
    [
        SampleConstraint("chi", 0.0),
        SampleConstraint("phi", 0.0),
        ReferenceConstraint("alpha_i", 5.0),  # was 0.0; now 5.0
    ],
    computed=g.modes["fixed_alpha_i_fixed_chi_fixed_phi"].computed,
    extras=dict(g.modes["fixed_alpha_i_fixed_chi_fixed_phi"].extras),
)
```

The `computed` field on {class}`~ad_hoc_diffractometer.mode.ConstraintSet` is
informational (documents which stages the solver computes) and does not
affect the calculation.

## Build a multi-constraint mode (six-circle)

For psic (N − 3 = 3), three constraints are needed:

```python
from ad_hoc_diffractometer import (
    ConstraintSet, BisectConstraint, SampleConstraint, DetectorConstraint
)

g = ahd.make_geometry("psic")

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
set on the geometry (e.g. ``azimuth`` for psi modes,
``surface_normal`` for surface modes) — they are considered stubs until
the prerequisite is met:

```python
g = ahd.make_geometry("fourcv")
g.mode_name = "fixed_psi"

# Without azimuth: not implemented
print(g.modes["fixed_psi"].is_implemented(g))  # False

# With azimuth: implemented
g.azimuth = (0, 0, 1)
print(g.modes["fixed_psi"].is_implemented(g))  # True
```

## Serialisation

{class}`~ad_hoc_diffractometer.mode.ConstraintSet` round-trips through `to_dict()` / `from_dict()`:

```python
import json
from ad_hoc_diffractometer import AdHocDiffractometer

g = ahd.make_geometry("fourcv")
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
