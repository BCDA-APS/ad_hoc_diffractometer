(declarative-geometry-schema)=
# Declarative Geometry Schema

Authoritative reference for the YAML schema used by
{func}`ad_hoc_diffractometer.geometry_loader.load_geometry_file` and
{func}`ad_hoc_diffractometer.geometry_loader.register_geometry_file`.
See the {ref}`how-to <howto-custom-geometry>` for a task-oriented
walkthrough.

The 10 YAML files shipped under
`src/ad_hoc_diffractometer/geometries/` are **demonstrations** of this
schema, not authoritative descriptions of any production
diffractometer.  Real instruments at real beamlines should ship their
own YAML — copy the closest demo file and adapt it.

---

## Top-level structure

Every declarative geometry file is a single YAML document (one mapping
at the top level) with the following keys.  The strict
{ref}`unknown-key policy <decl-strict-keys>` rejects any key not listed
here.

| Key | Required | Type | Purpose |
|---|---|---|---|
| `ad_hoc_diffractometer_geometry` | yes | mapping | Document marker (see below) |
| `name` | yes | non-empty string | Registry key (must equal filename stem for files in the package directory) |
| `documentation` | recommended | string | Free-form prose; the **first line is the one-line summary** (Python-docstring convention) |
| `basis` | optional | string or mapping | Coordinate basis (see {ref}`decl-basis`) |
| `parameters` | optional | mapping | Caller-tunable defaults; currently only `alpha_deg` |
| `kappa_chi_eq` | conditional | string or 3-vector | Equivalent-Eulerian chi axis; required for kappa geometries |
| `stages` | yes | non-empty list | Rotary stages (see {ref}`decl-stages`) |
| `modes` | yes | non-empty mapping | Diffraction modes (see {ref}`decl-modes`) |

---

(decl-marker)=
## Document marker

Every file must begin with the marker that identifies it as an
`ad_hoc_diffractometer` geometry declaration:

```yaml
ad_hoc_diffractometer_geometry:
    schema_revision: 1
```

The presence of the `ad_hoc_diffractometer_geometry` top-level key is
the unambiguous signal that the document is a geometry declaration.
The `schema_revision` integer selects which version of *this schema*
the file conforms to.  It is a fixed property of the schema, not a
per-file edit counter — users should treat the value verbatim and
only change it when migrating the file to a different declarative-
geometry schema revision.

| Loader constant | Value | Meaning |
|---|---|---|
| {data}`~ad_hoc_diffractometer.geometry_loader.KIND_KEY` | `"ad_hoc_diffractometer_geometry"` | The marker key name |
| {data}`~ad_hoc_diffractometer.geometry_loader.SUPPORTED_REVISIONS` | `frozenset({1})` | Schema revisions this build can parse |
| {data}`~ad_hoc_diffractometer.geometry_loader.CURRENT_REVISION` | `1` | Revision newly-authored YAML files SHOULD declare |

The schema revision is independent of the package version; schema
revisions are deliberate editorial events that change far less often
than package releases.

```{important}
When a string passed to
{func}`~ad_hoc_diffractometer.geometry_loader.load_geometry_file`
**lacks** the marker, the loader treats the string as a filesystem
path rather than as YAML text.  When the marker is **present** but
malformed (wrong type for the value, unsupported `schema_revision`,
unknown sub-key, …) the loader raises {class}`ValueError` immediately and
does **not** retry as a path — once the document declares itself a
geometry, the loader honors that declaration by reporting what is
wrong with it.
```

---

(decl-basis)=
## `basis`

Optional.  Specifies the coordinate basis used to resolve physical
direction names (`vertical`, `longitudinal`, `transverse`) in stage
axis fields and elsewhere.

Three accepted forms:

1. **Shorthand string** — one of:
   - `"BL"` → {data}`~ad_hoc_diffractometer.factories.BASIS_BL`
     (Busing & Levy 1967: vertical=Z, longitudinal=Y, transverse=X)
   - `"YOU"` → {data}`~ad_hoc_diffractometer.factories.BASIS_YOU`
     (You 1999: vertical=X, longitudinal=Y, transverse=Z)
   - `"DEFAULT"` → {data}`~ad_hoc_diffractometer.factories.BASIS_DEFAULT`
     (vertical=Y, longitudinal=Z, transverse=X)

2. **Mapping of signed-axis strings**

   ```yaml
   basis:
       vertical: '+z'
       longitudinal: '+y'
       transverse: '+x'
   ```

3. **Mapping of numeric unit vectors** — axis-aligned:

   ```yaml
   basis:
       vertical:     [0.0, 0.0, 1.0]
       longitudinal: [0.0, 1.0, 0.0]
       transverse:   [1.0, 0.0, 0.0]
   ```

   …or any orthonormal triple, including non-axis-aligned bases.  For
   example, a 45° rotation of the (vertical, longitudinal) pair within
   the X-Y plane:

   ```yaml
   basis:
       vertical:     [0.7071, 0.7071, 0.0]
       longitudinal: [-0.7071, 0.7071, 0.0]
       transverse:   [0.0, 0.0, 1.0]
   ```

The loader validates orthonormality (each vector has unit length and
the three are mutually orthogonal, both within tolerance `1e-9`).  A
defensive `|det([v, l, t])| ≈ 1` cross-check rejects degenerate or
co-planar inputs.  Either chirality of the
`(vertical, longitudinal, transverse)` ordering is accepted because
both BL and YOU conventions are physically meaningful right-handed
frames in their own canonical orderings.

### Precedence

A geometry's basis is resolved in this order (first match wins):

1. **Caller-supplied `basis=`** keyword to
   {func}`~ad_hoc_diffractometer.factories.make_geometry`,
   {func}`~ad_hoc_diffractometer.geometry_loader.load_geometry_file`, or
   {func}`~ad_hoc_diffractometer.geometry_loader.register_geometry_file`.
2. **YAML file's `basis:`** key (any of the three forms above).
3. **{data}`~ad_hoc_diffractometer.factories.BASIS_DEFAULT`** — the
   neutral fallback.

Omitting the YAML key emits a {class}`UserWarning` ("does not declare a
basis; using BASIS_DEFAULT …") so accidental omission is visible.  Real
instruments should declare `basis:` explicitly.

#### Worked example of the precedence rule

```python
import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.factories import BASIS_BL, BASIS_YOU

# (1) Caller wins over YAML default.
g = ahd.make_geometry("fourcv", basis=BASIS_YOU)
assert (g.basis["vertical"] == BASIS_YOU["vertical"]).all()

# (2) No kwarg → YAML's `basis: BL` wins.
g = ahd.make_geometry("fourcv")
assert (g.basis["vertical"] == BASIS_BL["vertical"]).all()
```

---

(decl-parameters)=
## `parameters`

Optional mapping of caller-tunable parameters with their defaults.
Currently the only recognized key is:

| Key | Type | Purpose |
|---|---|---|
| `alpha_deg` | float | Kappa tilt angle (degrees); used by `kappa_eulerian` axis specs |

When `parameters.alpha_deg` is declared, the loader also requires
{ref}`decl-kappa-chi-eq` and synthesizes a
{class}`~ad_hoc_diffractometer.kappa.KappaPseudoAngleConvention` from
the stages named `komega`, `kappa`, `kphi`.

```{note}
A caller can override `alpha_deg` at call time:
``ahd.make_geometry("kappa4cv", alpha_deg=55)``.
```

---

(decl-kappa-chi-eq)=
## `kappa_chi_eq`

Required when the file declares `parameters.alpha_deg`.  Specifies the
equivalent-Eulerian chi axis used by `kappa_axis_from_eulerian`.
Accepts the same forms as a stage axis: a signed-axis string
(`'+vertical'`, `'-x'`, …) or a length-3 numeric vector.

---

(decl-stages)=
## `stages`

A non-empty list.  Each entry is a mapping with **all four** keys:

| Key | Type | Purpose |
|---|---|---|
| `name` | string | Stage name (e.g. `omega`, `kappa`, `mu`) |
| `axis` | see below | Rotation axis (signed) |
| `parent` | string or null | Name of the stage this one sits on, or `null` for the lab frame |
| `role` | string | `"sample"`, `"detector"`, or any other role label |

### Axis grammar

Three accepted forms:

1. **Signed physical-direction string** — resolved against the active
   basis:
   ```yaml
   axis: -transverse
   axis: +longitudinal
   ```
2. **Signed Cartesian string** — resolved directly against `XHAT`,
   `YHAT`, `ZHAT`:
   ```yaml
   axis: '+x'
   axis: '-z'
   ```
3. **Length-3 numeric vector** (used as-is, sign included; no automatic
   normalisation other than what the
   {class}`~ad_hoc_diffractometer.stage.Stage` constructor performs):
   ```yaml
   axis: [0.0, 0.5, 0.866]
   ```
4. **Kappa-tilt mapping** — only meaningful when
   `parameters.alpha_deg` and `kappa_chi_eq` are declared.  The
   single argument is the outer (komega) axis:
   ```yaml
   axis: {kappa_eulerian: '+transverse'}
   ```
   The loader resolves this via
   {func}`~ad_hoc_diffractometer.kappa.kappa_axis_from_eulerian`.

```{important}
The sign of the axis vector encodes handedness: ``+nHat`` is
right-handed rotation about ``nHat``; ``-nHat`` is left-handed.  The
``rotation: left/right`` shorthand discussed in early issue drafts is
**not** part of the schema — use signed-axis strings instead.
```

---

(decl-modes)=
## `modes`

A non-empty mapping of mode-name → mode-spec.  Each mode-spec is a
mapping with these keys:

| Key | Required | Type | Purpose |
|---|---|---|---|
| `default` | optional | bool | When `true`, this is the geometry's default mode.  At most one mode per file may set this. |
| `constraints` | yes | list | Constraint specs (may be empty) |
| `computed` | optional | list of strings | Stage names the solver writes |
| `extras` | optional | mapping | Additional inputs/outputs (see below) |
| `cut_points` | optional | mapping | Per-stage SPEC `#G4` cut-points |

### Constraint specs

Each entry of `constraints:` is a mapping with a `type:` key selecting
one of five constraint classes.  Each type accepts a fixed set of
sibling keys; the {ref}`unknown-key policy <decl-strict-keys>` rejects
any others.

| `type` | Required keys | Constraint class |
|---|---|---|
| `bisect` | `stage1`, `stage2` | {class}`~ad_hoc_diffractometer.mode.BisectConstraint` |
| `virtual_bisect` | `stage1`, `stage2` | {class}`~ad_hoc_diffractometer.mode.VirtualBisectConstraint` |
| `sample` | `stage`, `value` | {class}`~ad_hoc_diffractometer.mode.SampleConstraint` |
| `detector` | `stage`, `value` | {class}`~ad_hoc_diffractometer.mode.DetectorConstraint` |
| `reference` | `name`, `value` | {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint` |

For `reference`, `name` is one of `"psi"`, `"alpha_i"`, `"beta_out"`,
`"a_eq_b"`, `"naz"`, `"omega"`; `value` is a float for the angular
constraints and the literal `true` for `a_eq_b`.  The `"omega"` name
selects the SPEC `OMEGA` pseudo-angle (the angle between Q and the
plane of the chi circle, SPEC `psic` `def OMEGA 'Q[6]'`); it applies
to psic-family geometries (those with a sample stage named `chi`) and
does not require any reference vector.

### `extras` and the `REQUIRED` sentinel

`extras` accepts arbitrary key→value pairs.  Two special values:

- The literal string `REQUIRED` is mapped to the
  {data}`~ad_hoc_diffractometer.mode.REQUIRED` sentinel; the solver
  checks for the sentinel before running and raises if the caller
  didn't supply a real value.
- YAML `null` → Python `None`.  Used for output slots that the solver
  populates after `forward()` runs.

```yaml
extras:
    n_hat: REQUIRED        # caller must supply; solver refuses if absent
    psi: null              # output slot; solver writes to it
```

---

(decl-strict-keys)=
## Strict unknown-key policy

The loader rejects any key it does not recognize at every nesting
level: top level, marker mapping, `parameters` mapping, each stage
entry, each mode entry, and each constraint spec.  The error names
the offending key, the containing context, and the accepted-key set,
so misspellings (`documantation:` for `documentation:`) surface
immediately.

Future schema revisions may relax this policy if a use case justifies
extension keys; revision 1 is strict.

---

## What `load_geometry_file()` accepts

{func}`~ad_hoc_diffractometer.geometry_loader.load_geometry_file` accepts
either a filesystem path **or** YAML text held in a Python string.  The
two cases work as follows:

- **A filesystem path** — passed as a {class}`pathlib.Path`, any
  {class}`os.PathLike`, or a path string.  The loader reads the file
  and parses its contents.  {class}`FileNotFoundError` if the file
  does not exist.
- **YAML text in memory** — passed as a string that begins (after
  optional leading whitespace) with the
  ``ad_hoc_diffractometer_geometry:`` schema marker.  The loader
  parses the string directly without touching the filesystem.

A string that does not contain the marker is interpreted as a path,
not as YAML text.  If neither attempt succeeds, the
{class}`FileNotFoundError` message names both ("first attempted as
YAML text … then attempted as a path …") so you can tell what went
wrong.

### Examples

```python
import pathlib
import ad_hoc_diffractometer as ahd

# (1) Path object
g = ahd.load_geometry_file(pathlib.Path("/lab/diffr.yml"))

# (2) Path string
g = ahd.load_geometry_file("/lab/diffr.yml")

# (3) YAML text in memory
g = ahd.load_geometry_file("""
ad_hoc_diffractometer_geometry:
    schema_revision: 1
name: my_diffr
basis: BL
stages:
    - {name: omega,  axis: -transverse, parent: null,  role: sample}
    - {name: ttheta, axis: -transverse, parent: null,  role: detector}
modes:
    bisecting:
        default: true
        constraints:
            - {type: bisect, stage1: omega, stage2: ttheta}
        computed: [omega, ttheta]
""")
```

---

## Run-time registration

To make a user-authored YAML file available via
{func}`~ad_hoc_diffractometer.factories.make_geometry`, register it once:

```python
import ad_hoc_diffractometer as ahd

name = ahd.register_geometry_file("/lab/diffr.yml")
g = ahd.make_geometry(name)        # name is the filename stem
g = ahd.make_geometry(name, basis=ahd.factories.BASIS_YOU)  # override
```

A second call with a colliding name raises {class}`ValueError`; pass
``name="my_lab_diffr"`` to register under a different name.  The
registry stores a wrapper that re-reads the file on each call, so
on-disk edits are picked up without re-registration.

---

## Migrating between schema revisions

No migrations are required at this time.  This section will be
populated when a new schema revision is introduced.

---

(decl-json-schema)=
## JSON Schema (revision 1)

A machine-readable JSON Schema is shipped alongside the package at
`src/ad_hoc_diffractometer/geometries/schema.json` and is accessible at
runtime via:

```python
from ad_hoc_diffractometer.geometry_loader import get_schema, get_schema_text

schema = get_schema()      # parsed dict
text = get_schema_text()   # raw JSON string
```

The loader **does not** consume this schema at runtime — its
hand-written validator enforces the same rules in pure Python and
produces context-rich error messages.  The JSON Schema exists for
three audiences:

1. **Editor tooling.**  VS Code, JetBrains, and other editors can
   validate `.yml` files against a JSON Schema for autocomplete and
   live linting.  Point your editor at the file above (or at the
   ``$id`` URL once published).
2. **Documentation.**  The page you are reading is generated from the
   same schema source.
3. **Human readers** who want a single dense source of truth.

The full schema:

```{literalinclude} ../../../src/ad_hoc_diffractometer/geometries/schema.json
:language: json
```

---

## See also

- {ref}`How-to: build a custom geometry <howto-custom-geometry>`
- {func}`~ad_hoc_diffractometer.geometry_loader.load_geometry_file`
- {func}`~ad_hoc_diffractometer.geometry_loader.register_geometry_file`
- {func}`~ad_hoc_diffractometer.factories.make_geometry`
- {data}`~ad_hoc_diffractometer.factories.BASIS_BL`
- {data}`~ad_hoc_diffractometer.factories.BASIS_YOU`
- {data}`~ad_hoc_diffractometer.factories.BASIS_DEFAULT`
