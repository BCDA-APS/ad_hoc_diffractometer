# AGENTS.md — ad_hoc_diffractometer

This file provides context for AI coding agents working on this project.

---

## Attribution

Always sign your work. Every commit message, pull request description, and
inline comment authored by an AI agent must end with:

```
Contributed by: <AgentName> (<model-id>)
```

Example: `Contributed by: OpenCode (argo/claudesonnet46)`

---

## Issue and PR workflow

1. **Before writing any code**, set the issue status to "In progress" on the
   GitHub project board:
   ```bash
   # Find the project item ID for the issue, then update its status field
   gh api graphql -f query='...'   # see project board queries below
   ```

2. **Open a PR** with a body that includes a closing reference as a bullet,
   followed by any remarks and the agent/model signature:
   ```
   - closes #N

   <optional additional context>

   Contributed by: OpenCode (argo/claudesonnet46)
   ```
   The `closes #N` keyword triggers GitHub automation: when the PR is merged
   the issue is closed and the project card moves to "Done" automatically.

   Every PR must also have its **milestone** and **project board** set,
   and its project board status must be set to **"In Progress"**.
   Issues and PRs may belong to **more than one project board** — add the
   card to every board that tracks the work.

   | Domain | Milestone | Project ID |
   |---|---|---|
   | Diffraction / calculation | `Priority 1`, `2`, or `3 — ...` | `PVT_kwHOACLKMM4BULjA` (P3), `PVT_kwHOACLKMM4BULi_` (P2), `PVT_kwHOACLKMM4BULi8` (P1) |
   | Documentation (Sphinx, content) | `Documentation` | `PVT_kwHOACLKMM4BUWg2` |

   Status field IDs are the same across all project boards:
   - Field ID: `PVTSSF_lAHOACLKMM4BULjAzhBVdT0` (P3), `PVTSSF_lAHOACLKMM4BUWg2zhBfRCE` (Documentation)
   - In Progress option ID: `47fc9ee4` (same on every board)

   ```bash
   # Set milestone when creating the PR
   gh pr create --milestone "Priority N — ..." ...   # or --milestone "Documentation"

   # Add the PR to the project board immediately after creation
   PR_NODE=$(gh api repos/OWNER/REPO/pulls/N --jq '.node_id')
   ITEM_ID=$(gh api graphql -f query="mutation {
     addProjectV2ItemById(input: {projectId: \"PROJECT_ID\", contentId: \"$PR_NODE\"})
     { item { id } } }" --jq '.data.addProjectV2ItemById.item.id')

   # Set the PR's project board status to "In Progress"
   # (STATUS_FIELD_ID and IN_PROGRESS_OPTION_ID are the same values used
   # for the issue status update in step 1 above)
   gh api graphql -f query="mutation {
     updateProjectV2ItemFieldValue(input: {
       projectId: \"PROJECT_ID\"
       itemId: \"$ITEM_ID\"
       fieldId: \"STATUS_FIELD_ID\"
       value: { singleSelectOptionId: \"IN_PROGRESS_OPTION_ID\" }
     }) { projectV2Item { id } } }"
   ```

3. **Never close issues manually** — let the merge automation do it.

   **When creating a new project board**, set the default view layout to
   **Board** (Kanban) so cards are visible by status column.  After creation,
   navigate to the project, open the view settings, and select *Board* as the
   layout type.  This cannot currently be set via the `gh` CLI; it must be
   done in the GitHub web UI immediately after the project is created.

4. **All code changes must be made on a feature branch**, never directly on
   `main` — the only exception is truly trivial changes (e.g. a one-word
   typo fix in a doc comment).  Branch names must start with the issue
   number followed by a short, hyphen-separated description of the topic:
   ```
   <issue-number>-<concise-topic>
   ```
   Examples: `2-motor-limits`, `1-wavelength-energy`, `9-diffraction-modes`

5. **On the feature branch**, if the completed issue corresponds to a section
   in `docs/roadmap.md`, check its box from `[ ]` to `[x]` as part of the
   same branch and include it in the PR.  The roadmap update is accepted when
   the PR is merged — no separate post-merge commit to `main` is needed.

Contributed by: OpenCode (argo/claudesonnet46)

---

## Project overview

`ad_hoc_diffractometer` is a Python package for describing multi-circle
diffractometer geometries used in X-ray and neutron crystallography.  It
provides:

- A class-based description of diffractometer stages (rotary axes) and
  their stacking order
- Predefined factory functions for standard synchrotron and laboratory
  diffractometer geometries
- Crystallographic lattice calculations (B matrix, reciprocal lattice)
- Caller-facing axis notation (+x, -z, vertical, lateral, longitudinal)
- A geometry registry (`list_geometries()`) for all predefined geometries
- Display precision control at package and instance level

The package is intended to grow toward a full diffraction calculation
engine (UB matrix, angle calculations, operating modes).  See
`docs/source/roadmap.md` for the planned feature list.

---

## Repository layout

```
diffractometer/                  # project root (git repo)
├── AGENTS.md                    # this file
├── pyproject.toml               # package metadata, ruff, isort, pytest config
├── docs/
│   ├── Makefile / make.bat      # sphinx build helpers
│   └── source/                  # Sphinx documentation source
│       ├── conf.py              # Sphinx configuration
│       ├── index.rst            # root toctree
│       ├── roadmap.md           # planned features — read before adding features
│       ├── install.md           # installation instructions
│       ├── api.rst              # AutoAPI stub
│       ├── changes.md           # includes CHANGES.md
│       ├── direct-lattice.md    # background material
│       ├── problem1.md / problem2.md / problem3.md
│       ├── problem1_solution.tex / .pdf
│       ├── fourcv_alignment_howto.ipynb
│       └── _static/
│           └── switcher.json    # version switcher stub
├── references/                  # journal articles and reference documents
│   ├── 1967 Busing and Levy a05492.pdf   # foundational four-circle paper
│   ├── 1999-JAppl-Cryst-32-614-623-H-You-psic-4S+2D/  # You (1999) psic
│   ├── 1993 J Appl Cryst 26 706 Lohmeier and Vlieg sixc/
│   ├── 2016 RefModuleMatSciMatEng_01215.pdf  # Walko — geometry survey
│   └── 2020-12-13-fourcc-alignment-7-id-c/  # real alignment session data
├── src/
│   └── ad_hoc_diffractometer/   # package source
│       ├── __init__.py          # public API
│       ├── constants.py         # XHAT, YHAT, ZHAT
│       ├── axes.py              # parse_axis(), axis_label(), kappa_axis()
│       ├── rotation.py          # rotation_matrix() — Rodrigues formula
│       ├── stage.py             # Stage class
│       ├── geometry.py          # AdHocDiffractometer class
│       ├── factories.py         # all geometry factories + registry
│       ├── lattice.py           # Lattice class, b_matrix(), standalone fns
│       └── display.py           # get/set_precision(), fmt()
└── tests/
    ├── test_<module>.py         # one file per source module (see testing instructions)
    ├── test_regression_issue_N.py  # cross-module regression tests (named by issue)
    └── helpers.py / conftest.py # shared test infrastructure
```

---

## Setup and install

```bash
pip install -e .                   # install in editable mode
python3 -m pre_commit install      # install pre-commit hooks (required once)
```

---

## Running tests

```bash
python3 -m pytest                        # run all tests
python3 -m pytest tests/test_lattice.py  # lattice tests only
python3 -m pytest -q                     # quiet summary
```

All tests must pass before committing.  Pre-commit runs automatically on
`git commit`.  To run it manually:

```bash
python3 -m pre_commit run --all-files
```

---

## Code style guidelines

- **Import alias**: `import ad_hoc_diffractometer as ahd`
- **Python**: ≥ 3.10 (uses `X | Y` union types)
- **Dependencies**: `numpy` only; no scipy, no sympy
- **Style**: ruff (E, W, F, UP, B rules) + ruff-format, line length 88;
  `E501` ignored (slight overruns tolerated)
- **Import sorting**: isort with `force_single_line = true`
- **No `geometry_` prefix** on factory functions
- **`display.fmt(value, digits)`** for all floating-point display;
  never use f-strings with hardcoded precision for user-facing output

---

## Testing instructions

**Test file naming**

- Each test file `tests/test_<module>.py` corresponds to exactly one source
  module `src/ad_hoc_diffractometer/<module>.py`.  Keep this one-to-one
  mapping when adding new modules.
- Do **not** create `test_<feature>.py` files that cover content from an
  existing module — add the tests to the matching `test_<module>.py`.
- Exception: purpose-specific files such as `test_regression_issue_N.py`
  are permitted for cross-module bug reports added later in the project's
  life, where the fix spans multiple modules and there is no single natural
  home.  These files should say so in their module docstring.

**Test content**

- Every parametrized test set uses `pytest.param(..., id="descriptive-id")`
- Every parametrized set includes a `context` parameter that is either
  `does_not_raise()` (from `contextlib.nullcontext`) or
  `pytest.raises(ExcType, match=re.escape("expected text"))`
- The `match` string is embedded directly in `pytest.raises(...)` —
  it is **not** passed as a separate parameter to the test function
- All test code runs inside `with context:` block
- Add tests for every new function, class, and error case
- Run `python3 -m pytest -q` after any change and fix all failures before
  committing

---

## Coordinate convention

The default coordinate system follows You (1999):

| Basis vector | Physical direction |
|---|---|
| `XHAT` (+x) | vertical (out of floor) |
| `YHAT` (+y) | longitudinal (along beam, toward equipment) |
| `ZHAT` (+z) | lateral (to our left facing the equipment) |

The Busing & Levy (1967) convention (used in `fourcv`, `fourch`,
`kappa4cv`, `kappa4ch`) has x=lateral, y=longitudinal, z=vertical.

Basis dicts (`_BASIS_YOU`, `_BASIS_BL`) in `factories.py` encode these.

---

## Axis sign convention

The sign of a stage's axis vector encodes its handedness:

- `+nHat` → right-handed rotation about nHat (thumb along +nHat, fingers curl positive)
- `-nHat` → left-handed rotation about nHat (equivalent to right-handed with negated angle)

This is the `parse_axis()` / `kappa_axis()` / `axis_label()` system in
`axes.py`.  Stage internal attributes always store a numpy array; the
string notation (`"+x"`, `"-z"`, `"vertical"`, `"-lateral"`) is the
caller-facing interface only.

---

## Geometry factories

All factory functions are decorated with `@register_geometry` and appear
in `list_geometries()`.  Naming convention:

| Suffix | Meaning |
|---|---|
| `v` | vertical scattering plane (synchrotron) — ttheta rotates about the lateral axis |
| `h` | horizontal scattering plane (laboratory) — ttheta rotates about the vertical axis |
| no suffix | detector convention unambiguous (psic, sixc) or compound (zaxis, s2d2) |

Current factories: `psic`, `fourcv`, `fourch`, `sixc`, `kappa4cv`,
`kappa4ch`, `kappa6c`, `zaxis`, `s2d2`, `fivec`.

Walko (2016) designation system is noted in docstrings (S3D1, S4D2, etc.).

The kappa axis is:
```python
kappa_axis = vertical * cos(alpha) + lateral * sin(alpha)
```
where `alpha` is measured from the vertical axis toward the lateral axis
(confirmed by ITC Vol. C Sec. 2.2.6 and Walko 2016).  Default alpha = 50°.

---

## Lattice class

`Lattice` deduces the crystal system from the minimum number of supplied
parameters.  The 7 systems and their minimum inputs:

| System | Minimum |
|---|---|
| cubic | a |
| tetragonal | a, c |
| hexagonal | a, c, gamma=120 (gamma must be explicit) |
| trigonal | a, alpha (alpha ≠ 90) |
| orthorhombic | a, b, c |
| monoclinic | a, b, c, beta |
| triclinic | a, b, c, alpha, beta, gamma |

Properties `cartesian_lattice_vectors`, `reciprocal_lattice_vectors`, and
`B` are lazy-computed and cached.  Any parameter setter invalidates the
cache and re-deduces the crystal system.

The B matrix follows the I16 convention: `(b1, b2, b3) = 2π * B.T`.

---

## Key references

| Citation | What it defines |
|---|---|
| Busing & Levy, Acta Cryst. 22, 457-464 (1967) | Four-circle geometry, B matrix, U matrix, UB matrix |
| Bloch, J. Appl. Cryst. 18, 33-36 (1985) | Z-axis geometry |
| Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987) | Five-circle geometry |
| Lohmeier & Vlieg, J. Appl. Cryst. 26, 706-716 (1993) | sixc surface geometry |
| Evans-Lutterodt & Tang, J. Appl. Cryst. 28, 318-326 (1995) | S2D2 geometry |
| You, J. Appl. Cryst. 32, 614-623 (1999) DOI:10.1107/S0021889899001223 | psic 4S+2D six-circle; axis sign conventions (mixed handedness) |
| ITC Vol. C, Sec. 2.2.6 (2006) DOI:10.1107/97809553602060000577 | Kappa 50° convention; normal-beam equatorial geometry |
| Walko, Ref. Module Mater. Sci. Mater. Eng. (2016) | Geometry survey; S/D designation system; kappa, zaxis, s2d2 |

---

## What is NOT yet implemented (see docs/source/roadmap.md)

- Wavelength / energy (Priority 1.1)
- Motor limits per stage (Priority 1.2)
- Kappa alpha queryable on instance (Priority 1.3)
- Motor-angle-to-phi-frame conversion (Priority 2.0) — prerequisite for:
  - U matrix computation (Priority 2.1)
  - UB matrix computation (Priority 2.2)
  - Orienting reflections (Priority 2.3)
- Diffraction modes as part of geometry constructor (Priority 3.1)
- Azimuthal reference vector (Priority 3.2)
- Detector geometry parameters (Priority 3.3)
- Alternative calculation engines: Q-space, d-spacing (Priority 3.4)
- Surface geometry: incidence/emergence angles, ψ scans (Priority 3.5, 3.6)
- Diffractometer inclination relative to beam (Priority 3.7)
