# AGENTS.md — ad_hoc_diffractometer

This file provides context for AI coding agents working on this project.

---

## Operating directives (pre-v1.0)

These rules apply for every issue and every PR until the project tags
v1.0.  They override anything in the rest of this file or in the
agent's general training that conflicts.

### No deprecation cycles

Rename freely.  Do not add forwarding aliases, `DeprecationWarning`
emitters, alias maps, `REMOVE-AT-V1.0` markers, or any other
deprecation-cycle infrastructure for any rename, signature change, or
removal.  The whole package is still pre-1.0; users have no API
stability contract yet.  Make the change cleanly, update every
in-tree call site, document the change in `CHANGES.md`, move on.

### Scope is the issue

The issue text defines the scope.  If something falls under the
issue's title or stated motivation, it is in scope by default.  Do
not unilaterally narrow the scope or split work into "follow-up"
issues.  If the scope is genuinely ambiguous, ask exactly one
clarifying question, then proceed.

### No volunteered follow-up issues

Do not file new issues, propose new issues, or build "future v1.0
cleanup" tracking unless the user asks.  Pre-1.0 cleanup is just
"edit the code now."

### Response length matches the request

When the user asks a brief question, answer briefly: a sentence or a
short list, no tables, no preamble, no recap.  When the user asks
for analysis or a plan, longer is fine.  Match the form of the
request.

### When in doubt, ask one short question

A single clarifying question is cheaper than building the wrong
thing and unwinding it.  Do not stack multiple questions; ask one,
wait, then proceed.

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

   Status field IDs and the corresponding option IDs are project-
   specific.  Discover them per board with the GraphQL query
   `node(id: "<PROJECT_ID>") { ... fields(first: 20) { ... } }` and
   reuse them where AGENTS.md examples below show `STATUS_FIELD_ID`,
   `IN_PROGRESS_OPTION_ID`, etc.  In Progress is present on every
   board; "In review" is present only on boards that have been
   configured with that status — leave the status at "In Progress"
   on boards where Todo / In Progress / Done are the only options.

   **Match the issue's boards.**  Before creating the PR, query the
   issue's `projectItems` to discover every board it belongs to.  Add the
   PR to **all** of those boards — not just the one implied by the domain
   table above.

   **Set "In review" when code is complete.**  Once the PR is pushed, set its
   project board status to **"In review"** (on boards that have that status;
   leave "In Progress" where only Todo / In Progress / Done exist).

   ```bash
   # Set milestone when creating the PR
   gh pr create --milestone "Example Milestone" ...

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

4. **Create a feature branch before making any file changes.**  This is the
   very first step after reading the issue — create the branch **before**
   editing, creating, or deleting any file.  All changes must be made on a
   feature branch, never directly on `main` — the only exception is truly
   trivial changes (e.g. a one-word typo fix in a doc comment).  If you
   have already modified files on `main` by mistake, stash the changes,
   create the branch, then pop the stash.  Branch names must start with the
   issue number followed by a short, hyphen-separated description of the
   topic:
   ```
   git checkout -b <issue-number>-<concise-topic>
   ```
   Examples: `2-motor-limits`, `1-wavelength-energy`, `9-diffraction-modes`

5. removed

6. **Monitor CI** after the PR is opened.  Watch for test failures, lint
   errors, or coverage gaps reported by the CI checks.  Push additional
   commits to the feature branch to resolve any failures — do **not** force-
   push unless explicitly requested.

7. **Resolve PR comments**.  Address reviewer feedback with further commits
   on the same branch.  Update the PR description if the scope changes
   materially.

8. **Merge** once all CI checks are green and the PR is approved (or
   self-approved if the repository allows it).  Use the merge strategy
   preferred by the repository (typically squash-merge or merge commit —
   follow the existing pattern in the commit history).  The `closes #N`
   keyword in the PR body closes the issue and moves the project card to
   "Done" automatically.

9. **Local cleanup** after a successful merge:
   ```bash
   git checkout main
   git pull --ff-only origin main
   git branch -d <feature-branch>
   ```
   Confirm the issue is closed on GitHub and the project card has moved to
   "Done" before starting the next item.

**Confirm before acting on implied direction.**  If a user's message could
be interpreted as either a question or an instruction, ask for clarification
before taking action (committing, pushing, opening PRs, etc.).

Contributed by: OpenCode (argo/claudesonnet46)

---

## Project overview

`ad_hoc_diffractometer` is a Python package for describing multi-circle
diffractometer geometries used in X-ray and neutron crystallography.  It
provides:

- A class-based description of diffractometer stages (rotary axes) and
  their stacking order
- Demo geometries for standard synchrotron and laboratory
  diffractometers
- Crystallographic lattice calculations (B matrix, reciprocal lattice)
- Caller-facing axis notation (+x, -z, vertical, transverse, longitudinal)
- A geometry registry (`list_geometries()`) for all predefined geometries
- Display precision control at package and instance level

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
│       ├── install.md           # installation instructions
│       ├── api.rst              # AutoAPI stub
│       ├── changes.md           # includes CHANGES.md
│       ├── direct-lattice.md    # background material
│       ├── problem1.md / problem2.md / problem3.md
│       ├── problem1_solution.tex / .pdf
│       ├── fourcv_alignment_howto.ipynb
│       └── _static/
│           └── switcher.json    # version switcher stub
├── references/                  # journal articles and reference documents (excluded from repo)
│   ├── 1967 Busing and Levy a05492.pdf   # foundational four-circle paper
│   ├── 1999-JAppl-Cryst-32-614-623-H-You-psic-4S+2D/  # You (1999) psic
│   ├── 1993 J Appl Cryst 26 706 Lohmeier and Vlieg sixc/
│   ├── 2016 RefModuleMatSciMatEng_01215.pdf  # Walko — geometry survey
│   └── 2020-12-13-fourcc-alignment-7-id-c/  # real alignment session data
├── src/
│   └── ad_hoc_diffractometer/   # package source
│       ├── __init__.py          # public API (Tier 1 names only)
│       ├── geometries/          # declarative YAML demo geometries (#267)
│       │   ├── __init__.py      # marker for importlib.resources
│       │   ├── schema.json      # JSON Schema for the declarative format
│       │   └── *.yml            # one file per declarative demo geometry
│       ├── geometry_loader.py   # YAML → AdHocDiffractometer loader (#267)
│       ├── factories.py         # geometry registry + shared definitions
│       ├── constants.py         # XHAT, YHAT, ZHAT
│       ├── axes.py              # parse_axis(), axis_label(), kappa_axis()
│       ├── rotation.py          # rotation_matrix() — Rodrigues formula
│       ├── stage.py             # Stage class
│       ├── diffractometer.py    # AdHocDiffractometer class
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

## Copyright handling

The repo uses three coordinated mechanisms to keep copyright text
consistent.  All three are wired into pre-commit and run automatically.

### 1. Per-file headers — `.copyright.txt` + `insert-license`

`.copyright.txt` is the single source of truth for the per-file header
text:

```text
Copyright (c) 2026-2026 UChicago Argonne, LLC
SPDX-License-Identifier: LicenseRef-UChicago-Argonne-LLC-License
```

The `Lucas-C/insert-license` pre-commit hook propagates that header to
every covered Python file on each commit.  Three hook instances cover:

* `src/ad_hoc_diffractometer/*.py` (excluding `_version.py`)
* `tests/test_*.py`
* `scripts/*.py`

To change the per-file header text project-wide, edit only
`.copyright.txt` and run `pre-commit run --all-files`.

### 2. Year-range bump — `scripts/update_copyright_year.py`

A **local** pre-commit hook (`update-copyright-year`) runs the
`scripts/update_copyright_year.py` script on every commit.  The script
rewrites the pattern `<START>-<OLD_END>` → `<START>-<CURRENT_YEAR>` in
each file listed in its `TARGET_FILES`:

* `.copyright.txt`
* `LICENSE` (line 1 only — the licence body is verbatim per ANL legal
  and is never edited)
* `docs/source/conf.py`

The script exits non-zero when it changes anything, so pre-commit fails
and the developer stages the rewrite before retrying the commit.  On
January 1 of any year, the next commit on any branch will bring every
year span forward.

### 3. Sphinx docs — static year range in `conf.py`

`docs/source/conf.py` declares `copyright` as a static year-range string,
not a build-time-dynamic `datetime.now().year` expression.  This is
intentional: the bump script is the single mechanism that keeps
`LICENSE`, `.copyright.txt`, and the rendered docs aligned, and a
dynamic expression would hide drift between them.

### What NOT to edit by hand

* The `Copyright (c) YYYY-YYYY ...` line on individual `.py` files —
  edit `.copyright.txt` and re-run pre-commit instead.
* The end year in `LICENSE` line 1, `docs/source/conf.py`, or
  `.copyright.txt` on January 1 — the bump script handles it.
* The body of `LICENSE` — verbatim per ANL legal.

### What to edit by hand

* The **start** year in any year range (project inception year);
  the script never touches it.
* Adding new files to `TARGET_FILES` if a new file legitimately
  contains a year range.

---

## Benchmark tests

The benchmark test suite (`tests/test_benchmark.py`) includes slow tests
that run `forward()`/`inverse()` across **all** registered geometries and
modes.  These are marked `@pytest.mark.slow_benchmark` and are **excluded
by default** from `pytest` via `addopts` in `pyproject.toml`.

Fast replacement tests (using a single monkeypatched geometry) maintain
100 % code coverage on `benchmark.py` during normal development.

**When to run full benchmarks locally.**  If your changes touch any module
in the hot path or geometry construction layer, run the full benchmark
suite before committing:

```bash
python3 -m pytest -m slow_benchmark --no-cov -q
```

Hot path: `forward.py`, `kappa.py`, `mode.py`, `orientation.py`,
`rotation.py`, `reference.py`, `surface.py`

Geometry construction: `geometry_loader.py`, `geometries/*.yml`,
`factories.py`, `diffractometer.py`, `stage.py`, `axes.py`,
`constants.py`

CI runs these automatically (via `.github/workflows/benchmark.yml`)
when any of these files change on `main` or in a pull request.

---

## Documentation build verification

A clean Sphinx build is **necessary but not sufficient** to confirm
that documentation changes render correctly.  Sphinx will emit "build
succeeded" even when reStructuredText markup leaks through into the
HTML as literal text.  Common failure modes:

- A role written inside a Napoleon field-name slot
  (e.g. ``:attr:`~module.Class.name` : float`` in a NumPy
  ``Attributes`` section) is pasted verbatim into the
  ``.. py:attribute::`` directive, so the rendered page shows the raw
  ``:attr:`~module.Class.name``` text instead of a link.
- Inline roles that should produce links (``:func:``, ``:class:``,
  ``:meth:``, ``:attr:``, ``:doc:``, ``:ref:``) silently render as
  literal text when their target cannot be resolved or when they
  appear in a context that does not parse RST roles (directive option
  values, certain table cells, etc.).
- MyST/Markdown files containing accidentally-mixed RST syntax
  produce no warning but display the literal markup.

**Required workflow for any change that touches a docstring, a
``docs/source/**`` file, or an AutoAPI Jinja2 template:**

1. Run ``make -C docs clean html`` and confirm zero warnings.
2. Open the affected rendered page(s) under ``docs/build/html/`` (or
   ``grep`` the HTML) and visually verify that every cross-reference
   appears as a real ``<a>`` link, **not** as literal ``:role:`...``` text.
3. For any class with a NumPy ``Attributes`` docstring section, also
   ``grep`` the AutoAPI page for ``&#58;attr&#58;`` /
   ``&lt;span class="pre"&gt;:`` patterns that would indicate role
   markup leaking into the HTML.

**Quick check** (from the repository root after rebuilding the docs).
Run **both** patterns — they catch different leak shapes:

```bash
# Pattern 1: literal :role:`...` (RST role text passed straight through).
# Catches leaks in .rst files and in MyST contexts where backticks were
# preserved verbatim.  Excludes docs/build/html/_modules/ (sphinx.ext.viewcode
# pages display source code verbatim, so role markup there is expected).
grep -rE ':(attr|func|class|meth|mod|doc|ref|obj|data|exc|term):`' \
     docs/build/html/ --include='*.html' \
     --exclude-dir=_modules \
     --exclude-dir=_static \
  && echo "FAIL: literal RST role leak — see above." \
  || echo "OK: no literal RST role leaks."

# Pattern 2: bare :role: prefix that survived even when MyST rewrote the
# backticked content into <code>.  This is the failure mode when an RST
# role (e.g. :attr:`~mod.Class.member`) appears in a MyST/Markdown file —
# the MyST parser strips the backticks but leaves the role prefix as
# literal text in front of an inline-code element.  Filter out legitimate
# href/title attributes and proper xref/reference link bodies; what
# remains is leaked role-prefix text.
grep -rhE ':(attr|func|class|meth|mod|doc|ref|obj|data|exc|term):' \
     docs/build/html/ --include='*.html' \
     --exclude-dir=_modules \
     --exclude-dir=_static \
  | grep -vE 'href=|title=|<title>|class="reference|class="xref' \
  | grep -E ':(attr|func|class|meth|mod|doc|ref|obj|data|exc|term):' \
  | head -10
# (above prints any leak lines; empty output = OK.)
```

This check belongs alongside the test suite in any pre-PR checklist
that involves documentation work.

In MyST/Markdown files use the brace form ``{role}\`target\``` (for
example ``{attr}`~mod.Class.member```) — never the colon form, which
is silently passed through as literal text.

---

## Generated artefacts

Several documentation files are **generated from source code or other
source files** and committed to git.  When the underlying source changes,
these artefacts must be rebuilt and included in the same branch.

### Geometry diagrams (SVG + HTML)

Interactive Plotly HTML figures and static SVG fallbacks for all 10
demo geometries live in
`docs/source/_static/geometries/<geometry>/<geometry>.{html,svg}`.

They are generated by `tools/generate_geometry_drawings.py`, which
instantiates a
`GeometryAxisFigure` (from `drawing.py`) for each registered geometry
and calls `write_html()` and `write_image()`.

**Rebuild command** (from repository root):

```bash
python tools/generate_geometry_drawings.py
```

**When to rebuild**: any change to `drawing.py`, `geometries/*.yml`,
`geometry_loader.py`, `factories.py`, `stage.py`, or `axes.py` that
affects stage names, axis labels, physical direction names, basis
dicts, or the drawing layout.

**Dependencies**: `plotly` and `kaleido` (for SVG export).  These are
not declared in `pyproject.toml` runtime deps — install them manually:

```bash
pip install plotly kaleido
```

### LaTeX / PDF (problem1_solution)

`docs/source/problem1_solution.tex` is a standalone LaTeX document
derived from the content of `docs/source/problem1.md`.  The PDF is
committed alongside it.

**Rebuild command** (from `docs/source/`):

```bash
pdflatex -interaction=nonstopmode problem1_solution.tex
pdflatex -interaction=nonstopmode problem1_solution.tex   # second pass
rm -f problem1_solution.{aux,log,out}                     # clean up
```

**When to rebuild**: any change to `problem1_solution.tex`.  If
`problem1.md` changes in a way that affects the worked solution, update
the `.tex` file to match, then rebuild the PDF.

### Jupyter notebook (fourcv_alignment_howto)

`docs/source/howto/fourcv_alignment_howto.ipynb` is an authored notebook
(not auto-generated).  It contains markdown cells describing the `fourcv`
geometry.  Edit it directly when terminology, axis labels, or geometry
descriptions change.  Sphinx executes it during `make html` via
`myst-nb`, so its code cells must remain runnable.

---

## Implementation philosophy

`ad_hoc_diffractometer` is a **pure-Python** package.  The runtime
dependencies beyond the Python Standard Library are **NumPy** (numerical
arrays and linear algebra) and **PyYAML** (parser for the declarative
geometry files in ``src/ad_hoc_diffractometer/geometries/``; see issue
#267).  There is no dependency on scipy, sympy, or any other scientific
library.  Keep it that way: every new algorithm must be implementable
with NumPy alone.

---

## Code style guidelines

- **Import alias**: `import ad_hoc_diffractometer as ahd`
- **Demo geometries**: `g = ahd.make_geometry("fourcv")` — the demo
  geometries are declarative YAML files under
  `ad_hoc_diffractometer/geometries/`, registered automatically by the
  loader at import time and discoverable via
  `ad_hoc_diffractometer.list_geometries()`
- **Tier 1 (top-level)**: Only ~23 names are exported from `__init__.py`
  (core classes, orientation, modes, registry).  All other names are
  accessed via their submodule: `ahd.display.fmt()`,
  `ahd.radiation.energy_to_wavelength()`, `ahd.conversions.hkl_to_d()`, etc.
- **Python**: ≥ 3.10 (uses `X | Y` union types)
- **Dependencies**: `numpy` and `pyyaml` only; no scipy, no sympy
- **Style**: ruff (E, W, F, UP, B rules) + ruff-format, line length 88;
  `E501` ignored (slight overruns tolerated)
- **Import sorting**: isort with `force_single_line = true`
- **No `geometry_` prefix** on demo-geometry function names
- **`display.fmt(value, digits)`** for all floating-point display;
  never use f-strings with hardcoded precision for user-facing output
- **US English spellings** — see the dedicated section
  [US English spelling (mandatory)](#us-english-spelling-mandatory) below.

---

## US English spelling (mandatory)

**Every word you write — in code, comments, docstrings, error
messages, log messages, exception text, type names, identifiers,
filenames, branch names, commit messages, PR descriptions, issue
comments, YAML strings, JSON descriptions, and Markdown / reST
documentation — must use US (American) English spelling.**  This rule
is non-negotiable and is checked at review time.

This applies even when paraphrasing or summarizing material that uses
British spellings (e.g. literature references): retype the prose into
US English when incorporating it into the project.  Verbatim quotations
inside quotation marks may preserve the source's spelling; everything
else must conform.

### Mandatory substitutions

| Use (US) | NOT (British) |
|---|---|
| `analyze`, `analyzer`, `analyzed`, `analyzing` | `analyse`, `analyser`, `analysed`, `analysing` |
| `behavior`, `behaviors` | `behaviour`, `behaviours` |
| `center`, `centers`, `centered`, `centering` | `centre`, `centres`, `centred`, `centring` |
| `color`, `colors`, `colored`, `coloring` | `colour`, `colours`, `coloured`, `colouring` |
| `characterize`, `characterized`, `characterizing` | `characterise`, `characterised`, `characterising` |
| `customize`, `customized` | `customise`, `customised` |
| `defense`, `defensive` | `defence`, `defencive` |
| `emphasize`, `emphasized` | `emphasise`, `emphasised` |
| `factorize`, `factorized` | `factorise`, `factorised` |
| `favor`, `favorite` | `favour`, `favourite` |
| `fiber`, `fibers` | `fibre`, `fibres` |
| `fulfill`, `fulfilled` | `fulfil`, `fulfilled` (British single-l) |
| `generalize`, `generalized`, `generalizing` | `generalise`, `generalised`, `generalising` |
| `gray` | `grey` |
| `honor`, `honors`, `honored`, `honoring` | `honour`, `honours`, `honoured`, `honouring` |
| `initialize`, `initialized`, `initializing` | `initialise`, `initialised`, `initialising` |
| `labeled`, `labeling` | `labelled`, `labelling` |
| `liter`, `liters` | `litre`, `litres` |
| `maximize`, `maximized` | `maximise`, `maximised` |
| `meter`, `meters` (the unit) | `metre`, `metres` |
| `millimeter`, `millimeters` | `millimetre`, `millimetres` |
| `minimize`, `minimized`, `minimizing` | `minimise`, `minimised`, `minimising` |
| `modeled`, `modeling` | `modelled`, `modelling` |
| `normalize`, `normalized`, `normalizing` | `normalise`, `normalised`, `normalising` |
| `optimize`, `optimized`, `optimizing` | `optimise`, `optimised`, `optimising` |
| `organization` | `organisation` |
| `parameterize`, `parameterized` | `parameterise`, `parameterised` |
| `polarization`, `polarized`, `polarizer` | `polarisation`, `polarised`, `polariser` |
| `practice` (verb) | `practise` |
| `program` | `programme` (in software contexts) |
| `realize`, `realized` | `realise`, `realised` |
| `recognize`, `recognized`, `recognizing` | `recognise`, `recognised`, `recognising` |
| `serialize`, `serialized`, `serializing` | `serialise`, `serialised`, `serialising` |
| `signaled`, `signaling` | `signalled`, `signalling` |
| `specialized`, `specializing` | `specialised`, `specialising` |
| `stabilize`, `stabilized` | `stabilise`, `stabilised` |
| `summarize`, `summarized` | `summarise`, `summarised` |
| `synthesize`, `synthesized` | `synthesise`, `synthesised` |
| `theater` | `theatre` |
| `traveled`, `traveling`, `traveler` | `travelled`, `travelling`, `traveller` |
| `utilize`, `utilized` | `utilise`, `utilised` |
| `visualize`, `visualized`, `visualizing` | `visualise`, `visualised`, `visualising` |

The list above is illustrative, not exhaustive.  Whenever you are
unsure, prefer the spelling used by the *Merriam-Webster* dictionary
or the OED's "US" tag.  The rule applies to **every variant** of the
listed roots (`-ize / -ization / -izing / -ized` etc.).

### Enforcement check

Before opening a pull request — and again after addressing reviewer
comments — run this check on every file you have added or modified
(replace `$FILES` with `git diff --name-only main...HEAD` or an
equivalent list, **excluding `AGENTS.md` itself** because the
substitution table above intentionally contains British spellings as
negative examples):

```bash
grep -EnH "\b(honour|honours|honoured|honouring|colour|colours|coloured|colouring|behaviour|behaviours|analyse|analysed|analyser|analysers|analyses|analysing|polariser|polarisers|polarisation|polarised|polarising|normalise|normalised|normalises|normalising|minimise|minimised|minimises|minimising|maximise|maximised|maximises|maximising|optimise|optimised|optimises|optimising|generalise|generalised|generalises|generalising|recognise|recognised|recognises|recognising|characteris|millimetre|centre|centres|centred|centring|customise|customised|emphasise|emphasised|factorise|factorised|fibre|fibres|favour|favourite|grey|initialise|initialised|labelled|labelling|litre|litres|metre|metres|modelled|modelling|organisation|parameterise|parameterised|realise|realised|serialise|serialised|signalled|signalling|specialised|specialising|stabilise|stabilised|summarise|summarised|synthesise|synthesised|theatre|travelled|travelling|traveller|utilise|utilised|visualise|visualised|visualising|programme)\b" $FILES
```

The check must produce **no output** (other than expected hits inside
the AGENTS.md substitution table, which should be excluded from
`$FILES`).  Any hit elsewhere must be corrected before review.  This
check belongs alongside `pytest`, `pre-commit`, and the documentation
grep checks in your pre-PR checklist.

### Pre-existing British spellings

A small number of British spellings remain in pre-existing files that
were authored before this rule was tightened.  Do not silently fix
them in unrelated PRs — open a dedicated cleanup PR (or include them
explicitly in a PR whose scope already touches those files) so the
diff is reviewable.

---

## CHANGES.md style

`CHANGES.md` is a user-facing changelog, not a design document.  Each
bullet is **a few words plus the issue/PR reference**.  The issue/PR
reference carries the full load of explanation — readers who want
details follow the link.  The bullet itself names the change in the
fewest words that identify it.

**Hard length cap: ≤ 12 words before the ` (#NNN)` reference.**  If
you find yourself writing a second clause, a semicolon, an "and", a
parenthetical aside, or a "users must …" instruction in the bullet
itself, the bullet is too long — cut it down.  All of that material
belongs in the referenced issue or PR, not in the changelog.

**Use v0.6.0 and earlier as the style reference.**  Several recent
entries (v0.7.0, v0.8.0) drifted toward multi-sentence bullets with
implementation detail and benchmark numbers; do not imitate them.
v0.6.0 is the canonical pattern: each bullet is a short headline
followed by ` (#NNN)` or ` (#NNN, #MMM)`.

Examples of acceptable length:

```
- Rotation composition and `ub_identity` corrected. (#280)
- B matrix for non-cubic lattices now matches BL1967. (#280)
- Doc version-switcher dropdown now shows every published version. (#273)
- Switch to declarative YAML for demo geometries. (#267)
```

Examples of bullets that are **too long** and must be cut:

```
- Compose stage rotations outermost-leftmost (BL1967); place crystal a*
  along the beam at ub_identity; saved UB matrices must be re-derived. (#280)
- Default UB matrix (ub_identity) now places the crystal a-axis along
  the beam direction; sessions relying on the previous default must
  explicitly set U or measure orienting reflections. (#NNN)
```

Both of the "too long" examples pack the issue body into the
changelog.  The issue/PR already has the body; the bullet should not
duplicate it.

- One bullet per change.  A few words followed by the issue/PR
  reference.  No second sentence, no rationale, no cross-validation
  numbers, no comparison to other packages, no internal implementation
  detail, no benchmark numbers, no migration instructions.
- Use the existing subsection headings: `Breaking changes`, `Behavior
  change`, `Added`, `Changed`, `Fixed`.  Omit any that are empty.
- **Sort bullets alphabetically within each `###` subsection.**
  Case-insensitive; treat surrounding backticks, asterisks, and
  punctuation as transparent so the leading *word* is what sorts.
  This makes diff review predictable and keeps related bullets near
  each other across releases.  Applies going forward (post-v0.11.1);
  do not silently re-sort historical release sections in unrelated
  PRs — open a dedicated cleanup PR if you want to backfill.
- A correctness fix that changes computed numerical output gets one
  short `Behavior change` bullet naming what changed.  Migration
  instructions (e.g. "re-derive saved UB matrices") belong in the
  referenced issue/PR, not in the bullet.
- Do **not** add bullets for test-suite changes (new tests,
  parametrization edits, test-only HKL substitutions), docstring
  fixes, or internal refactoring.  A `Fixed` entry implicitly
  includes its regression test.
- The release header is `## Release vX.Y.Z` followed by
  `Released YYYY-MM-DD.` and nothing else.  No prose paragraph.
- A single-fix release should be three or four lines of bullets, not
  more (see v0.9.0 for the model: two bullets, one line each).

When in doubt, err on the side of fewer words.  If a reviewer can
look at the bullet and identify what the change is about (well
enough to decide whether to click the issue link), it is long
enough.

Contributed by: OpenCode (argo/claudeopus47)

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
| `ZHAT` (+z) | transverse (to our left facing the equipment) |

The Busing & Levy (1967) convention (used in `fourcv`, `fourch`,
`kappa4cv`, `kappa4ch`) has x=transverse, y=longitudinal, z=vertical.

Basis dicts (`BASIS_YOU`, `BASIS_BL`) in `factories.py` encode these.
A third constant, `BASIS_DEFAULT` (vertical=YHAT, longitudinal=ZHAT,
transverse=XHAT), is the neutral fallback used by the declarative
geometry loader when a YAML file omits the `basis:` key.  It is
deliberately distinct from `BASIS_YOU` and `BASIS_BL` so that the
package does not appear to espouse either literature convention as a
"project default."

---

## Axis sign convention

The sign of a stage's axis vector encodes its handedness:

- `+nHat` → right-handed rotation about nHat (thumb along +nHat, fingers curl positive)
- `-nHat` → left-handed rotation about nHat (equivalent to right-handed with negated angle)

This is the `parse_axis()` / `kappa_axis()` / `axis_label()` system in
`axes.py`.  Stage internal attributes always store a numpy array; the
string notation (`"+x"`, `"-z"`, `"vertical"`, `"-transverse"`) is the
caller-facing interface only.

---

## Demo geometries

All demo geometries are declarative YAML files in
`src/ad_hoc_diffractometer/geometries/`, parsed at import time by
`geometry_loader.py` and registered in the geometry registry; they
appear in `list_geometries()`.  Access them as
`ahd.make_geometry("fourcv")`, etc.

Naming convention:

| Suffix | Meaning |
|---|---|
| `v` | vertical scattering plane (synchrotron) — ttheta rotates about the transverse axis |
| `h` | horizontal scattering plane (laboratory) — ttheta rotates about the vertical axis |
| no suffix | detector convention unambiguous (psic, sixc) or compound (zaxis, s2d2) |

Current demo geometries: `psic`, `fourcv`, `fourch`, `sixc`, `kappa4cv`,
`kappa4ch`, `kappa6c`, `zaxis`, `s2d2`, `fivec`.

Walko (2016) designation system is noted in docstrings (S3D1, S4D2, etc.).

The kappa axis is:
```python
kappa_axis = vertical * cos(alpha) + transverse * sin(alpha)
```
where `alpha` is measured from the vertical axis toward the transverse axis
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
