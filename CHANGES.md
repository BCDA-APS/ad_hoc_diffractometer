# Change History

User-facing changes by release.  See the
[development roadmap](https://github.com/prjemian/ad_hoc_diffractometer/blob/main/roadmap.md)
for planned future work and
[GitHub Issues](https://github.com/prjemian/ad_hoc_diffractometer/issues)
for the full issue tracker.

## Release v0.5.0

Released 2026-04-20.

### Breaking changes

- `__init__.py` trimmed to ~23 Tier 1 exports; preset factories moved
  to `ahd.presets.*` (#185)
- `BisectingMode` / `FixedAngleMode` replaced by `ConstraintSet` with
  `BisectConstraint`, `SampleConstraint`, `DetectorConstraint`,
  `ReferenceConstraint` (#148)

### Added

- 57 SVG axis diagrams and graphviz stage coupling diagrams (#163)
- Constraint system: `ConstraintSet`, `BisectConstraint`,
  `SampleConstraint`, `DetectorConstraint`, `ReferenceConstraint`;
  `EwaldSphereViolation`, `ConstraintViolation` exceptions (#148)
- Diffraction modes for all 10 geometries
  (#149, #150, #151, #152, #154, #155, #156)
- Double-diffraction 4D simultaneous solver for 4 geometries (#184)
- `drawing.py`: `geometry_dot()`, `draw_stage_axis()`,
  `draw_geometry_axes()` (#163)
- HOWTO: diffractometer inclination (#164)
- Kappa virtual-angle dispatcher: `constant_omega`/`chi`/`phi` modes;
  `kappa.py` module (#174)
- Psi-constant validation filter for 6 geometries (#184)
- Qaz detector constraint: `lifting_detector_*` modes for psic,
  kappa6c (#177)
- Quick Start guide; Diátaxis User Guide restructure (#146)
- Surface diffraction forward solvers for sixc, zaxis, s2d2 (#175)
- `tools/generate_geometry_drawings.py` with `list_geometries()`
  auto-discovery (#163)

### Changed

- Constraint/mode docs overhauled across all geometry pages (#173)
- Doc build timestamp in Chicago local time (#165)
- Geometry diagrams use SVG instead of PNG (#163)
- `presets.py` split from `factories.py` (#185)

### Fixed

- `references.md` broken links (#152)
- Unused variable and formatting cleanups (#172, #178)

## Release v0.4.1

Released 2026-04-13.

Patch release: bug fix to `FixedAngleMode` and geometry documentation improvements.

### Breaking changes

- `FixedAngleMode` now holds the stage at its **current angle** (read from
  the geometry at solve time) rather than the value stored at construction.
  Preset the stage angle with `g.set_angle(name, value)` before calling
  `forward()`.  The `value` argument to `FixedAngleMode(...)` is now an
  initial default only.  This also applies to the `frozen_angles` dict on
  `BisectingMode` (e.g. `mu`/`nu` in psic and kappa6c).

### Fixed

- Geometry reference pages: reST `.. list-table::` directives inside MyST
  `.md` files rendered as raw source text in the browser — replaced with
  MyST-native pipe tables. (#139)
- `FixedAngleMode`: stored angle value was used by the solver instead of
  the stage's current angle, making it impossible to change the fixed value
  without replacing the mode object. (#142)

### Changed

- Geometry reference pages enriched: each mode now shows specific motor
  names, class AutoAPI links, GitHub source links, and plain-language
  constraint descriptions. (#140)
- Pre-built geometries table in `factories.py` module docstring sorted to
  match the geometry reference page order
  (fourcv → fourch → fivec → psic → sixc → kappas → surface). (#141)
- Geometry page titles use a consistent pattern:
  Eulerian Four-Circle, Eulerian Five-Circle, Eulerian Six-Circle,
  Kappa Four-Circle, Kappa Six-Circle, Z-Axis, S2D2. (#140)
- Geometry reference toctree reordered: fivec placed between four-circle
  and six-circle Eulerianss. (#140)
- `howto/modes.md` and `howto/forward.md` updated to show `g.set_angle()`
  pattern for fixed-angle modes. (#142)
- v1.0 criteria added to `roadmap.md`: #122 (diffraction modes) required;
  #114 (sphere of confusion) deferred. (#139)

## Release v0.4.0

Released 2026-04-13.

Major documentation release plus several API improvements.  Introduces a
complete User Guide (Concepts, How-to Guides, Geometry Reference, Glossary,
References), AutoAPI fixes, and a set of API quality improvements.

### Breaking changes

- `ReflectionList.setor1()` renamed to `setor0()` and `setor2()` renamed to
  `setor1()` to match the SPEC convention (primary = or0, secondary = or1).
  Update all call sites: `setor1("x")` → `setor0("x")`;
  `setor2("x")` → `setor1("x")`. (#120)
- `hkl_trajectory()`, `psi_trajectory()`, and `trajectory_plan()` are now
  **generators** (they `yield` one point dict at a time instead of returning a
  list).  Callers that iterate with `for pt in result:` are unaffected.
  Callers that index the result directly (e.g. `result[0]`) must wrap the
  call: `list(hkl_trajectory(...))`. (#115)
- `engines.py` renamed to `conversions.py`; all imports from
  `ad_hoc_diffractometer.engines` must be updated to
  `ad_hoc_diffractometer.conversions`. (#108)

### Added

- **User Guide**: Concepts overview, How-to Guides (wavelength, lattice,
  orient, forward/inverse, modes, trajectory, lattice refinement,
  serialization, crystal alignment), Geometry Reference (one page per
  pre-built geometry), Glossary, and References page. (#93, #94, #95)
- **API Reference**: module dependency graph (graphviz), per-layer card
  grid, AutoAPI module-prefix hidden via CSS. (#92, #102)
- `ahd.wh(geometry)` and `ahd.pa(geometry)` top-level convenience
  functions (SPEC-familiar status commands). (#131)
- `AdHocDiffractometer.stages_by_role(role)` — query stages with any
  arbitrary role string (enables analysers, polarisers, slits, etc.). (#129)
- How-to guide: save/restore diffractometer configuration to JSON or YAML
  via `to_dict()` / `from_dict()`. (#132)
- `sphinx-tabs` extension; coordinate-convention tabs in `concepts.md`
  showing You (1999), BL (1967), NeXus, and Hkl conventions side by side
  with the project default (You 1999) as the first tab. (#93)
- `dollarmath` MyST extension — inline `$...$` math now renders in `.md`
  files. (#120)

### Changed

- `Stage.role` is now documented as an **open string** — any value is
  accepted, not just `"sample"` and `"detector"`. (#129)
- Physical direction definitions (vertical, longitudinal, lateral)
  throughout the docs and docstrings now use unambiguous, observer-free
  language: vertical = opposite to gravity; longitudinal = chosen direction
  in the plane perpendicular to vertical; lateral = completes the
  right-handed system. (#128)
- `engines.py` renamed to `conversions.py` (see Breaking changes). (#108)
- Coordinate convention section in `concepts.md` reframed: observable
  physical directions first; Cartesian assignment shown in tabs. (#124, #125)
- AutoAPI: `sig-prename` (module-path prefix) hidden via `custom.css`
  so class/function signatures show short names only. (#102)
- `myst_enable_extensions` now includes `dollarmath` for inline math
  in Markdown files. (#120)
- Notebooks rendered by `myst-nb` (no pandoc required). (#95)
- `fourcv_alignment_howto.ipynb` moved from `docs/source/` to
  `docs/source/howto/`. (#94)
- `roadmap.md` moved from `docs/source/` to repo root (development-facing
  document). (#93)

### Fixed

- All 89 Sphinx AutoAPI warnings eliminated (duplicate descriptions,
  reST markup, ambiguous cross-references, module-prefix in signatures).
  (#92, #102)
- `switcher.json` automation: `latest` is always first; versions sorted
  newest-first. (#98)
- Inline `$...$` math in `concepts.md` now renders correctly online. (#120)

## Release v0.3.1

Released 2026-04-12.

Maintenance release: no user-facing API changes.  Exercises the new
versioned-docs and PyPI publishing automation introduced since v0.3.0.

### Added

- `docs_backfill.yml`: manual-trigger GitHub Actions workflow to build and
  publish documentation for any older tagged release and update
  `switcher.json` accordingly (#91)
- `pypi.yml`: GitHub Actions workflow to build and publish the package to
  PyPI on every version tag, using OIDC Trusted Publishing (no API token
  required) (#97)

### Changed

- `docs.yml`: on every tag push, automatically update
  `latest/_static/switcher.json` on the `gh-pages` branch — new stable
  releases are marked `preferred` and supersede pre-release entries; stale
  pre-release directories are pruned from `gh-pages` (#91)
- README now shows documentation and license badges (#84)

## Release v0.3.0

Released 2026-04-11.

### Fixed

- Stage axes corrected for `fourcv`, `fourch`, `kappa4cv`, `kappa4ch`;
  `v`/`h` suffix now consistently describes the scattering plane (#66)
- `two_theta` stage name renamed to `ttheta` throughout (#66)
- `Lattice.B` corrected to include the 2π factor per BL1967/SPEC convention (#78)
- Duplicate geometry name in entry-point registry now raises `ValueError` (#71)

### Added

- All factory functions accept `basis` as an optional argument; `BASIS_YOU`
  and `BASIS_BL` exported publicly; stage axes expressed as
  `VERTICAL`/`LATERAL`/`LONGITUDINAL` local aliases (#67)
- `logging.getLogger(__name__)` in all modules; silent by default (#61)
- 100% line and branch coverage enforced via `pytest-cov` (#60)
- Dynamic versioning via `hatch-vcs`; `ahd.__version__` exposed (#59)
- `to_dict()` / `from_dict()` on all major classes; JSON round-trip (#52)
- `g.wh()` and `g.pa()` status methods with `print=True` default (#51)
- `azimuthal_reference` property and `psi()` method (You 1999 definition) (#11)
- Entry-point plugin support for geometry factories;
  `GEOMETRY_ENTRY_POINT_GROUP` constant (#37)
- Diffraction modes / operating constraints: `mode_name`, `BisectingMode`,
  `FixedAngleMode`, `FixedAngleMode`; geometry accepts active mode (#9)
- `forward()`: hkl → motor angles; numerical multi-solution solver with
  seed-based bisection; returns all valid in-limits solutions (#35)
- Surface geometry: `alpha_i()`, `alpha_f()`, `q_components()`,
  `is_specular()`, `is_evanescent()` (#13)
- Detector geometry parameters: `detector_distance`, `detector_tilt`,
  `detector_offset` (#10)
- Diffractometer inclination: `inclination_matrix` property,
  `set_inclination()` (#15)
- Reciprocal-space trajectory computation: `hkl_trajectory()`,
  `psi_trajectory()`, `trajectory_plan()`, `NEAREST_ANGLES` solution-key
  callable; BL1967 operational ψ definition (#14)
- Alternative calculation engines — coordinate conversion functions:
  `hkl_to_Q()`, `Q_to_hkl()`, `Q_to_d()`, `d_to_Q_mag()`, `hkl_to_d()`,
  `d_to_two_theta()`, `two_theta_to_d()`, `hkl_to_two_theta()`,
  `two_theta_to_Q_mag()`, `Q_mag_to_two_theta()` (#12)
- X-ray radiation constants and conversions (NIST CODATA 2022):
  `HC_KEV_ANGSTROM` = 12.398 419 843 320 026 keV·Å (exact),
  `HC_KEV_ANGSTROM_UNCERTAINTY` = 0; `XRAY_LINES` named emission line
  wavelengths (Cu, Mo, Ag, Co Kα/Kα1/Kα2);
  `wavelength_to_energy()`, `energy_to_wavelength()`,
  `wavelength_to_wavenumber()`, `wavenumber_to_wavelength()`;
  `energy` and `wavenumber` lazy properties on `AdHocDiffractometer` (#21)
- Neutron source support: `source_type` property (`"xray"` default,
  `"neutron"` for fixed-wavelength reactor neutrons);
  `neutron_wavelength_to_energy()`, `neutron_energy_to_wavelength()`;
  `NEUTRON_MEV_ANGSTROM2` = 81.804 210 235 2 meV·Å² (CODATA 2022),
  `NEUTRON_MEV_ANGSTROM2_UNCERTAINTY`; `energy_units` property (`"keV"` or
  `"meV"`); `SOURCE_TYPES`; spallation/TOF out of scope (#8)
- Sphinx documentation framework: `docs/source/` with `conf.py`,
  `index.rst`, AutoAPI, myst-parser, pydata-sphinx-theme; GitHub Actions
  workflow `docs.yml` publishes to GitHub Pages on push to `main` or tag
  (#57)

### Changed

- `wh()` and `pa()` refactored from property to method with `print=True`
  keyword argument; `status.py` removed from the package API and moved to
  `references/` (#51)
- `spec.py` removed from the package API and moved to `references/` (#76)
- `summary()` on `AdHocDiffractometer` now reports energy (keV or meV) and
  its units alongside wavelength (#21, #8)
- Test files restructured to enforce one-to-one `test_<module>.py` naming;
  serialisation tests redistributed from `test_export.py` into their
  natural per-module homes and rewritten as parametrized functions (#88)

## Release v0.2

Released 2026-04-10.

### Added

- `Reflection` / `ReflectionList` classes: named reflections with hkl,
  motor angles, and wavelength; `setor1()` / `setor2()` designation;
  angle-key validation against geometry stages (#7)
- `Sample` / `SampleDict` classes: named samples with lattice, reflection
  list, U and UB matrices; active-sample guard on `SampleDict` (#25)
- Tolerance-aware `__eq__` for `Lattice` and `Reflection` (#29)
- `orientation.py` — UB matrix computation module:
  - `angles_to_phi_vector()`: motor angles → Q in the phi frame (#4)
  - `ub_identity()`: set U = I, UB = B
  - `ub_from_one_reflection()`: provisional UB from one reflection (#31)
  - `ub_from_two_reflections_bl1967()`: BL1967 eqs. 23-27 (#5)
  - `ub_from_three_reflections_bl1967()`: BL1967 eqs. 29-31 (#6)
- `AdHocDiffractometer.inverse()`: motor angles → (h, k, l) (#34)
- `Sample.parent`: back-reference from sample to geometry (#31)
- `spec.py` — SPEC fourc `#G1` line support (#26):
  - `FourcG1` named-tuple; `parse_fourc_g1()`; `emit_fourc_g1()`
  - `g1_to_sample()`: populate geometry from a `#G1` line
  - `sample_to_g1()`: emit current state as a `#G1` line
  - Verified against three historical lines from Align4Pete.spec
- `refinement.py` — lattice and orientation refinement (#32, #33):
  - `refine_lattice_bl1967()`: iterative least-squares (BL1967
    §Refinement); finite-difference Jacobian; no scipy required
  - `refine_lattice_simplex()`: Nelder-Mead derivative-free minimisation;
    uses scipy when available, pure-numpy fallback otherwise
  - `refine_all=False` (default): refines only the free parameters for
    the current crystal system, enforcing symmetry constraints at every
    iteration; `refine_all=True` treats all six parameters as independent
- `status.py` — SPEC-style status commands (#38):
  - `wh(geometry)`: terse position report — current HKL, λ, motor table
  - `pa(geometry)`: verbose parameters — geometry, reflections, lattice
    constants (real + reciprocal), λ; modelled on Align4Pete.log output

### Changed

- `AdHocDiffractometer` now holds a `SampleDict` (`.samples`) and an
  active-sample property (`.sample`), replacing the earlier single-sample
  design
- All new public symbols exported from `__init__.py` and listed in
  `__all__`

## Release v0.1

Released 2026-04-09.

### Added

- `AdHocDiffractometer` class: ordered rotary stages, basis validation,
  sample and detector rotation matrices, `summary()`, `check_limits()`
- `Stage` class: rotation axis, parent stacking, motor angle limits
  (`limits`, `in_limits()`)
- `wavelength` attribute on `AdHocDiffractometer` (Å, default `None`)
- `kappa_alpha_deg` property on `AdHocDiffractometer` (`None` for
  non-kappa geometries)
- Geometry factories: `psic`, `fourcv`, `fourch`, `sixc`, `kappa4cv`,
  `kappa4ch`, `kappa6c`, `zaxis`, `s2d2`, `fivec`
- `list_geometries()`, `get_geometry()`, `make_geometry()` registry API
- `Lattice` class: 7 crystal systems, B matrix, reciprocal and Cartesian
  lattice vectors, lazy computation, display precision
- `display.py`: `get_precision()`, `set_precision()`, `fmt()`
- GitHub Actions CI: Python 3.10–3.13 matrix; 3.14-dev as allowed failure
- Pre-commit hooks: ruff (lint + format), isort, standard file checks

### Changed

- Factory suffix convention: underscore removed from v/h detector suffix
  (`fourc_v` → `fourcv`, `kappa4c_h` → `kappa4ch`, etc.)
