# Change History

User-facing changes by release.  For future work planning, see [GitHub
Issues](https://github.com/BCDA-APS/ad_hoc_diffractometer/issues) for the full
issue tracker.  The initial project development roadmap is documented here:
[roadmap](https://github.com/BCDA-APS/ad_hoc_diffractometer/blob/main/roadmap.md).

## Unreleased


## Release v0.11.4

Released 2026-06-30

### Fixed

- psic `fixed_incidence_*` / `fixed_emergence_*` solutions now reproduce the requested hkl. (#307)
- psic `fixed_psi_vertical` no longer returns empty at the natural ψ. (#306)


## Release v0.11.3

Released 2026-06-24

### Fixed

- `ReferenceConstraint.evaluate()` computes real residuals; verified post-solve. (#304)


## Release v0.11.2

Released 2026-06-22

### Changed

- Adopt pre-v1.0 no-deprecation policy in AGENTS.md; pre-commit hook enforces. (#299)
- Rename `a_eq_b` reference constraint to `specular`. (#299)
- Rename `alpha_eq_beta_*` modes to `specular_*`. (#299)
- Rename `alpha_i` / `beta_out` reference constraints to `incidence` / `emergence`. (#299)
- Rename `azimuthal_reference` property and kwarg to `azimuth`. (#298)
- Rename `exit_angle` function to `emergence_angle`. (#299)
- Rename `fixed_alpha_i_*` / `fixed_beta_out_*` modes to `fixed_incidence_*` / `fixed_emergence_*`. (#299)
- Rename `geometry.alpha_i()` / `.alpha_f()` methods to `.incidence()` / `.emergence()`. (#299)
- Rename `surface.alpha_i()` / `.alpha_f()` functions to `.incidence()` / `.emergence()`. (#299)


## Release v0.11.1

Released 2026-06-09

### Added

- `AdHocDiffractometer.required_reference_vector` property. (#294)
- `ConstraintSet.with_constraint_values()` one-call value override. (#293)
- `register_geometry_yaml` in-memory registration entry point. (#288)
- Warn on `extras['n_hat']` assignment. (#294)

### Changed

- Disambiguate n̂ / surface_normal / azimuthal_reference / n_axis. (#294)
- Document how to override fixed-axis constraint values. (#293)

### Fixed

- Kappa equivalent-Eulerian chi axis now matches fourcv/fourch/psic. (#284)
- Populate alpha_i/beta_out/psi/omega output extras after `forward()`. (#292)

### Maintenance

- Unify copyright automation across hklpy2 family. (#290)


## Release v0.11.0

Released 2026-05-19


### Behavior change

- Rotation composition and `ub_identity` corrected. (#280)
- B matrix for non-cubic lattices now matches BL1967. (#280)

### Added

- Package now available on conda-forge.
- `ad_hoc_diffractometer.reference.natural_psi(g, h, k, l)`. (#278)

### Changed

- `forward()` warns when a `fixed_psi*` mode target ψ is unreachable. (#278)

### Fixed

- Doc version-switcher dropdown now shows every published version. (#273)
- Surface solver now respects detector pin in psic horizontal modes. (#279)

## Release v0.10.1

Released 2026-05-07

### Breaking changes

- Repository moved to <https://github.com/BCDA-APS/ad_hoc_diffractometer>

### Added

- Redirect docs site root to `latest/` on GitHub Pages. (#272)

## Release v0.10.0

Released 2026-05-06.

### Breaking changes

- Remove redundant psic modes. (#243)
- Switch to declarative YAML for demo geometries. (#267)

### Behavior change

- psic `fixed_chi_horizontal` switch `chi=90°` (kinematically infeasible) to `chi=0°`. (#259)
- psic `fixed_phi`/`fixed_chi` modes lock scattering-plane axis. (#243)

### Added

- `BASIS_DEFAULT` neutral basis used when a YAML file omits `basis:`. (#267)
- Declarative YAML geometry definitions and a loader. (#267)
- Glossary entries for zone-mode work (Azimuth, Basis, Constraint, Cut point, Double diffraction, Extras, Zone). (#266)
- psic and kappa6c zone modes (You §6, SPEC `setmode 5`). (#262)
- `register_geometry_file()` and `load_geometry_file()`. (#267)
- SPEC `OMEGA` pseudo-angle (`Q[6]`) as a `ReferenceConstraint` name. (#264)
- psic `fixed_omega_vertical` and `fixed_omega_horizontal` modes. (#264)
- psic `fixed_alpha_i_fixed_chi_fixed_phi` mode. (#264)
- psic `lifting_detector_eta` mode. (#264)

### Changed

- `pyyaml` is now a runtime dependency. (#267)
- psic `lifting_detector_phi` and `lifting_detector_mu` drop `qaz=90` and fix every sample stage except the named one. (#264)
- psic `fixed_psi_vertical` drops the bisect; pins `nu = 0` and `mu`. (#264)
- psic `fixed_psi_horizontal` drops the bisect; pins `delta = 0` and `eta`. (#264)

### Fixed

- RST markup leak in the inclination how-to rendered HTML. (#257)

## Release v0.9.3

Released 2026-05-04.

### Changed

- Re-focus the documentation for diffractometer-users (#255)

## Release v0.9.2

Released 2026-05-01.

### Fixed

- Restore correct kappa axis direction in kappa4cv, kappa4ch, kappa6c presets per published references (#252)

### Added

- Basis-invariance tests for all preset geometries (#252)

## Release v0.9.1

Released 2026-04-29.

### Behavior change

- Kappa pseudoangle conversions are now geometry-aware; re-derive saved kappa UB matrices. (#241)

### Added

- `KappaPseudoAngleConvention` and axis-aware kappa pseudoangle helpers. (#241)
- Cross-validation regression suite against `hkl_soleil` (libhkl). (#242)

### Fixed

- Kappa bisecting and fixed-virtual-angle solvers return solutions for physically reachable reflections. (#241)

## Release v0.9.0

Released 2026-04-28.

### Behavior change

- B / UB / `forward()` / `inverse()` results change for non-orthogonal
  lattices; re-derive saved UB matrices. (#237, #238)

### Fixed

- B matrix column layout for non-orthogonal lattices. (#237, #238)

## Release v0.8.0

Released 2026-04-28.

Performance and correctness release for the `forward()` solver.  The
slow-benchmark suite drops from ~561 s to ~56 s (~10× faster overall)
on top of the 3-4× speedup already delivered in v0.6.0.  See #217 for
the umbrella plan and the seven sub-issues #221–#227 for the
individual optimisations.

### Breaking changes

- Kappa **bisecting** modes now enforce *true* virtual bisecting
  (`omega_virtual = ttheta/2`, Walko 2016 eq. [16]) instead of the
  literal motor approximation `komega = ttheta/2`.  Affected modes:
  `kappa4cv.bisecting`, `kappa4ch.bisecting`,
  `kappa6c.bisecting_vertical`, `kappa6c.fixed_mu`, `kappa6c.fixed_nu`.
  Some `(h, k, l)` reflections previously accepted as "bisecting" are
  no longer accessible in these modes because their old solutions did
  not satisfy the physical bisecting condition.  Horizontal kappa
  bisect modes (`bisecting_horizontal`, `fixed_delta`) are unchanged.
  (#226, #235)

### Added

- `VirtualBisectConstraint` class — subclass of `BisectConstraint` that
  enforces a constraint on a *virtual* (Eulerian) angle computed from a
  kappa motor triple via `kappa_to_eulerian()`.  Exported from the
  top-level package namespace. (#226, #235)
- `rotation.rotation_matrix_and_derivative_normalized()` — returns both
  `R` and `dR/dθ` from a single trig evaluation, enabling closed-form
  Jacobians without redundant work. (#222, #229)
- `ForwardContext.jacobian_analytic()` — exact analytic Jacobian for
  the forward solver, computed by chain rule through the sample
  rotation matrix product. (#222, #229)
- Private functions and methods are now published in the AutoAPI
  reference (algorithm documentation); private attributes/data/
  properties remain hidden. (#231, #232)

### Changed (performance)

- **Analytic bisecting solver** for standard Eulerian geometries:
  `(chi, phi)` derived directly from `atan2`/`acos` formulas applied
  to the cached `Z_prefix` and `Q_lab` vectors.  **58–86× speedup**
  on bisecting modes for fourcv, fourch, psic, sixc, fivec.  Newton
  retained as fallback for kappa and degenerate cases. (#223, #230)
- **Decomposed double-diffraction solver**: 4D Newton-Raphson
  replaced by sequential subproblems (closed-form detector angle,
  1D outer-angle scan with sign-change detection, analytic 2D
  inner solve, scalar Ewald filter with 1D Newton refinement).
  **2.9–24×** per geometry; the slow-benchmark suite drops from
  561 s to 56 s (~10× overall). (#225, #234)
- **Analytic qaz detector solver**: 26-seed Newton loop replaced by
  closed-form solution.  Detector-angle solve goes from ~4.7 ms to
  ~0.3 µs (~13 500×); end-to-end `forward()` improves ~6 % since
  the inner sample-angle solver dominates. (#224, #233)
- **Analytic 1-free-angle solver**: when the only free sample stage
  rotates about a cardinal axis (±X, ±Y, ±Z), the 3×1 system
  reduces to a single `atan2`.  **8.4× speedup** on representative
  fourcv cases.  Falls back to Newton for tilted/kappa axes. (#227,
  #236)
- **Analytic Jacobian** in the generic Newton solver: closed-form
  derivative of the Rodrigues rotation matrix replaces the
  finite-difference Jacobian, eliminating `n_free` extra residual
  evaluations per Newton step.  **1.4–2.4× per Newton step** across
  every Newton path. (#222, #229)
- **Early termination in `solve_kappa_virtual()`**: matches the
  `_MAX_STALE` / `_MAX_SOLUTIONS` pattern from `_solve_bisecting`.
  **2–6×** speedup on kappa virtual-angle modes. (#221, #228)

### Fixed

- Kappa bisecting modes were returning solutions whose virtual omega
  did not equal ttheta/2 (the literal `komega = ttheta/2` only
  coincides with true bisecting at `kappa = 0`).  Now physically
  correct.  See breaking changes above. (#226, #235)
- Two docstring formatting issues in `forward.py` that produced RST
  warnings once private docstrings became visible in AutoAPI. (#232)

## Release v0.7.0

Released 2026-04-27.

### Breaking changes

- Standardized 14 mode names across 6 geometries to use the consistent
  `fixed_AXIS` naming convention (#220):
  - `constant_omega` → `fixed_omega` (fourcv, fourch, kappa4cv, kappa4ch)
  - `psi_constant` → `fixed_psi` (fourcv, fourch, kappa4cv, kappa4ch)
  - `constant_chi` → `fixed_chi` (kappa4cv, kappa4ch)
  - `constant_phi` → `fixed_phi` (kappa4cv, kappa4ch)
  - `constant_omega_noncoplanar` → `fixed_omega_noncoplanar` (fivec)
  - `mu_fixed` → `fixed_mu` (s2d2)
- No backward-compatibility aliases — mode names from v0.6.0 will raise
  `KeyError` if used directly

## Release v0.6.0

Released 2026-04-26.

### Breaking changes

- Renamed `LATERAL` to `TRANSVERSE` throughout (#193)
- Renamed `geometry.py` to `diffractometer.py` (#213)
- Removed legacy drawing functions; replaced by `StageAxisFigure`
  and `GeometryAxisFigure` (#209)

### Added

- Benchmark tool (`benchmark.py`) (#194)
- Complete psic mode set with vertical/horizontal suffixes and
  surface modes (#208)
- Interactive Plotly geometry diagrams (#197)
- HOWTO: custom geometry (#192), basis vectors (#191)
- Diagnostic logging in `drawing.py` (#207)
- Arc-direction handedness tests for all preset stages (#207)

### Changed

- `forward()` performance ~3-4x speedup (#195)
- Parent stage shown in geometry doc tables (#198)
- V/L/T axis labels in geometry diagrams (#204)
- Version switcher includes all releases; CI preserves live
  `switcher.json` on deploy (#216)

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
  arbitrary role string (enables analyzers, polarizers, slits, etc.). (#129)
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
- Physical direction definitions (vertical, longitudinal, transverse)
  throughout the docs and docstrings now use unambiguous, observer-free
  language: vertical = opposite to gravity; longitudinal = chosen direction
  in the plane perpendicular to vertical; transverse = completes the
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
  `VERTICAL`/`TRANSVERSE`/`LONGITUDINAL` local aliases (#67)
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
    constants (real + reciprocal), λ; modeled on Align4Pete.log output

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
