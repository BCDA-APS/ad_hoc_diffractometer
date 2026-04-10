# Change History

User-facing changes by release.  See [docs/roadmap.md](docs/roadmap.md)
for planned future work and
[GitHub Issues](https://github.com/prjemian/ad_hoc_diffractometer/issues)
for the full issue tracker.

## Unreleased

### Added

- Export/restore full diffractometer settings (#52): ``to_dict()`` /
  ``from_dict()`` on ``Lattice``, ``Reflection``, ``ReflectionList``,
  ``Sample``, and ``AdHocDiffractometer``; JSON-serialisable; top-level
  dict includes ``_meta`` with software name, version, and timestamp;
  complete round-trip preserves all stages, samples, UB matrices,
  reflections, and settings
- `AdHocDiffractometer.wh` and `AdHocDiffractometer.pa` properties (#51):
  access the terse/verbose status strings as ``g.wh`` / ``g.pa`` without
  needing to import the module-level functions; module-level ``wh(g)`` and
  ``pa(g)`` remain as thin wrappers for backward compatibility
- Azimuthal reference vector and ψ angle (#11):
  - ``AdHocDiffractometer.azimuthal_reference``: stores the reference
    direction as Miller indices (h, k, l); default ``None``; validated as
    a non-zero 3-vector
  - ``AdHocDiffractometer.psi(angles=None)``: computes the azimuthal angle
    ψ (You 1999 eqs. 10-11); ψ = 0 when the reference lies in the
    scattering plane; raises ``ValueError`` when reference ‖ Q
  - ``pa()`` shows the azimuthal reference; ``wh()`` shows a Psi line
- Entry-point extensibility for geometry factories (#37): all 10 built-in
  factories declared in ``pyproject.toml`` under the group
  ``"ad_hoc_diffractometer.geometries"``; third-party packages can
  contribute additional geometries without modifying this package.
  ``list_geometries()`` and ``get_geometry()`` discover plugins
  automatically.  ``GEOMETRY_ENTRY_POINT_GROUP`` constant exported.

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
