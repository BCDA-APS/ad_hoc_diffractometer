# Change History

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
