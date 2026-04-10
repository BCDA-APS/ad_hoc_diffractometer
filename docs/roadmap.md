# ad_hoc_diffractometer — Development Roadmap

This document records parameters and features that are not yet implemented,
grouped by category and priority.  It is a living document; check items off
as they are completed.

---

## Status key

- [ ] not started
- [~] in progress
- [x] done

---

## Already implemented

- [x] Rotation axis description: signed basis vector notation (+x, -z, etc.)
- [x] Physical direction names (vertical, lateral, longitudinal) resolved
      against a caller-supplied basis dict
- [x] `kappa_axis(alpha_deg, basis)` for tilted kappa axes
- [x] `Stage` class: name, axis, parent, role, angle
- [x] `AdHocDiffractometer` class: stacking order, basis validation,
      sample/detector rotation matrices, summary
- [x] `@register_geometry` decorator and `list_geometries()` registry
- [x] Geometry factories: psic, fourcv, fourch, sixc, kappa4cv, kappa4ch,
      kappa6c, zaxis, s2d2, fivec
- [x] `Lattice` class: 7 crystal systems, lazy B matrix, reciprocal vectors,
      Cartesian lattice vectors, display precision
- [x] `display.py`: package-level and per-instance display precision
- [x] Pre-commit hooks: ruff, isort, standard checks

---

## Priority 1 — Near-term (needed for diffraction calculations)

### 1.1 Wavelength on AdHocDiffractometer ([#1](https://github.com/prjemian/ad_hoc_diffractometer/issues/1))

Units are fixed as Å to match unit cell edge lengths.  Energy conversions,
wave number, and named laboratory lines are out of scope here — see #21.

- [x] Add `wavelength: float | None` attribute to `AdHocDiffractometer`; default `None`
- [x] Validate on assignment: if not `None`, must be `> 0`; raise `ValueError` otherwise
- [x] `summary()` reports `λ = {value} Å` when set, `λ not set` when `None`
- [x] Display value via `display.fmt()`
- [x] Tests: default is `None`, valid assignment, invalid (≤ 0), display in `summary()`

### 1.2 Motor limits per stage ([#2](https://github.com/prjemian/ad_hoc_diffractometer/issues/2))

- [x] Add `limits: tuple[float, float]` attribute to `Stage`
      (min_angle, max_angle in degrees); default (-180, 180)
- [x] Add `Stage.in_limits(angle_deg)` method
- [x] Add `AdHocDiffractometer.check_limits(**angles)` method that
      verifies all supplied angles are within their stage limits
- [x] Validate limits at construction: min < max
- [x] Tests covering valid, invalid, and boundary cases

### 1.4 Split unit tests into module-specific files ([#17](https://github.com/prjemian/ad_hoc_diffractometer/issues/17))

- [x] Create `tests/test_axes.py`, `test_rotation.py`, `test_stage.py`, `test_geometry.py`, `test_factories.py`, `test_display.py`
- [x] Migrate all tests from `tests/test_diffractometer.py` into the appropriate new file
- [x] Remove `tests/test_diffractometer.py` once migration is complete
- [x] Verify full test suite still passes

### 1.5 GitHub Actions workflow for unit testing ([#18](https://github.com/prjemian/ad_hoc_diffractometer/issues/18))

- [x] Create `.github/workflows/tests.yml`
- [x] Trigger on `push` and `pull_request` to `main`
- [x] Matrix over all supported Python versions (`3.10`, `3.11`, `3.12`, `3.13`); `3.14-dev` included as an allowed failure for early warning
- [x] Install with `pip install -e .[dev]` and run `python -m pytest`

### 1.3 Kappa alpha queryable on the geometry instance ([#3](https://github.com/prjemian/ad_hoc_diffractometer/issues/3))

- [x] Store `kappa_alpha_deg` as an attribute on the `AdHocDiffractometer`
      instance returned by kappa factory functions (kappa4cv, kappa4ch,
      kappa6c); currently baked only into the axis vector
- [x] Add a `kappa_alpha_deg` property on the instance; `None` for non-kappa geometries
- [x] Tests verifying the value is correct and matches the axis vector

---

## Priority 2 — Medium-term (needed for UB formalism)

These items require motor-angle-to-phi-frame vector conversion, which must
be implemented first.

### 2.0 Motor-angle-to-phi-frame conversion ([#4](https://github.com/prjemian/ad_hoc_diffractometer/issues/4))

- [x] Function `angles_to_phi_vector(geometry, **motor_angles)`
      that computes the scattering vector in the phi-axis frame from a set
      of motor angles
- [x] This is the missing link needed for U and UB computation
- [x] Reference: Busing & Levy (1967), You (1999)

### 2.1 `ub_from_two_reflections_bl1967` (BL1967 eqs. 23-27) ([#5](https://github.com/prjemian/ad_hoc_diffractometer/issues/5))

- [x] Implement in `orientation.py`
- [x] Inputs: two `Reflection` objects + known lattice (B matrix)
- [x] Algorithm: Gram-Schmidt orthonormal triples Tc (crystal) and Tφ (phi
      frame); U = Tφ @ Tc.T; UB = U @ B
- [x] Sets `sample.U` and `sample.UB` in-place; returns UB

### 2.2 `ub_from_three_reflections_bl1967` (BL1967 eqs. 29-31) ([#6](https://github.com/prjemian/ad_hoc_diffractometer/issues/6))

- [x] Implement in `orientation.py`
- [x] Inputs: three `Reflection` objects; no prior lattice needed
- [x] Algorithm: UB = Hφ @ H⁻¹ (direct matrix inversion); U = UB @ B⁻¹
- [x] Sets `sample.U` and `sample.UB` in-place; returns UB

### 2.3 Orienting reflections data structure ([#7](https://github.com/prjemian/ad_hoc_diffractometer/issues/7))

- [x] `Reflection` dataclass: `name`, `hkl`, `angles` (keyed by stage name),
      `wavelength`, `geometry_name`; normalisation, validation, `__eq__`
- [x] `ReflectionList` class: ordered dict of named reflections with
      dict-like interface (`__getitem__`, `__delitem__`, `__contains__`,
      `__len__`, `__iter__`), `add()`, `remove()`, `clear()`
- [x] `ReflectionList.setor1()` / `setor2()`: designate primary/secondary
      orienting reflections; moving a reflection between slots clears the old
- [x] `ReflectionList.orienting_reflections` property: `[]`, `[or1]`, or
      `[or1, or2]`
- [x] Angle keys validated against the geometry's stage names at `add()` time
- [x] `AdHocDiffractometer.reflections` holds the `ReflectionList`
- [x] `AdHocDiffractometer.add_reflection()` convenience wrapper (inherits
      geometry wavelength)

### 2.4 SPEC #G1 format and reflection round-trip tests ([#26](https://github.com/prjemian/ad_hoc_diffractometer/issues/26))

- [x] Parse and emit the SPEC #G1 line format (hkl + motor angles + wavelength)
- [x] Tests: store a reflection → compute UB → verify UB @ or1_hkl ∥ Q_phi(or1)
- [x] Verify compatibility with alignment data in
      `references/2020-12-13-fourcc-alignment-7-id-c/`
- [x] Depends on: #7 (2.3), #5 (2.1 U matrix), #6 (2.2 UB matrix), #26 (2.5)

### 2.5 Sample dict on AdHocDiffractometer ([#25](https://github.com/prjemian/ad_hoc_diffractometer/issues/25))

- [x] New `Sample` class (`sample.py`): `name`, `lattice` (Lattice),
      `reflections` (ReflectionList), `U` (3×3 or None), `UB` (3×3 or None);
      `__eq__`, `__repr__`
- [x] `SampleDict` class: guarded ordered dict; rejects non-Sample values,
      blocks remove/replace/pop/clear of the active sample; `_samples` is a
      read-only property so the container itself cannot be replaced
- [x] Default sample: `"test"`, cubic lattice a=1 Å, empty reflections,
      U=None, UB=None
- [x] `AdHocDiffractometer.samples` — `SampleDict`, initialised with `"test"`
- [x] `AdHocDiffractometer.sample` — property for the active sample (get/set
      by name or object)
- [x] `AdHocDiffractometer.add_sample()` / `remove_sample()` — guarded CRUD
- [x] `AdHocDiffractometer.add_reflection()` delegates to the active sample's
      `ReflectionList`
- [x] U and UB matrices live on `Sample`, not on the geometry
- [x] `Lattice.__eq__` added (compares all six parameters, exact — see 2.6)
- [x] 37 new tests in `test_sample.py`

### 2.6 Tolerance-aware `__eq__` for Lattice and Reflection ([#29](https://github.com/prjemian/ad_hoc_diffractometer/issues/29))

- [x] Add `precision_atol(digits=None)` to `display.py`: returns
      `0.5 * 10**(-digits)` — half a unit in the last displayed decimal place
- [x] Add `allclose(a, b, atol=None, digits=None)` to `display.py`: wraps
      `np.allclose` with `rtol=0` and the precision-derived tolerance
- [x] `Lattice.__eq__`: pack six parameters into arrays → `np.allclose`
- [x] `Reflection.__eq__`: hkl and angle values → `np.allclose`; name and
      `geometry_name` → exact string comparison
- [x] Export `precision_atol` and `allclose` from `__init__.py`
- [x] Tests in `test_display.py`, `test_lattice.py`, `test_reflection.py`

### 2.7 `ub_from_one_reflection` — initial UB from one reflection ([#31](https://github.com/prjemian/ad_hoc_diffractometer/issues/31))

- [x] `Sample.parent` attribute: back-reference to owning
      `AdHocDiffractometer`; set by `add_sample()` and default sample
      construction; cleared to `None` by `remove_sample()`; excluded from
      `__eq__`; shown in `__repr__`
- [x] Create `orientation.py` module
- [x] `ub_identity(sample)`: sets U=I, UB=B; returns UB
- [x] `ub_from_one_reflection(sample, reflection, reference_hkl,
      reference_stage)`: Rodrigues rotation from crystal direction to lab
      axis; resolves `reference_stage` from Stage, str, or None+parent;
      sets `sample.U` and `sample.UB` in-place; returns UB
- [x] Edge cases: parallel (U=I), anti-parallel (perpendicular axis)
- [x] Export `ub_identity` and `ub_from_one_reflection` from `__init__.py`
- [x] Tests in `test_orientation.py` and `test_sample.py`

### 2.8 `refine_lattice_bl1967` — least-squares lattice refinement ([#32](https://github.com/prjemian/ad_hoc_diffractometer/issues/32))

- [x] Create `refinement.py` module (separate from `orientation.py`)
- [x] Implement `refine_lattice_bl1967(sample, reflections, ...)` using
      BL1967 §Refinement; simultaneous cell + orientation least-squares
- [x] `refine_all=False` (default): refine only the free params for the
      current crystal system; `refine_all=True`: all six params free
- [x] Returns result dict; updates `sample.lattice`, `sample.U`,
      `sample.UB` in-place

### 2.9 `refine_lattice_simplex` — derivative-free lattice refinement ([#33](https://github.com/prjemian/ad_hoc_diffractometer/issues/33))

- [x] Implement `refine_lattice_simplex(sample, reflections, ...)` in
      `refinement.py` using Nelder-Mead simplex; derivative-free
      alternative to BL1967 least squares
- [x] `refine_all` option (same semantics as `refine_lattice_bl1967`)
- [x] Same return dict structure as `refine_lattice_bl1967`

---

## Priority 3 — Later (instrument-specific, operational)

### 3.1 Diffraction mode / operating constraints ([#9](https://github.com/prjemian/ad_hoc_diffractometer/issues/9))

Modes are a first-class concern of the geometry description, not an
afterthought.  A mode specifies which degrees of freedom are
constrained during a diffraction calculation, and which are free.
The user should be able to declare available modes as part of the
geometry definition passed to `AdHocDiffractometer`, analogously to
how stages are declared.

**Interface design questions to resolve:**

- [ ] **Mode declaration in the geometry dict / constructor**: the
      `AdHocDiffractometer` constructor (or the factory functions)
      should accept a `modes` argument — a list or dict of named
      mode objects — so that the available modes are part of the
      geometry description, not set separately after construction.
      Example sketch:
      ```python
      AdHocDiffractometer(
          name="psic",
          stages=[...],
          modes={
              "bisecting":   BisectingMode(),
              "fixed_chi":   FixedAngleMode(stage="chi", value=90.0),
              "fixed_phi":   FixedAngleMode(stage="phi", value=0.0),
          },
          default_mode="bisecting",
      )
      ```
- [ ] **Mode object interface**: each mode must declare:
      - which stages it constrains (and to what values or relationships)
      - which stages remain free
      - any additional constraints (e.g. bisecting: ω = 2θ/2;
        reference-vector: Q || N)
      - a method `constrain(geometry, hkl) -> dict[stage_name, angle]`
        that computes the constrained angle set for a given (h, k, l)
- [ ] **Active mode**: `AdHocDiffractometer` holds a reference to the
      currently active mode; switching modes changes which constraints
      are applied during angle calculations
- [ ] **Geometry-specific modes**: factory functions should pre-populate
      the modes dict with the canonical modes for that geometry.
      For example, psic() should include bisecting, fixed-chi, fixed-mu,
      and reference-vector modes as defined by You (1999).
- [ ] **Mode name attribute**: current mode name accessible as a property
- [ ] **Frozen angle values**: stages held fixed in a given mode
      (matches SPEC #G0 frozen motor values; see spec_G_lines.md)
- [ ] **Cut points**: SPEC-style branch selection per axis determining
      which of the multiple valid angle solutions is chosen
      (matches SPEC #G4; see spec_G_lines.md)
- [ ] Reference: You (1999) section on operating modes and constraints;
      Walko (2016) section 3.3 "Operating Modes";
      Busing & Levy (1967) section on angle settings.
- [ ] Depends on: 2.0, 2.2 (angle calculations require UB matrix).

### 3.2 Azimuthal reference vector ([#11](https://github.com/prjemian/ad_hoc_diffractometer/issues/11))

- [ ] A reference direction (surface normal, crystal axis, or arbitrary
      vector) used to define the azimuthal angle ψ
- [ ] Needed for surface diffraction modes and azimuthal scans
- [ ] Reference: Busing & Levy (1967); You (1999) eqs. 10-11

### 3.3 Detector geometry parameters ([#10](https://github.com/prjemian/ad_hoc_diffractometer/issues/10))

- [ ] Sample-to-detector distance (mm or m)
- [ ] Detector tilt / offset angles (correction for non-ideal alignment)
- [ ] Relevant primarily when using area detectors rather than point detectors

### 3.4 Alternative calculation engines ([#12](https://github.com/prjemian/ad_hoc_diffractometer/issues/12))

The hkl (Miller index) formalism is the primary calculation engine
implemented here, following Busing & Levy (1967) and You (1999).  Users
have requested additional calculation engines:

- [ ] **Q-space engine**: calculate in terms of the scattering vector
      Q = (Qx, Qy, Qz) in Å⁻¹ rather than (h, k, l) in reciprocal
      lattice units.  Useful when the crystal structure is unknown or
      when doing diffuse scattering or liquid/amorphous studies.
      Relationship: Q = 2π · UB · h, so Q and hkl are interconvertible
      given UB.
- [ ] **d-spacing / 2θ engine**: express positions in terms of d-spacing
      (Å) and/or 2θ (degrees) directly.  Requires wavelength (item 1.1).
- [ ] **Reciprocal lattice unit (rlu) engine**: express Q in units of
      the reciprocal lattice vectors rather than Å⁻¹, i.e. fractional
      Miller indices without requiring integer h, k, l.
- [ ] Design consideration: engines should be pluggable — a common
      interface that can be swapped without changing the diffractometer
      geometry description.

### 3.5 Surface geometry: incidence and emergence angles ([#13](https://github.com/prjemian/ad_hoc_diffractometer/issues/13))

Surface diffraction geometries (zaxis, s2d2, sixc, psic in surface mode)
require calculation and control of additional angles beyond the standard
hkl mapping:

- [ ] **Angle of incidence (αi)**: angle between the incident beam and
      the sample surface.  For zaxis this is the alpha motor angle;
      for s2d2 it is the mu motor angle; for psic it must be computed
      from the motor angles and a surface normal reference vector.
- [ ] **Angle of emergence / exit (αf)**: angle between the diffracted
      beam and the sample surface.  A compound function of the detector
      motor angles and the surface normal.
- [ ] **In-plane / out-of-plane decomposition**: decompose Q into
      components parallel (Q‖) and perpendicular (Q⊥) to the sample
      surface.  Requires a surface normal reference vector (item 3.2).
- [ ] **Critical angle and evanescent wave conditions**: for GIXD,
      the incidence angle relative to the critical angle determines
      whether the beam is in total external reflection (surface-sensitive)
      or bulk-penetrating.  May be out of scope for this package but
      worth noting.
- [ ] **Specular condition**: αi = αf constraint, useful as an operating
      mode for reflectometry.
- [ ] Reference: Lohmeier & Vlieg (1993); You (1999) eqs. 10-11;
      Vlieg et al., J. Appl. Cryst. 20, 330-337 (1987).
- [ ] Raised by users; closely related to item 3.2 (azimuthal reference
      vector) and item 3.4 (diffractometer inclination).

### 3.6 Scans about an arbitrary reciprocal-space vector ([#14](https://github.com/prjemian/ad_hoc_diffractometer/issues/14))

Diffractometers are often used to scan motor angles such that the
scattering vector Q traces a specific path in reciprocal space relative
to a chosen reference direction.  Two important classes:

- [ ] **Arbitrary hkl vector scans**: given a reciprocal-space direction
      (h, k, l) and a starting point, compute the sequence of motor
      angles that keep Q on a specified trajectory — radially (along Q),
      transversely (perpendicular to Q), or along a crystal axis.
      Examples: L-scan, H-scan, radial scan, transverse scan.
      Requires the UB matrix (item 2.2) and a mode selection (item 3.1).
- [ ] **ψ (psi) scans**: rotation of the sample about the scattering
      vector Q while keeping the reflection condition satisfied.  The
      azimuthal angle ψ is the angle of rotation about Q (or d*).
      For a given hkl, ψ is varied by adjusting ω, χ, φ (or their
      kappa equivalents) while holding 2θ fixed.  Used to:
      - measure anisotropic absorption (ψ-scan absorption correction)
      - study crystal symmetry and orientation
      - probe the azimuthal dependence of diffracted intensity
        (resonant scattering, magnetic scattering)
      ψ = 0 is conventionally defined when the reference vector N
      (item 3.2) lies in the scattering plane.
- [ ] **Reciprocal-space trajectory planning**: given start and end
      points in (h, k, l) or (Qx, Qy, Qz), compute the motor-angle
      path and flag any portions that exceed motor limits (item 1.2)
      or enter inaccessible regions.
- [ ] Reference: Busing & Levy (1967), section on ψ scans;
      You (1999), azimuthal angle definition;
      ITC Vol. C, Sec. 2.2.6 (ψ scan via ω, χ, φ adjustment).
- [ ] Depends on: 2.0, 2.2, 3.1, 3.2.

### 3.7 Diffractometer inclination with respect to the incident beam ([#15](https://github.com/prjemian/ad_hoc_diffractometer/issues/15))

- [ ] Some instruments (or experimental configurations) mount the entire
      diffractometer at a non-zero angle relative to the incident beam
      direction — for example, tilted to access a specific range of
      incidence angles or to accommodate a grazing-incidence geometry
- [ ] This is distinct from the individual motor angles; it is a property
      of the overall instrument mounting
- [ ] Representation options to consider:
      - A single tilt angle (rotation about a specified axis)
      - A full 3×3 rotation matrix describing the lab-frame orientation
        of the diffractometer coordinate system relative to the beam
      - Euler angles or a quaternion
- [ ] The inclination would modify the effective basis vectors seen by the
      beam, and would need to be folded into the sample and detector
      rotation matrix products
- [ ] Reference: relevant for grazing-incidence X-ray diffraction (GIXD)
      and for instruments where the diffractometer is not aligned with
      the beam axis by default
- [ ] Raised by users; priority to be determined once Priority 1 and 2
      items are complete

### 3.7b `wh()` and `pa()` status commands ([#38](https://github.com/prjemian/ad_hoc_diffractometer/issues/38))

- [x] `wh(geometry)` — terse one-screen status: current HKL (via `inverse()`),
      wavelength, and motor-angle table (SPEC-style column names)
- [x] `pa(geometry)` — verbose parameter listing: geometry name, orienting
      reflections with angles and hkl, real- and reciprocal-space lattice
      constants, wavelength
- [x] Both functions return a `str`; graceful fallback when UB or wavelength
      not set; modelled on Align4Pete.log SPEC output
- [x] Exported from `__init__.py`; 35 tests in `test_status.py`

### 3.8 Energy / wave-number conversions and named radiation lines ([#21](https://github.com/prjemian/ad_hoc_diffractometer/issues/21))

Useful for reporting and cross-checking but not load-bearing for diffraction
angle calculations.  Depends on #1 (wavelength on `AdHocDiffractometer`).

- [ ] Compute photon energy E (keV) lazily: E = 12.39842 / λ (Å)
- [ ] Compute wave number k (Å⁻¹) lazily: k = 2π / λ
- [ ] `summary()` optionally reports E alongside λ
- [ ] Named constants for common laboratory lines:
      Cu Kα ≈ 1.5406 Å, Mo Kα ≈ 0.7107 Å, Ag Kα ≈ 0.5594 Å, Co Kα ≈ 1.7902 Å
- [ ] Tests covering conversions, lazy recomputation, and named constant values
- [ ] Related to: #8 (neutron source — independent, different formula)

### 3.9 Neutron radiation source support ([#8](https://github.com/prjemian/ad_hoc_diffractometer/issues/8))

Moved from Priority 1.  Wavelength storage (#1) is the only prerequisite;
this issue is independent of the X-ray energy conversion issue (#21).

- [ ] Add a `source_type` parameter (or separate subclasses) to distinguish
      `"xray"` and `"neutron"` radiation
- [ ] For neutrons: compute E (meV) lazily via de Broglie: E = 81.8042 / λ²
- [ ] `summary()` reports the correct energy units for each source type
- [ ] Validate that energy/wavelength conversions are not mixed across source types
- [ ] Tests covering both source types, unit correctness, and invalid cross-type usage
- [ ] Related to: #21 (X-ray energy conversion — independent, different formula)

---

## Notes

- Walko (2016) designation system: S = sample axes, D = detector axes,
  number = count, parentheses = shared base.  All implemented geometries
  are catalogued in factories.py.
- The distinction between laboratory (horizontal scattering plane, vertical
  detector) and synchrotron (vertical scattering plane, lateral detector)
  is encoded in the _v / _h suffix convention.
- The SPEC #G line format is documented in
  references/2020-12-13-fourcc-alignment-7-id-c/spec_G_lines.md.
- Threading and RunEngine diagnostics are in threading_diagnostics.md.
