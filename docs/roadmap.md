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
- [x] Geometry factories: psic, fourc_v, fourc_h, sixc, kappa4c, kappa4c_h,
      kappa6c, zaxis, s2d2, fivec
- [x] `Lattice` class: 7 crystal systems, lazy B matrix, reciprocal vectors,
      Cartesian lattice vectors, display precision
- [x] `display.py`: package-level and per-instance display precision
- [x] Pre-commit hooks: ruff, isort, standard checks

---

## Priority 1 — Near-term (needed for diffraction calculations)

### 1.1 Wavelength / energy ([#1](https://github.com/prjemian/ad_hoc_diffractometer/issues/1))

- [ ] New class or module `radiation.py` (or attribute on a future
      `Experiment` class)
- [ ] Store wavelength λ (Å) as primary quantity; **default: Cu Kα, λ = 1.5406 Å**
- [ ] Compute photon energy E (keV) lazily: E = hc/λ = 12.39842 / λ
- [ ] Compute wave number k (Å⁻¹) lazily: k = 2π/λ
- [ ] Validate: λ > 0; reasonable range warning (e.g. 0.01 – 10 Å)
- [ ] Display precision via `display.fmt()`
- [ ] `__str__` reports λ and E
- [ ] Provide named constants for common laboratory lines
      (Cu Kα ≈ 1.5406 Å, Mo Kα ≈ 0.7107 Å, Ag Kα ≈ 0.5594 Å, Co Kα ≈ 1.7902 Å)
- [ ] Tests covering valid inputs, invalid inputs, lazy recomputation, and default value

### 1.1a Neutron radiation source support ([#8](https://github.com/prjemian/ad_hoc_diffractometer/issues/8))

- [ ] Add a `source_type` parameter (or separate subclasses) to distinguish `"xray"` and
      `"neutron"` radiation
- [ ] For neutrons: store λ (Å) as primary; compute E (meV) lazily via de Broglie:
      E = 81.8042 / λ²
- [ ] Ensure `__str__` reports the correct energy units for each source type
- [ ] Validate that energy/wavelength conversions are not mixed across source types
- [ ] Tests covering both source types, unit correctness, and invalid cross-type usage
- [ ] See also: 1.1 (X-ray wavelength/energy, issue #1)

### 1.2 Motor limits per stage ([#2](https://github.com/prjemian/ad_hoc_diffractometer/issues/2))

- [x] Add `limits: tuple[float, float]` attribute to `Stage`
      (min_angle, max_angle in degrees); default (-180, 180)
- [x] Add `Stage.in_limits(angle_deg)` method
- [x] Add `AdHocDiffractometer.check_limits(**angles)` method that
      verifies all supplied angles are within their stage limits
- [x] Validate limits at construction: min < max
- [x] Tests covering valid, invalid, and boundary cases

### 1.4 Split unit tests into module-specific files ([#17](https://github.com/prjemian/ad_hoc_diffractometer/issues/17))

- [ ] Create `tests/test_axes.py`, `test_rotation.py`, `test_stage.py`, `test_geometry.py`, `test_factories.py`, `test_display.py`
- [ ] Migrate all tests from `tests/test_diffractometer.py` into the appropriate new file
- [ ] Remove `tests/test_diffractometer.py` once migration is complete
- [ ] Verify full test suite still passes

### 1.5 GitHub Actions workflow for unit testing ([#18](https://github.com/prjemian/ad_hoc_diffractometer/issues/18))

- [ ] Create `.github/workflows/tests.yml`
- [ ] Trigger on `push` and `pull_request` to `main`
- [ ] Matrix over all supported Python versions (`3.10`, `3.11`, `3.12`, `3.13`)
- [ ] Install with `pip install -e .[dev]` and run `python -m pytest`

### 1.3 Kappa alpha queryable on the geometry instance ([#3](https://github.com/prjemian/ad_hoc_diffractometer/issues/3))

- [ ] Store `kappa_alpha_deg` as an attribute on the `AdHocDiffractometer`
      instance returned by kappa factory functions (kappa4c, kappa4c_h,
      kappa6c); currently baked only into the axis vector
- [ ] Add a `kappa_alpha` property or metadata dict to the instance
- [ ] Tests verifying the value is correct and matches the axis vector

---

## Priority 2 — Medium-term (needed for UB formalism)

These items require motor-angle-to-phi-frame vector conversion, which must
be implemented first.

### 2.0 Motor-angle-to-phi-frame conversion ([#4](https://github.com/prjemian/ad_hoc_diffractometer/issues/4))

- [ ] Function `angles_to_phi_vector(geometry, hkl_or_q, **motor_angles)`
      that computes the scattering vector in the phi-axis frame from a set
      of motor angles
- [ ] This is the missing link needed for U and UB computation
- [ ] Reference: Busing & Levy (1967), You (1999)

### 2.1 U matrix (orientation matrix) ([#5](https://github.com/prjemian/ad_hoc_diffractometer/issues/5))

- [ ] Implement Busing & Levy (1967) two-reflection algorithm (eqs. 23-27)
- [ ] Requires: two orienting reflections (hkl + motor angles), B matrix
- [ ] Returns orthonormal U (3×3)
- [ ] Depends on 2.0

### 2.2 UB matrix ([#6](https://github.com/prjemian/ad_hoc_diffractometer/issues/6))

- [ ] Compute UB = U @ B
- [ ] Also implement Busing & Levy (1967) three-reflection direct method
      (eqs. 29-31) for UB without known lattice parameters
- [ ] Depends on 2.1

### 2.3 Orienting reflections ([#7](https://github.com/prjemian/ad_hoc_diffractometer/issues/7))

- [ ] Data structure to store primary and secondary orienting reflections:
      hkl (Miller indices), motor angles at which the reflection was found,
      wavelength used
- [ ] Match the SPEC #G1 format (see spec_G_lines.md)
- [ ] Tests verifying round-trip: reflection → U → predicted angles

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
