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

### 1.1 Wavelength / energy

- [ ] New class or module `radiation.py` (or attribute on a future
      `Experiment` class)
- [ ] Store wavelength λ (Å) as primary quantity
- [ ] Compute photon energy E (keV) lazily: E = hc/λ = 12.39842 / λ
- [ ] Compute wave number k (Å⁻¹) lazily: k = 2π/λ
- [ ] Validate: λ > 0; reasonable range warning (e.g. 0.01 – 10 Å)
- [ ] Display precision via `display.fmt()`
- [ ] `__str__` reports λ and E
- [ ] Tests covering valid inputs, invalid inputs, lazy recomputation

### 1.2 Motor limits per stage

- [ ] Add `limits: tuple[float, float]` attribute to `Stage`
      (min_angle, max_angle in degrees); default (-180, 180)
- [ ] Add `Stage.in_limits(angle_deg)` method
- [ ] Add `AdHocDiffractometer.check_limits(**angles)` method that
      verifies all supplied angles are within their stage limits
- [ ] Validate limits at construction: min < max
- [ ] Tests covering valid, invalid, and boundary cases

### 1.3 Kappa alpha queryable on the geometry instance

- [ ] Store `kappa_alpha_deg` as an attribute on the `AdHocDiffractometer`
      instance returned by kappa factory functions (kappa4c, kappa4c_h,
      kappa6c); currently baked only into the axis vector
- [ ] Add a `kappa_alpha` property or metadata dict to the instance
- [ ] Tests verifying the value is correct and matches the axis vector

---

## Priority 2 — Medium-term (needed for UB formalism)

These items require motor-angle-to-phi-frame vector conversion, which must
be implemented first.

### 2.0 Motor-angle-to-phi-frame conversion

- [ ] Function `angles_to_phi_vector(geometry, hkl_or_q, **motor_angles)`
      that computes the scattering vector in the phi-axis frame from a set
      of motor angles
- [ ] This is the missing link needed for U and UB computation
- [ ] Reference: Busing & Levy (1967), You (1999)

### 2.1 U matrix (orientation matrix)

- [ ] Implement Busing & Levy (1967) two-reflection algorithm (eqs. 23-27)
- [ ] Requires: two orienting reflections (hkl + motor angles), B matrix
- [ ] Returns orthonormal U (3×3)
- [ ] Depends on 2.0

### 2.2 UB matrix

- [ ] Compute UB = U @ B
- [ ] Also implement Busing & Levy (1967) three-reflection direct method
      (eqs. 29-31) for UB without known lattice parameters
- [ ] Depends on 2.1

### 2.3 Orienting reflections

- [ ] Data structure to store primary and secondary orienting reflections:
      hkl (Miller indices), motor angles at which the reflection was found,
      wavelength used
- [ ] Match the SPEC #G1 format (see spec_G_lines.md)
- [ ] Tests verifying round-trip: reflection → U → predicted angles

---

## Priority 3 — Later (instrument-specific, operational)

### 3.1 Diffraction mode / operating constraints

- [ ] Mode name attribute on `AdHocDiffractometer`
      (e.g. "bisecting", "fixed-chi", "fixed-phi", "reference-vector")
- [ ] Frozen angle values: which stages are held fixed in a given mode
      (matches SPEC #G0 frozen motor values; see spec_G_lines.md)
- [ ] Cut points: SPEC-style branch selection per axis
      (matches SPEC #G4; see spec_G_lines.md)

### 3.2 Azimuthal reference vector

- [ ] A reference direction (surface normal, crystal axis, or arbitrary
      vector) used to define the azimuthal angle ψ
- [ ] Needed for surface diffraction modes and azimuthal scans
- [ ] Reference: Busing & Levy (1967); You (1999) eqs. 10-11

### 3.3 Detector geometry parameters

- [ ] Sample-to-detector distance (mm or m)
- [ ] Detector tilt / offset angles (correction for non-ideal alignment)
- [ ] Relevant primarily when using area detectors rather than point detectors

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
