# SPEC #G Control Line Interpretation — fourc geometry

> **Geometry-specific warning**: The `#G` line format — number of fields,
> field order, and field meaning — is **specific to each SPEC diffractometer
> geometry**.  The interpretation below applies only to the **fourc**
> geometry (Eulerian four-circle, Busing & Levy convention).  Other
> geometries (sixc, psic, kappa, surf, etc.) use different `#G` formats
> and field counts.  Do not apply this table to data files from other
> geometries without verifying the field layout independently.

Decoded from `Align4Pete.spec` by cross-referencing with `Align4Pete.log`
and the Busing & Levy (1967) four-circle diffractometer formalism.

The `#G` lines appear in the header of every scan in a SPEC data file.
They record the complete diffractometer geometry state at the time the
scan was started.  There are four lines: `#G0`, `#G1`, `#G3`, `#G4`.
(`#G2` is not used by the fourc geometry.)

---

## #G0 — Diffractometer mode and scan settings (27 values)

```
#G0 0 0 0 0 0 1 0 0 0 0 0 0 50 0 0.1 0 68 68 50 -1 1 1 3.13542 3.13542 0 463.6 838.8
```

| Index | Value | Meaning |
|-------|-------|---------|
|  0 | 0     | Diffractometer mode (0 = omega=0 / bisecting) |
|  1 | 0     | Sector (0 = default) |
|  2 | 0     | Freeze flag (0 = not frozen) |
|  3 | 0     | Omega-sector flag |
|  4 | 0     | F_ALPHA (frozen alpha angle, degrees) |
|  5 | 1     | F_BETA  (frozen beta angle, degrees) |
|  6 | 0     | F_OMEGA |
|  7 | 0     | F_AZIMUTH |
|  8 | 0     | F_THETA |
|  9 | 0     | F_PHI |
| 10 | 0     | F_CHI_Z |
| 11 | 0     | F_PHI_Z |
| 12 | 50    | Number of scan points (last scan) |
| 13 | 0     | CUT_AZI |
| 14 | 0.1   | Count time (seconds) |
| 15 | 0     | Omega-sector |
| 16 | 68    | Beam divergence parameter (horizontal) |
| 17 | 68    | Beam divergence parameter (vertical) |
| 18 | 50    | Scan points (repeated) |
| 19 | -1    | Scan range start |
| 20 |  1    | Scan range stop |
| 21 |  1    | Scan step |
| 22 | 3.13542 | Q value |
| 23 | 3.13542 | Q value |
| 24 | 0     | (unused) |
| 25 | 463.6 | Beamline parameter (energy / flux related) |
| 26 | 838.8 | Beamline parameter (energy / flux related) |

This line is **identical across all scans** in this file — the mode and
sector never changed during the session.

---

## #G1 — Lattice constants and orienting reflections (34 values)

Three distinct variants appear in this file, showing the progression of
the orientation as Walko refined the crystal alignment.

```
#G1 4.785 4.785 12.991 90 90 120
    1.516237713 1.516237713 0.483656786 90 90 60
    0 0 6   1 0 0
    41.94188 20.97 90 0   0 0
    60 30 0 0   0 0
    1.549802558 1.549802558   0 0
```

| Index | Meaning |
|-------|---------|
|  0–2  | Direct lattice **a, b, c** (Å) |
|  3–5  | Direct lattice **α, β, γ** (degrees) |
|  6–8  | Reciprocal lattice **a\*, b\*, c\*** (Å⁻¹, with 2π factor) |
|  9–11 | Reciprocal lattice **α\*, β\*, γ\*** (degrees) |
| 12–14 | Primary orienting reflection **h₀, k₀, l₀** |
| 15–17 | Secondary orienting reflection **h₁, k₁, l₁** |
| 18–21 | Primary angles: **2θ₀, θ₀, χ₀, φ₀** (degrees) |
| 22–23 | Unused (always zero) |
| 24–27 | Secondary angles: **2θ₁, θ₁, χ₁, φ₁** (degrees) |
| 28–29 | Unused (always zero) |
|    30 | **λ₀** — wavelength for primary reflection (Å) |
|    31 | **λ₁** — wavelength for secondary reflection (Å) |
| 32–33 | Unused (always zero) |

### Progression of #G1 across the session

| Scans | Primary hkl | Primary 2θ, θ, χ, φ | Secondary hkl | Secondary 2θ, θ, χ, φ | λ₀, λ₁ |
|-------|-------------|----------------------|---------------|------------------------|---------|
| 1–7   | (0,0,6) | 41.94188, 20.97, 90, 0 | (1,0,0) | 60, 30, 0, 0 | 1.5498, 1.5498 |
| 8–... | (0,0,6) | 41.9394, 20.3654, 89.32, 0 | (1,0,0) | 60, 30, 0, 0 | 1.5498, 1.5498 |
| later | (0,0,6) | 41.9394, 20.3654, 89.32, 0 | (1,0,4) | 35.39, 17.64, 50.89, 29.95 | 1.5498, 1.5498 |

The shift in the primary reflection angles between scans 1–7 and scan 8
onward reflects Walko entering the measured (006) peak position to replace
the initial calculated estimate.  The appearance of (1,0,4) as the
secondary reflection later in the session indicates a second reflection
was measured and entered.

---

## #G3 — UB matrix (9 values, row-major)

```
#G3 1.516237713 0.7581188565 2.961543674e-17
    5.523479432e-21 1.313100377 -7.934919163e-06
    -3.366724039e-16 2.15428495e-05 0.483656786
```

The 9 values are the UB matrix stored in row-major order:

```
UB = [[G3[0], G3[1], G3[2]],
      [G3[3], G3[4], G3[5]],
      [G3[6], G3[7], G3[8]]]
```

SPEC's UB matrix **includes the 2π factor**, so:

    UB @ hkl  gives the scattering vector q in Å⁻¹ (with 2π)
    |UB @ hkl| = 2π / d_hkl

Verified: `|UB @ [0,0,6]| = 6 × c* = 6 × 0.4837 = 2.902 Å⁻¹` ✓

Two distinct UB matrices appear in this file:

**Before** entering the measured (006) position (scans 1–7):
```
UB = [[ 1.516237713,   0.758118857,   ~0         ],
      [ ~0,            1.313100377,   ~0         ],
      [ ~0,            ~0,             0.483656786]]
```
This is nearly diagonal — consistent with the fake/calculated orientation
where the crystal is assumed to be perfectly aligned.

**After** entering the measured (006) position (scan 8 onward):
```
UB = [[ 1.516130941,   0.758065471,   0.005739700],
      [ 0.000189797,   1.313122226,  -0.005101257],
      [-0.017992647,   0.004854270,   0.483595823]]
```
The small but non-zero off-diagonal terms indicate the crystal's axes are
slightly misaligned from the diffractometer axes — the c-axis is not
exactly parallel to the φ-axis.

---

## #G4 — Current reciprocal space position and motor limits (26 values)

```
#G4 -1.571735014e-05 3.143470028e-05 5.999775333
    1.549802558
    20.97094 20.96931
    0 -180
    0 0 0 0 0 0 0 0
    -180 -170 -180 -180
    -180 -180 -180 -180 -180
    0
```

| Index | Meaning |
|-------|---------|
|  0–2  | Current reciprocal space position **H, K, L** (Miller indices) |
|   3   | Current wavelength **λ** (Å) |
|   4   | Current **α** — angle of incidence (degrees) |
|   5   | Current **β** — exit angle (degrees) |
|   6   | Current azimuthal angle (degrees) |
|   7   | CUT_AZI — azimuthal cut point (degrees) |
|  8–15 | Frozen motor values: F_ALPHA, F_BETA, F_OMEGA, F_AZIMUTH, F_THETA, F_PHI, F_CHI_Z, F_PHI_Z |
|  16   | CUT_TTH — cut point for 2θ (degrees) |
|  17   | CUT_TH  — cut point for θ  (degrees) |
|  18   | CUT_CHI — cut point for χ  (degrees) |
|  19   | CUT_PHI — cut point for φ  (degrees) |
| 20–24 | Additional cut points (azimuthal and extended motor limits) |
|  25   | Unused (always zero) |

Cut points define which branch of the multi-valued angle solutions SPEC
selects.  For example, CUT_TH = -170 means θ is kept in [-170°, +190°).
The defaults throughout this session are CUT_TTH = -180, CUT_TH = -170,
CUT_CHI = -180, CUT_PHI = -180.

The H, K, L values at indices 0–2 are the reciprocal space coordinates
of the **starting position** of the scan, not the scan centre.  These
are computed from the current motor angles and the UB matrix.

---

## Notes

- This entire document describes the **fourc** geometry only.  Field
  counts and meanings differ for other SPEC geometries (sixc, psic,
  kappa, surf, etc.).
- `#G2` does not appear in fourc geometry files.
- All `#G` lines are written at the **start of each scan**, capturing the
  complete geometry state at that moment.
- The UB matrix in `#G3` is recomputed by SPEC whenever an orienting
  reflection is entered or the lattice constants are changed.
- The reciprocal lattice parameters in `#G1` indices 6–11 are derived
  quantities computed from the direct lattice constants — they are stored
  for convenience and are redundant with indices 0–5.
