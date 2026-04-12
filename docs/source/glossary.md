(glossary)=
# Glossary

Key terms used throughout `ad_hoc_diffractometer`, in alphabetical order.

```{glossary}

B matrix
   The matrix that encodes the reciprocal lattice and maps Miller indices
   **h** = (h, k, l)ᵀ to the scattering vector in Cartesian crystal-frame
   coordinates: **Q**_c = B **h**.  Constructed from unit-cell parameters
   (a, b, c, α, β, γ) following Busing & Levy (1967), eq. 3.
   See {class}`~ad_hoc_diffractometer.lattice.Lattice` and
   {doc}`direct-lattice`.

Bisecting condition
   A diffraction mode constraint in which the sample stage angle equals half
   the detector angle (ω = 2θ/2, or equivalently η = δ/2 for psic).  Places
   the sample symmetrically between the incident and diffracted beams.
   See {class}`~ad_hoc_diffractometer.mode.BisectingMode`.

d-spacing
   The interplanar spacing d_{hkl} for the reflection (h, k, l), related to
   the scattering vector magnitude by |**Q**| = 2π / d_{hkl} and to the
   Bragg angle by 2d sin θ = λ.

Diffraction mode
   A named set of constraints that describes how
   {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.forward`
   computes motor angles: which stages are free, which are held fixed, and
   which are coupled to other stages.
   See {class}`~ad_hoc_diffractometer.mode.DiffractionMode` and {doc}`howto/modes`.

Forward problem
   Given Miller indices (h, k, l) and a UB matrix, find the motor angles
   that satisfy the Bragg condition.  Returns a list of 0 to ~12 solutions
   depending on geometry and active diffraction mode.
   See {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.forward`
   and {doc}`howto/forward`.

Inverse problem
   Given motor angles and a UB matrix, find the unique Miller indices (h, k, l)
   in the Bragg condition.  Unlike the forward problem, the inverse always
   has a unique solution once UB is established.
   See {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.inverse`.

Miller indices
   Integer (or rational) triplet (h, k, l) that indexes a set of crystal
   lattice planes.  Used throughout as the reciprocal-space coordinate of a
   Bragg reflection.

Monochromatic radiation
   Radiation of a single, fixed wavelength λ.  All diffraction calculations
   in this package assume monochromatic radiation — a fundamental assumption
   of the single-wavelength four-circle formalism.
   Spallation (time-of-flight) sources are out of scope.

Orienting reflection
   A Bragg reflection measured at a known (h, k, l) whose measured motor
   angles are used to determine the U matrix.  Typically two reflections
   (or1 and or2) are used to compute UB via Busing & Levy (1967), eqs. 23–27.
   See {func}`~ad_hoc_diffractometer.ub_from_two_reflections_bl1967` and
   {doc}`howto/orient`.

Phi frame
   The Cartesian reference frame attached to the innermost sample stage (φ
   for Eulerian geometries, kφ for kappa geometries).  The UB matrix maps
   Miller indices to scattering vectors expressed in the phi frame.

Q (scattering vector)
   The momentum transfer vector **Q** = **k**_f − **k**_i, where **k**_i and
   **k**_f are the incident and scattered wave vectors.  For elastic scattering
   |**Q**| = (4π/λ) sin θ = 2π / d_{hkl}.

Sphere of confusion
   The small three-dimensional volume swept out by the nominal sample position
   as all diffractometer stages rotate.  Caused by mechanical tolerances;
   ideally a point but in practice a sphere of finite radius.

Stage
   A single rotary axis of the diffractometer.  Characterised by its axis
   direction (signed unit vector), parent stage, role (sample or detector),
   and optional motor limits.
   See {class}`~ad_hoc_diffractometer.stage.Stage`.

Stack
   An ordered chain of stages in which each stage sits mechanically on the
   one below it (its parent).  The combined rotation matrix is the ordered
   product of individual stage matrices, from the floor-mounted (outermost)
   stage to the innermost.

Two-theta
   The total scattering angle 2θ between the incident and diffracted beams,
   related to d-spacing and wavelength by Bragg's law: 2d sin θ = λ.

U matrix
   The orthonormal orientation matrix that relates the Cartesian crystal
   frame to the phi-axis frame.  U encodes how the crystal is mounted on the
   diffractometer; it is determined from one or more orienting reflections.
   See {doc}`howto/orient` and {doc}`problem2`.

UB matrix
   The product U × B.  Maps Miller indices directly to the phi-axis frame
   scattering vector: **Q**_φ = UB **h**.  Can be determined from measured
   reflections without prior knowledge of the unit-cell parameters.
   See {doc}`howto/orient`.

Vertical / horizontal scattering plane
   The plane containing the incident beam, the sample, and the diffracted
   beam.  *Vertical* — the scattering plane is vertical, the two-theta arm
   rotates about the lateral (horizontal) axis; typical for synchrotrons.
   *Horizontal* — the scattering plane is horizontal, two-theta rotates
   about the vertical axis; typical for laboratory instruments.

ψ (psi) angle
   Two definitions appear in the literature:
   (1) **You (1999)**: azimuthal angle of a reference vector about **Q** —
   a crystal-orientation diagnostic, constant for a given (hkl, UB).
   (2) **Busing & Levy (1967)**: angle of sample rotation about **Q** relative
   to a reference orientation — the quantity physically varied in a ψ scan.
   See {meth}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer.psi` and
   {func}`~ad_hoc_diffractometer.psi_trajectory`.
```
