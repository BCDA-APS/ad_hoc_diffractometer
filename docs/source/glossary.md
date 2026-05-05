(glossary)=
# Glossary

Key terms used throughout `ad_hoc_diffractometer`, in alphabetical order.

```{glossary}

Azimuth
   An angle measured around a chosen reference axis, in the plane
   perpendicular to that axis.  In `ad_hoc_diffractometer` the relevant
   azimuth angles are **ψ** (rotation of a reference vector around **Q**;
   see *ψ angle*) and **qaz**/**naz** (lab-frame azimuth of **Q** or of the
   surface normal in the plane spanned by the vertical and transverse
   axes; You 1999, eq. 18).  See
   {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint` and
   {class}`~ad_hoc_diffractometer.mode.DetectorConstraint`.

B matrix
   The matrix that encodes the reciprocal lattice and maps Miller indices
   **h** = (h, k, l)ᵀ to the scattering vector in Cartesian crystal-frame
   coordinates: **Q**_c = B **h**.  Constructed from unit-cell parameters
   (a, b, c, α, β, γ) following Busing & Levy (1967), eq. 3.
   See {class}`~ad_hoc_diffractometer.lattice.Lattice` and
   {doc}`direct-lattice`.

Basis
   The three-vector dictionary that maps the physical direction names
   ``"vertical"``, ``"longitudinal"``, ``"transverse"`` to Cartesian unit
   vectors.  Each demo geometry carries a basis attribute; the package
   ships three named choices,
   {data}`~ad_hoc_diffractometer.factories.BASIS_YOU` (You 1999:
   vertical=+x, longitudinal=+y, transverse=+z),
   {data}`~ad_hoc_diffractometer.factories.BASIS_BL` (Busing & Levy
   1967: transverse=+x, longitudinal=+y, vertical=+z), and
   {data}`~ad_hoc_diffractometer.factories.BASIS_DEFAULT` — a neutral
   third assignment (vertical=+y, longitudinal=+z, transverse=+x) used
   by the declarative geometry loader as the fallback when a YAML file
   omits the ``basis:`` key.  Switching the basis re-expresses every
   stage axis in a different lab frame without changing the physics;
   see {doc}`problem2` for the basis-invariance case study.

Bisecting condition
   A diffraction-mode constraint in which a sample stage angle is driven
   to half the detector angle (e.g. ``omega = ttheta/2`` for fourcv,
   ``eta = delta/2`` for psic).  This places the sample symmetrically
   between the incident and diffracted beams.  Encoded as a
   {class}`~ad_hoc_diffractometer.mode.BisectConstraint`; on kappa
   geometries the physically correct version is the
   {class}`~ad_hoc_diffractometer.mode.VirtualBisectConstraint`, which
   bisects on the *virtual* Eulerian omega pseudoangle rather than the
   literal kappa motor.

Constraint
   An equation that fixes one of the geometry's free degrees of freedom
   during a forward calculation.  Four kinds are recognized:
   {class}`~ad_hoc_diffractometer.mode.SampleConstraint` (fix one sample
   angle), {class}`~ad_hoc_diffractometer.mode.DetectorConstraint` (fix
   one detector angle, including the special ``"qaz"`` pseudo-angle),
   {class}`~ad_hoc_diffractometer.mode.BisectConstraint` (relate one
   sample stage to one detector stage), and
   {class}`~ad_hoc_diffractometer.mode.ReferenceConstraint` (constrain a
   pseudo-angle between **Q** and an external reference vector).  A
   {class}`~ad_hoc_diffractometer.mode.ConstraintSet` is a fully-defined
   collection that resolves N − 3 free angles for an N-stage geometry.
   See {doc}`howto/constraints`.

Cut point
   A SPEC #G4 convention that selects the canonical 360°-window of an
   angle solution.  A cut-point ``C`` for stage X means returned values
   for that stage lie in ``[C, C + 360°)``.  Cut-points may be set
   per-geometry on the diffractometer or per-mode on a
   {class}`~ad_hoc_diffractometer.mode.ConstraintSet`; mode-level
   cut-points override geometry-level cut-points when both are present.

d-spacing
   The interplanar spacing d_{hkl} for the reflection (h, k, l), related to
   the scattering vector magnitude by |**Q**| = 2π / d_{hkl} and to the
   Bragg angle by 2d sin θ = λ.

Declarative geometry
   A diffractometer geometry described in a YAML file rather than as a
   Python factory function.  Loaded by
   {func}`~ad_hoc_diffractometer.geometry_loader.load_geometry_file` and
   {func}`~ad_hoc_diffractometer.geometry_loader.register_geometry_file`.
   Every declarative geometry file must declare the
   ``ad_hoc_diffractometer_geometry`` schema marker (see *Schema marker*
   below) and follow the format documented in
   {ref}`declarative-geometry-schema`.  The 10 demo geometries shipped
   under ``src/ad_hoc_diffractometer/geometries/`` are written in this
   form; user instruments should follow the same pattern.  See
   {ref}`howto-custom-geometry`.

Demo geometry
   One of the 10 YAML files shipped under
   ``src/ad_hoc_diffractometer/geometries/``: ``fourcv``, ``fourch``,
   ``psic``, ``sixc``, ``kappa4cv``, ``kappa4ch``, ``kappa6c``,
   ``zaxis``, ``s2d2``, ``fivec``.  These files are **demonstrations**
   of the declarative geometry schema, not authoritative descriptions
   of any production diffractometer.  Real instruments at real
   beamlines should copy the closest demo and adapt it.  Listed by
   {func}`~ad_hoc_diffractometer.factories.list_geometries`;
   instantiated by
   {func}`~ad_hoc_diffractometer.factories.make_geometry`.

Diffraction mode
   A named set of constraints that describes how
   {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward`
   computes motor angles: which stages are free, which are held fixed, and
   which are coupled to other stages.  Implemented as a
   {class}`~ad_hoc_diffractometer.mode.ConstraintSet`; a few specialized
   modes (double-diffraction, zone) are dispatched by their
   {attr}`~ad_hoc_diffractometer.mode.ConstraintSet.extras` schema rather
   than by their constraint composition.  See {doc}`howto/modes`.

Double diffraction
   A diffraction mode that solves for motor angles at which a primary
   reflection (h₁, k₁, l₁) and a secondary reflection (h₂, k₂, l₂)
   simultaneously satisfy the Ewald sphere condition.  Useful for
   measuring multi-beam interference effects.  The secondary reflection
   is supplied through the mode
   {attr}`~ad_hoc_diffractometer.mode.ConstraintSet.extras` dict
   (``h2``, ``k2``, ``l2``).  Available on psic and kappa6c as
   ``double_diffraction_vertical`` and ``double_diffraction_horizontal``;
   on the four-circle geometries simply as ``double_diffraction``.

Eulerian cradle
   An Eulerian cradle is a mechanical device used to hold and rotate objects
   (typically spheres or other symmetric bodies) in three-dimensional space
   along multiple axes simultaneously. It's named after the mathematician
   Leonhard Euler and his work on rotational motion. The device consists of
   a series of nested, rotating rings or gimbals mounted on each other. Each
   ring can rotate independently around a different axis (typically the x, y,
   and z axes). This allows an object mounted at the center to be oriented in
   any direction without mechanical singularities or constraints.

Extras
   The mode-attached dictionary of additional inputs and outputs needed
   beyond (h, k, l) for a forward calculation.  Inputs that the caller
   must supply (a secondary reflection for double diffraction, the
   reciprocal-lattice vectors that span a zone plane, a surface normal
   for surface modes) are pre-populated with the
   {data}`~ad_hoc_diffractometer.mode.REQUIRED` sentinel; outputs that
   the solver writes back (computed ψ, in-plane residual, etc.) are
   pre-populated with ``None``.  Stored as
   {attr}`ConstraintSet.extras <ad_hoc_diffractometer.mode.ConstraintSet>`.

Forward problem
   Given Miller indices (h, k, l) and a UB matrix, find the motor angles
   that satisfy the Bragg condition.  Returns a list of 0 to ~12 solutions
   depending on geometry and active diffraction mode.
   See {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.forward`
   and {doc}`howto/forward`.

Inverse problem
   Given motor angles and a UB matrix, find the unique Miller indices (h, k, l)
   in the Bragg condition.  Unlike the forward problem, the inverse always
   has a unique solution once UB is established.
   See {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.inverse`.

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

Schema marker (declarative geometry)
   The top-level mapping
   ``ad_hoc_diffractometer_geometry: {revision: <int>}`` that every
   declarative geometry YAML file must declare.  Its presence
   identifies the document as a geometry declaration; the integer
   ``revision`` selects which schema version the loader applies.
   Independent of the package version: schema revisions are deliberate
   editorial events that change far less often than package releases.
   The currently supported revision is
   {data}`~ad_hoc_diffractometer.geometry_loader.CURRENT_REVISION`.
   See {ref}`decl-marker`.

Sphere of confusion
   The small three-dimensional volume swept out by the nominal sample position
   as all diffractometer stages rotate.  Caused by mechanical tolerances;
   ideally a point but in practice a sphere of finite radius.

Stage
   A single rotary axis of the diffractometer.  Characterized by its axis
   direction (signed unit vector), parent stage, role (sample or detector),
   and optional motor limits.
   See {class}`~ad_hoc_diffractometer.stage.Stage`.

Serialization
   The process of converting the complete diffractometer state to a Python
   dict (via `to_dict()`) and restoring it (via `from_dict()`).  The dict
   contains only JSON-compatible types and can be saved to JSON (stdlib) or
   YAML (`pyyaml`).  All major classes support this round-trip.
   See {doc}`howto/serialize`.

Stack
   An ordered chain of stages in which each stage sits mechanically on the
   one below it (its parent).  The combined rotation matrix is the ordered
   product of individual stage matrices, from the base-mounted (outermost)
   stage to the innermost.

Two-theta
   The total scattering angle 2θ between the incident and diffracted
   beams, related to d-spacing and wavelength by Bragg's law:
   2 d sin θ = λ.  In simple geometries the detector arm rotates by 2θ
   directly, and the corresponding stage carries that name (``ttheta``
   on fourcv, fourch, kappa4cv, kappa4ch, fivec).  In multi-detector-
   axis geometries the detector position is a compound of two motors
   and 2θ is *computed* rather than directly motorised — for example
   psic, sixc, and kappa6c reach the angle ``2θ`` via the (nu, delta)
   pair (You 1999, eq. 17), and zaxis via (gamma, delta, alpha)
   (Walko 2016, eq. 17).  ``2θ`` itself is geometry-independent; the
   *stage name* that holds it (or computes to it) is not.

U matrix
   The orthonormal orientation matrix that relates the Cartesian crystal
   frame to the phi-axis frame.  U encodes how the crystal is mounted on the
   diffractometer; it is determined from one or more orienting reflections.
   U depends on the basis assignment for the diffractometer axes: see
   {doc}`concepts` for the definition and {doc}`problem2` for the basis-
   invariance case study.  See {doc}`howto/orient` for the orienting
   procedure.

UB matrix
   The product U × B.  Maps Miller indices directly to the phi-axis frame
   scattering vector: **Q**_φ = UB **h**.  Can be determined from measured
   reflections without prior knowledge of the unit-cell parameters.
   See {doc}`howto/orient`.

Vertical / horizontal scattering plane
   The plane containing the incident beam, the sample, and the diffracted
   beam.  *Vertical* — the scattering plane is vertical, the two-theta arm
   rotates about the transverse (horizontal) axis; typical for synchrotrons.
   *Horizontal* — the scattering plane is horizontal, two-theta rotates
   about the vertical axis; typical for laboratory instruments.

Zone (crystallographic)
   A set of crystal lattice planes that all share a common direction
   in direct space — the **zone axis** ``[u v w]``.  Equivalently, in
   reciprocal space, a set of reflections whose scattering vectors
   (the plane normals **g** = h **a*** + k **b*** + l **c***) are all
   perpendicular to that zone axis and therefore lie in a single
   plane — the *zone plane* — that passes through the reciprocal
   origin.  Any two non-parallel reciprocal-lattice vectors ``z0``
   and ``z1`` in the zone plane span it and uniquely identify the
   zone, and the *Weiss zone law* ``h u + k v + l w = 0``
   characterises membership in terms of the zone-axis indices.

Zone (mode)
   A diffraction mode (You 1999, §6; SPEC ``setmode 5``) named after
   the crystallographic zone above: the scattering vector **Q** is
   confined to the plane spanned by two reciprocal-lattice vectors
   ``z0`` and ``z1`` supplied through the mode
   {term}`extras <Extras>` dict.  A ``forward(h, k, l)`` call on a
   zone mode verifies that the requested **Q** lies in that plane and
   returns bisecting solutions; off-plane requests yield an empty
   list with a warning.  Available on psic and kappa6c as
   ``zone_vertical`` and ``zone_horizontal``; see the per-geometry
   pages for usage.

ψ (psi) angle
   Two definitions appear in the literature:
   (1) **You (1999)**: azimuthal angle of a reference vector about **Q** —
   a crystal-orientation diagnostic, constant for a given (hkl, UB).
   (2) **Busing & Levy (1967)**: angle of sample rotation about **Q** relative
   to a reference orientation — the quantity physically varied in a ψ scan.
   See {meth}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer.psi` and
   {func}`~ad_hoc_diffractometer.scan.psi_trajectory`.
```
