# Case Study: Describing a Diffractometer

This case study presents the diffractometer geometry problem that started the
`ad_hoc_diffractometer` project.  The equipment described here is a six-circle
diffractometer matching the You (1999) psic geometry; the analysis leads
directly to the demo geometries in {mod}`~ad_hoc_diffractometer.presets` and
the {class}`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer` class.

For a worked example showing that the choice of right-handed basis
applied to this equipment is arbitrary — different basis assignments
produce U/UB matrices that differ by a fixed rotation, but the
physical motor-angle ↔ (h, k, l) conversion is invariant — see
[Case Study: Choice of Basis and the UB Matrix](problem2.md).

---

## The equipment

We have a piece of mechanical equipment consisting of rotary stages only
(no translational stages).  The rotational axes of all stages ideally
coincide at a single point of intersection — the sample position.  In
practice, engineering tolerances cause this point to expand into a small
3-D volume known as the *sphere of confusion*.

The equipment rotates the sample in three-dimensional space and positions a
detector to observe the scattered radiation.  All stage angles are described
at their zero-degree positions.

## Reference frame

The reference frame is defined by three mutually perpendicular directions.
The positive sense of each is chosen to form a right-handed system:

- **vertical**: opposite to the direction of gravitational acceleration
  (i.e., upward).
- **longitudinal**: a chosen direction in the plane perpendicular to
  vertical, conventionally aligned with the nominal incident beam direction
  projected onto that plane, positive toward the equipment.  This is a
  property of the instrument installation, not of the beam itself.
- **transverse**: orthogonal to both vertical and longitudinal, with positive
  sense chosen so that the three directions form a right-handed system
  (vertical × longitudinal).

These directions are properties of the physical setup — they do not depend
on how basis vectors are assigned to Cartesian axes.  The mapping of
vertical, longitudinal, and transverse to x, y, z (or any other symbol) is a
separate choice made by the caller and encoded in a basis dict.  Two common
choices — You (1999) and Busing & Levy (1967) — are applied to this equipment
in {doc}`problem2`, where they are shown to produce U/UB matrices that
differ by a fixed rotation without changing the underlying physics.

## Equipment description

The equipment consists of two independent stage stacks that share a common
axis at their base.

**Stack 1 — detector stack** (2 stages)

- Stage 1 (base)
  - Axis of rotation: vertical
  - Sign of rotation: right-handed (consistent with coordinate system)
- Stage 2 (sits on stage 1)
  - Axis of rotation: transverse
  - Positive rotation: from longitudinal toward vertical
  - The detector is mounted on a long radial arm pointing at the sample

**Stack 2 — sample stack** (4 stages)

- Stage 1 (base)
  - Axis of rotation: vertical, colinear with stack 1 stage 1
  - Sign of rotation: same as stack 1 stage 1
- Stage 2 (sits on stage 1)
  - Axis of rotation: transverse
  - Sign of rotation: same as stack 1 stage 2
- Stage 3 (sits on stage 2)
  - Axis of rotation: longitudinal
  - Sign of rotation: right-handed (consistent with coordinate system)
- Stage 4 (sits on stage 3)
  - Axis of rotation: vertical
  - Sign of rotation: right-handed (consistent with coordinate system)

## Questions

1. Assign basis vectors (xHat, yHat, zHat) to each axis of the reference
   frame.  Describe the orientation of each stack and stage in terms of these
   basis vectors, including the sign of rotation.  Describe the steps to
   compute the orientation matrix U.

2. Is it possible to make different assignments of the basis vectors?  What
   are the resulting stage orientation vectors?  How does the U matrix differ?

Both questions are worked out in
[Case Study: Choice of Basis and the UB Matrix](problem2.md).  The
goal there is to show — using two different right-handed basis
assignments applied to the same physical equipment described above —
that the choice of basis is arbitrary: it changes the numerical
values of the axis vectors and the U/UB matrices by a fixed
rotation, but the physical conversion of motor angles ↔ (h, k, l)
is invariant.  The equipment also maps exactly to the six-circle
psic geometry.
