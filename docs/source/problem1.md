# Case Study: Describing a Diffractometer

This case study presents the diffractometer geometry problem that started the
`ad_hoc_diffractometer` project.  The equipment described here is a six-circle
diffractometer matching the You (1999) psic geometry; the analysis leads
directly to the {mod}`~ad_hoc_diffractometer.factories` factory functions and
the {class}`~ad_hoc_diffractometer.geometry.AdHocDiffractometer` class.

For a full analysis of how this geometry maps to the You (1999) convention,
see [You (1999) Geometry](problem2.md).

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

The reference frame is defined relative to the floor:

- **vertical**: positive direction is normal to the floor, pointing upward.
- **longitudinal**: positive direction is along the line of sight toward the
  equipment, normal to the vertical.
- **lateral**: positive direction is normal to both vertical and longitudinal,
  pointing to our left when facing the equipment.

## Equipment description

The equipment consists of two independent stage stacks that share a common
axis at their base.

**Stack 1 — detector stack** (2 stages)

- Stage 1 (base)
  - Axis of rotation: vertical
  - Sign of rotation: right-handed (consistent with coordinate system)
- Stage 2 (sits on stage 1)
  - Axis of rotation: lateral
  - Positive rotation: from longitudinal toward vertical
  - The detector is mounted on a long radial arm pointing at the sample

**Stack 2 — sample stack** (4 stages)

- Stage 1 (base)
  - Axis of rotation: vertical, colinear with stack 1 stage 1
  - Sign of rotation: same as stack 1 stage 1
- Stage 2 (sits on stage 1)
  - Axis of rotation: lateral
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

The answers to both questions are worked out in [You (1999) Geometry](problem2.md),
which shows how this equipment maps exactly to the six-circle psic geometry.
