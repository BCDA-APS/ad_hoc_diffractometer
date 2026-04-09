# Our equipment

We have a piece of equipment we'll use to examine an object.  The equipment will rotate the object in a 3-D space.  The equipment will also position a detector for our observations.

## Reference frame

We'll describe our reference frame for now with reference only to the floor.  

- vertical: The positive vertical direction is normal to the plane of the floor
  and points out of the floor.  
- longitudinal: The positive longitudinal direction is along our line of sight
  towards the equipment and normal to the vertical direction.
- lateral: The lateral direction is normal to both the vertical and
  longitudinal directions, which positive direction is to our left.
  
We can see the equipment by looking in the positive longitudinal direction,
normal to the lateral-vertical plane.

## Equipment Description

There is a piece of mechanical equipment consisting of various rotary stages.
No translational stages are involved. The rotational axes of all stages coincide
ideally at a single point of intersection in a cartesian space. The object to be
examined is mounted at this point of intersection.  

In practice, due to engineering precision, this point of intersection devolves
into a 3-D volume termed a "sphere of confusion".

The axes and sign for each of the rotary stages will be described when all
stages are at their 0 degree positions.

In its default orientation, with each axis at its 0 degree position, the
equipment can be subdivided into two independent stacks of stages.

The stack 1 sits on the floor and has two rotary stages.

- Stack 1
  - sits on the floor
  - stage 1
    - axis of rotation in the vertical direction
    - sign of motion consistent with coordinate system
  - stage 2
    - sits on stage 1
    - axis of rotation in the lateral direction
    - positive rotation (from the 0 degree position): positive longitudinal
      towards positive vertical
    - our detector is mounted on a long arm radial from this stage looking at the point of intersection
-Stack 2
  - sits on the floor
  - stage 1
    - axis of rotation in the vertical direction
    - axis is colinear with stack 1 stage 1
    - same sign of motion
  - stage 2
    - sits on stage 1
    - axis of rotation in the lateral direction
    - same sign of motion as stack 1 stage 2
  - stage 3
    - sits on stage 2
    - axis of rotation in the longitudinal direction
    - sign of motion consistent with coordinate system
  - stage 4
    - sits on stage 3
    - axis of rotation in the vertical direction
    - sign of motion consistent with coordinate system

## Problem 1

- Assign basis vectors (xHat, yHat, zHat) to each of the axes in our reference frame.
- Describe the orientation of each stack and stage in terms of these basis
  vectors, including the sign of rotation.
- Describe the steps to compute orientation matrix, U.

## Problem 2

- Is it possible to make different assignments of the basis vectors?
- What are the stage orientation vectors?
- How is the U matrix different?
