.. _howto:

How-to Guides
=============

Step-by-step guides for common tasks.  Each guide assumes you have
installed the package and are familiar with the
:ref:`concepts <concepts>`.

.. toctree::
   :hidden:

   wavelength
   lattice
   orient
   forward
   modes
   constraints
   trajectory
   refine_lattice
   serialize
   fourcv_alignment_howto

.. icons: https://fonts.google.com/icons

.. grid:: 2

   .. grid-item-card:: :material-outlined:`settings_input_antenna;3em` Set Wavelength / Energy
      :link: wavelength
      :link-type: doc

      Set the radiation wavelength or energy on a geometry, and convert
      between wavelength, energy, d-spacing, and Q.

   .. grid-item-card:: :material-outlined:`grain;3em` Define the Sample Lattice
      :link: lattice
      :link-type: doc

      Specify unit-cell parameters for any crystal system and inspect
      the resulting B matrix and reciprocal lattice.

   .. grid-item-card:: :material-outlined:`explore;3em` Orient a Crystal
      :link: orient
      :link-type: doc

      Compute the U and UB matrices from one, two, or three orienting
      reflections.

   .. grid-item-card:: :material-outlined:`calculate;3em` Solve the Forward Problem
      :link: forward
      :link-type: doc

      Find the motor angles that satisfy the Bragg condition for a
      given reflection (hkl → motor angles).

   .. grid-item-card:: :material-outlined:`tune;3em` Switch Diffraction Modes
      :link: modes
      :link-type: doc

      Choose which stages are free, fixed, or coupled during a
      forward calculation.

   .. grid-item-card:: :material-outlined:`rule;3em` Work with Constraints
      :link: constraints
      :link-type: doc

      Understand the constraint framework: DOF rule, constraint categories,
      custom modes, and the extras dict for advanced modes.

   .. grid-item-card:: :material-outlined:`route;3em` Plan a Trajectory
      :link: trajectory
      :link-type: doc

      Compute motor-angle sequences along a path through reciprocal
      space, or for a ψ scan.

   .. grid-item-card:: :material-outlined:`straighten;3em` Refine Lattice Constants
      :link: refine_lattice
      :link-type: doc

      Refine unit-cell parameters from measured Bragg peak positions using
      Busing & Levy (1967) least-squares or the Nelder-Mead simplex method,
      with and without crystal-system constraints.

   .. grid-item-card:: :material-outlined:`save;3em` Save and Restore Configuration
      :link: serialize
      :link-type: doc

      Save the complete diffractometer state to JSON or YAML and restore
      it in a later session.

   .. grid-item-card:: :material-outlined:`align_horizontal_center;3em` Align a Crystal
      :link: fourcv_alignment_howto
      :link-type: doc

      Full worked example: orient a sapphire crystal on a four-circle
      diffractometer (APS 7-ID-C session).
