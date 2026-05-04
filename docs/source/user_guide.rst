.. _user_guide:

==========
User Guide
==========

.. toctree::
   :hidden:
   :caption: Overview

   about

.. toctree::
   :hidden:
   :caption: Tutorials

   quick_start

.. toctree::
   :hidden:
   :caption: How-to Guides

   howto/index

.. toctree::
   :hidden:
   :caption: Reference

   geometries/index
   glossary
   references

.. toctree::
   :hidden:
   :caption: Explanation

   concepts
   direct-lattice
   problem1
   problem2

.. icons: https://fonts.google.com/icons

Tutorials
---------

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rocket_launch;3em` Quick Start
      :link: quick_start
      :link-type: doc

      Build a four-circle diffractometer step by step — without a factory
      function — and run your first forward calculation.

How-to Guides
-------------

.. grid:: 2

   .. grid-item-card:: :material-outlined:`directions_run;3em` How-to Guides
      :link: howto/index
      :link-type: doc

      Step-by-step guides: set wavelength, define a lattice, orient a
      crystal, solve the forward problem, switch modes, plan trajectories,
      and align a crystal end-to-end.

Reference
---------

.. grid:: 2

   .. grid-item-card:: :material-outlined:`format_list_bulleted;3em` Prebuilt Geometries
      :link: geometries/index
      :link-type: doc

      Eulerian (4, 5, and 6 circle) · Kappa (4 and 6 circle) ·
      Surface and special-purpose

   .. grid-item-card:: :material-regular:`menu_book;3em` API Reference
      :link: api
      :link-type: doc

      Complete auto-generated reference for every public class,
      function, and constant.

   .. grid-item-card:: :material-outlined:`abc;3em` Glossary
      :link: glossary
      :link-type: doc

      Alphabetised definitions of key terms.

   .. grid-item-card:: :material-outlined:`library_books;3em` References
      :link: references
      :link-type: doc

      All literature citations — geometry papers, physical constants,
      and numerical methods.

Explanation
-----------

.. grid:: 2

   .. grid-item-card:: :material-outlined:`architecture;3em` Concepts
      :link: concepts
      :link-type: doc

      Coordinate conventions, axis sign convention, B/U/UB matrices,
      diffraction modes, and the ψ angle.

   .. grid-item-card:: :material-outlined:`functions;3em` Direct Lattice
      :link: direct-lattice
      :link-type: doc

      Vector mathematics and lattice vector conventions in crystallography.

   .. grid-item-card:: :material-outlined:`science;3em` Case Study
      :link: problem1
      :link-type: doc

      The diffractometer problem that started this project.

   .. grid-item-card:: :material-outlined:`calculate;3em` Choice of Basis and the UB Matrix
      :link: problem2
      :link-type: doc

      Two basis assignments applied to the same case-study equipment;
      the U/UB matrices differ by a fixed rotation but the physical
      angle ↔ (h, k, l) conversion is invariant.

Overview
--------

.. grid:: 2

   .. grid-item-card:: :material-outlined:`info;3em` About ad_hoc_diffractometer
      :link: about
      :link-type: doc

      What the package does, what you need to run it, and where it
      fits in your workflow.
