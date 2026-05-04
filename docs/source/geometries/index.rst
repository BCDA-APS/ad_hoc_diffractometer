.. _geometries:

Demonstration Geometries
========================

Each geometry demonstrated below returns a fully-configured
:class:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer` instance.
Geometries are grouped below by their chi-circle mechanism.

.. _quick-start-demo:

Quick start
-----------

Pick one of the demo geometries below — for example the
:ref:`six-circle psic <geometry-psic>` instrument — set a wavelength,
attach a sample lattice, and print a summary of the diffractometer:

.. code-block:: python

   import ad_hoc_diffractometer as ahd

   # Use the six-circle psic demo geometry
   g = ahd.psic()
   g.wavelength = 1.0  # Å

   # Define the sample lattice (cubic silicon)
   g.sample.lattice = ahd.Lattice(a=5.431)

   # Show a summary of the diffractometer
   print(g.summary())

To build a four-circle diffractometer step by step — choosing a basis,
stacking stages, defining diffraction modes, and running a forward
calculation — without starting from a demo geometry, see the
:doc:`Quick Start guide </quick_start>`.

.. toctree::
   :hidden:

   fourcv
   fourch
   fivec
   psic
   sixc
   kappa4cv
   kappa4ch
   kappa6c
   zaxis
   s2d2

.. icons: https://fonts.google.com/icons

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Eulerian

      The chi circle is a full Eulerian cradle.

      - :ref:`geometry-fourch`
      - :ref:`geometry-fourcv`
      - :ref:`geometry-fivec`
      - :ref:`geometry-psic`
      - :ref:`geometry-sixc`

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Kappa

      The kappa stage is a replacement for the traditional chi circle.

      - :ref:`geometry-kappa4ch`
      - :ref:`geometry-kappa4cv`
      - :ref:`geometry-kappa6c`

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Inclination

      Geometries designed for surface diffraction or with fully decoupled sample
      and detector axes.

      - :ref:`geometry-s2d2`
      - :ref:`geometry-zaxis`
