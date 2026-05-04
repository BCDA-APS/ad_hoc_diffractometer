.. _geometries:

Demonstration Geometries
========================

Each geometry demonstrated below returns a fully-configured
:class:`~ad_hoc_diffractometer.diffractometer.AdHocDiffractometer` instance.
Geometries are grouped below by their chi-circle mechanism.

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
