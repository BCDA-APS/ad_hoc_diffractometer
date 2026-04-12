(geometries)=
Geometry Reference
==================

Each factory function returns a fully configured
:class:`~ad_hoc_diffractometer.geometry.AdHocDiffractometer` instance.
Geometries are grouped below by their chi-circle mechanism.

.. toctree::
   :hidden:

   fourcv
   fourch
   psic
   sixc
   fivec
   kappa4cv
   kappa4ch
   kappa6c
   zaxis
   s2d2

.. icons: https://fonts.google.com/icons

Eulerian four-circle
--------------------

The chi circle is a full Eulerian cradle; omega and 2θ share the same
rotation axis.

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rotate_right;3em` fourcv
      :link: fourcv
      :link-type: doc

      Vertical scattering plane — synchrotron convention.
      ω and 2θ rotate about the lateral axis.
      ``ahd.fourcv()``

   .. grid-item-card:: :material-outlined:`rotate_right;3em` fourch
      :link: fourch
      :link-type: doc

      Horizontal scattering plane — laboratory convention.
      ω and 2θ rotate about the vertical axis.
      ``ahd.fourch()``

Eulerian six-circle
-------------------

A four-circle Eulerian sample stack with an additional mu (or shared alpha)
base stage and a two-axis detector arm.

.. grid:: 3

   .. grid-item-card:: :material-outlined:`rotate_right;3em` psic
      :link: psic
      :link-type: doc

      You (1999) 4S+2D. Four sample stages (mu, eta, chi, phi) and two
      detector stages (nu, delta).
      ``ahd.psic()``

   .. grid-item-card:: :material-outlined:`rotate_right;3em` sixc
      :link: sixc
      :link-type: doc

      Lohmeier & Vlieg (1993) surface geometry. Sample and detector share
      an alpha base stage.
      ``ahd.sixc()``

   .. grid-item-card:: :material-outlined:`rotate_right;3em` fivec
      :link: fivec
      :link-type: doc

      Vlieg et al. (1987). fourcv on a vertical mu base; sample and
      detector coupled through mu.
      ``ahd.fivec()``

Kappa
-----

The chi circle is replaced by a kappa axis tilted at α = 50° from the
vertical, giving a larger accessible volume and fewer mechanical
obstructions.

.. grid:: 3

   .. grid-item-card:: :material-outlined:`rotate_right;3em` kappa4cv
      :link: kappa4cv
      :link-type: doc

      Four-circle kappa, vertical scattering plane (synchrotron).
      ``ahd.kappa4cv()``

   .. grid-item-card:: :material-outlined:`rotate_right;3em` kappa4ch
      :link: kappa4ch
      :link-type: doc

      Four-circle kappa, horizontal scattering plane (laboratory).
      ``ahd.kappa4ch()``

   .. grid-item-card:: :material-outlined:`rotate_right;3em` kappa6c
      :link: kappa6c
      :link-type: doc

      Six-circle kappa with psic-style outer axes (mu, nu).
      ``ahd.kappa6c()``

Surface / special
-----------------

Geometries designed for surface diffraction or with fully decoupled
sample and detector axes.

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rotate_right;3em` zaxis
      :link: zaxis
      :link-type: doc

      Bloch (1985) Z-axis geometry. Surface normal parallel to Z.
      Sample and detector share an alpha base.
      ``ahd.zaxis()``

   .. grid-item-card:: :material-outlined:`rotate_right;3em` s2d2
      :link: s2d2
      :link-type: doc

      Evans-Lutterodt & Tang (1995). Two fully independent sample axes
      (mu, Z) and two detector axes (nu, delta).
      ``ahd.s2d2()``
