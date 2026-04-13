.. _geometries:

Prebuilt Geometries
===================

Each factory function returns a fully configured
:class:`~ad_hoc_diffractometer.geometry.AdHocDiffractometer` instance.
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

Eulerian four-circle
--------------------

The chi circle is a full Eulerian cradle; omega and 2θ share the same
rotation axis.

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Synchrotron
      :link: fourcv
      :link-type: doc

      Vertical scattering plane — ω and 2θ rotate about the lateral axis.

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Laboratory
      :link: fourch
      :link-type: doc

      Horizontal scattering plane — ω and 2θ rotate about the vertical axis.

Eulerian five- and six-circle
-----------------------------

A four-circle Eulerian sample stack extended with one or two additional
base or detector stages.

.. grid:: 3

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Five-circle (fivec)
      :link: fivec
      :link-type: doc

      fourcv on a vertical mu base; sample and detector
      coupled through mu. Vlieg et al. (1987).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` 4S+2D (psic)
      :link: psic
      :link-type: doc

      Four sample stages (mu, eta, chi, and phi) and two independent detector
      stages (nu, delta). You (1999).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Shared base (sixc)
      :link: sixc
      :link-type: doc

      Sample and detector share a common alpha base stage.
      Lohmeier & Vlieg (1993).

Kappa
-----

The chi circle is replaced by a kappa axis tilted at α = 50° from the
vertical, giving a larger accessible volume and fewer mechanical
obstructions.

.. grid:: 3

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Four-circle synchrotron
      :link: kappa4cv
      :link-type: doc

      Vertical scattering plane.

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Four-circle laboratory
      :link: kappa4ch
      :link-type: doc

      Horizontal scattering plane.

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Six-circle
      :link: kappa6c
      :link-type: doc

      Psic-style outer axes (mu, nu) with kappa inner sample stages.

Surface / special
-----------------

Geometries designed for surface diffraction or with fully decoupled
sample and detector axes.

.. grid:: 2

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Z-axis
      :link: zaxis
      :link-type: doc

      Surface normal parallel to Z. Sample and detector share an alpha
      base stage. Bloch (1985).

   .. grid-item-card:: :material-outlined:`rotate_right;3em` Decoupled axes
      :link: s2d2
      :link-type: doc

      Two fully independent sample axes (mu, Z) and two detector axes
      (nu, delta). Evans-Lutterodt & Tang (1995).
