# Copyright (c) 2026 Pete R. Jemian <prjemian+ad_hoc_diffractometer@gmail.com>
# SPDX-License-Identifier: CC-BY-4.0
"""
Regression tests for issue #267 — declarative geometries.

For every demo geometry that has been migrated from a Python factory to
a declarative YAML file, this module verifies that the loader-built
geometry is **structurally equivalent** to a hand-built reference
geometry (the legacy Python factory).  This is the primary safety net
for the migration.

The current commit migrates only ``fourcv``.  Subsequent commits will
extend the parametrisation as the remaining nine demo geometries are
moved to YAML.
"""

from __future__ import annotations

import numpy as np
import pytest

import ad_hoc_diffractometer as ahd
from ad_hoc_diffractometer.diffractometer import AdHocDiffractometer
from ad_hoc_diffractometer.factories import BASIS_BL
from ad_hoc_diffractometer.factories import BASIS_YOU
from ad_hoc_diffractometer.factories import KAPPA_ALPHA_DEFAULT
from ad_hoc_diffractometer.kappa import KappaPseudoAngleConvention
from ad_hoc_diffractometer.kappa import kappa_axis_from_eulerian
from ad_hoc_diffractometer.mode import REQUIRED
from ad_hoc_diffractometer.mode import BisectConstraint
from ad_hoc_diffractometer.mode import ConstraintSet
from ad_hoc_diffractometer.mode import DetectorConstraint
from ad_hoc_diffractometer.mode import ReferenceConstraint
from ad_hoc_diffractometer.mode import SampleConstraint
from ad_hoc_diffractometer.mode import VirtualBisectConstraint
from ad_hoc_diffractometer.stage import Stage

# ---------------------------------------------------------------------------
# Hand-built reference geometries
# ---------------------------------------------------------------------------


def _reference_fourcv() -> AdHocDiffractometer:
    """Hand-built fourcv (the pre-#267 ``presets.fourcv`` body verbatim)."""
    basis = BASIS_BL
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -TRANSVERSE, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent=None, role="detector"),
    ]
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("omega", "ttheta")],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["omega", "chi", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name="fourcv",
        stages=stages,
        basis=BASIS_BL,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(vertical scattering plane, transverse ttheta, synchrotron)"
        ),
        modes=modes,
        default_mode="bisecting",
    )


def _reference_fourch() -> AdHocDiffractometer:
    """Hand-built fourch (the pre-#267 ``presets.fourch`` body verbatim)."""
    basis = BASIS_BL
    VERTICAL = basis["vertical"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("omega", -VERTICAL, parent=None, role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -VERTICAL, parent="chi", role="sample"),
        Stage("ttheta", -VERTICAL, parent=None, role="detector"),
    ]
    modes = {
        "bisecting": ConstraintSet(
            [BisectConstraint("omega", "ttheta")],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["omega", "chi", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["chi", "phi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["omega", "chi", "phi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name="fourch",
        stages=stages,
        basis=basis,
        description=(
            "Busing & Levy (1967) four-circle Eulerian diffractometer "
            "(horizontal scattering plane, vertical ttheta, laboratory)"
        ),
        modes=modes,
        default_mode="bisecting",
    )


def _reference_psic() -> AdHocDiffractometer:
    """Hand-built psic (the pre-#267 ``presets.psic`` body verbatim)."""
    basis = BASIS_YOU
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("eta", -TRANSVERSE, parent="mu", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="eta", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -TRANSVERSE, parent="nu", role="detector"),
    ]
    modes = {
        # ── Vertical scattering plane ───────────────────────────────────
        "bisecting_vertical": ConstraintSet(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
        ),
        "fixed_phi_vertical": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "delta"],
        ),
        "fixed_chi_vertical": ConstraintSet(
            [
                SampleConstraint("chi", 90.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "phi", "delta"],
        ),
        "fixed_alpha_i_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_beta_out_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("beta_out", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "alpha_eq_beta_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("a_eq_b", True),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_psi_vertical": ConstraintSet(
            [
                BisectConstraint("eta", "delta"),
                SampleConstraint("mu", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "fixed_alpha_i_fixed_chi_fixed_phi": ConstraintSet(
            [
                SampleConstraint("chi", 0.0),
                SampleConstraint("phi", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["mu", "eta", "nu", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_omega_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
                ReferenceConstraint("omega", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"omega": None},
        ),
        "double_diffraction_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        # ── Horizontal scattering plane ────────────────────────────────
        "bisecting_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
        ),
        "fixed_phi_horizontal": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "nu"],
        ),
        "fixed_chi_horizontal": ConstraintSet(
            [
                SampleConstraint("chi", 0.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "phi", "nu"],
        ),
        "fixed_alpha_i_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_beta_out_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("beta_out", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "alpha_eq_beta_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("a_eq_b", True),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_psi_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("eta", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "fixed_omega_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
                ReferenceConstraint("omega", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"omega": None},
        ),
        "double_diffraction_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        # ── Zone modes ─────────────────────────────────────────────────
        "zone_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["eta", "chi", "phi", "delta"],
            extras={"z0": REQUIRED, "z1": REQUIRED, "in_plane_residual": None},
        ),
        "zone_horizontal": ConstraintSet(
            [
                SampleConstraint("eta", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "chi", "phi", "nu"],
            extras={"z0": REQUIRED, "z1": REQUIRED, "in_plane_residual": None},
        ),
        # ── Lifting detector (out-of-plane) ────────────────────────────
        "lifting_detector_phi": ConstraintSet(
            [
                SampleConstraint("phi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["phi", "nu", "delta"],
        ),
        "lifting_detector_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("eta", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["mu", "nu", "delta"],
        ),
        "lifting_detector_eta": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("chi", 0.0),
                SampleConstraint("phi", 0.0),
            ],
            computed=["eta", "nu", "delta"],
        ),
    }
    return AdHocDiffractometer(
        name="psic",
        stages=stages,
        basis=BASIS_YOU,
        description=(
            "You (1999) 4S+2D six-circle diffractometer "
            "(transverse detector, vertical scattering plane, synchrotron)"
        ),
        modes=modes,
        default_mode="bisecting_vertical",
    )


def _reference_kappa4cv(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
) -> AdHocDiffractometer:
    """Hand-built kappa4cv (the pre-#267 ``presets.kappa4cv`` body verbatim).

    Reproduces the legacy factory exactly, including the synthesized
    :class:`KappaPseudoAngleConvention`.
    """
    basis = BASIS_BL
    TRANSVERSE = basis["transverse"]
    VERTICAL = basis["vertical"]
    kax = kappa_axis_from_eulerian(+TRANSVERSE, +VERTICAL, alpha_deg)
    stages = [
        Stage("komega", -TRANSVERSE, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -TRANSVERSE, parent="kappa", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent=None, role="detector"),
    ]
    convention = KappaPseudoAngleConvention(
        n_komega=-TRANSVERSE,
        n_kappa=kax,
        n_kphi=-TRANSVERSE,
        n_chi_eq=+VERTICAL,
    )
    modes = {
        "bisecting": ConstraintSet(
            [VirtualBisectConstraint("omega", "ttheta")],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_kphi": ConstraintSet(
            [SampleConstraint("kphi", 0.0)],
            computed=["komega", "kappa", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction": ConstraintSet(
            [],
            computed=["komega", "kappa", "kphi", "ttheta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
    }
    return AdHocDiffractometer(
        name="kappa4cv",
        stages=stages,
        basis=basis,
        description=(
            f"Four-circle kappa diffractometer, vertical scattering plane "
            f"(synchrotron). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        kappa_pseudo_angle_convention=convention,
        modes=modes,
        default_mode="bisecting",
    )


def _reference_kappa6c(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
) -> AdHocDiffractometer:
    """Hand-built kappa6c (the pre-#267 ``presets.kappa6c`` body verbatim).

    Reproduces the legacy factory exactly, including the synthesized
    :class:`KappaPseudoAngleConvention`.  Same kappa-axis convention
    as kappa4cv (tilted from +TRANSVERSE toward +VERTICAL); the only
    structural difference is the addition of the psic-style ``mu``
    sample outer axis and ``nu`` detector outer axis.
    """
    basis = BASIS_YOU
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    kax = kappa_axis_from_eulerian(+TRANSVERSE, +VERTICAL, alpha_deg)
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("komega", -TRANSVERSE, parent="mu", role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -TRANSVERSE, parent="kappa", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -TRANSVERSE, parent="nu", role="detector"),
    ]
    convention = KappaPseudoAngleConvention(
        n_komega=-TRANSVERSE,
        n_kappa=kax,
        n_kphi=-TRANSVERSE,
        n_chi_eq=+VERTICAL,
    )
    modes = {
        "bisecting_vertical": ConstraintSet(
            [
                VirtualBisectConstraint("omega", "delta"),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "bisecting_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("komega", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
        ),
        "fixed_kphi": ConstraintSet(
            [
                SampleConstraint("kphi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "delta"],
        ),
        "fixed_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                VirtualBisectConstraint("omega", "delta"),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "fixed_nu": ConstraintSet(
            [
                DetectorConstraint("nu", 0.0),
                VirtualBisectConstraint("omega", "delta"),
                SampleConstraint("mu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
        ),
        "fixed_delta": ConstraintSet(
            [
                DetectorConstraint("delta", 0.0),
                BisectConstraint("mu", "nu"),
                SampleConstraint("komega", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
        ),
        "lifting_detector_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("komega", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["mu", "nu", "delta"],
        ),
        "lifting_detector_kphi": ConstraintSet(
            [
                SampleConstraint("kphi", 0.0),
                SampleConstraint("mu", 0.0),
                DetectorConstraint("qaz", 90.0),
            ],
            computed=["kphi", "nu", "delta"],
        ),
        "fixed_psi_vertical": ConstraintSet(
            [
                VirtualBisectConstraint("omega", "delta"),
                SampleConstraint("mu", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "fixed_psi_horizontal": ConstraintSet(
            [
                BisectConstraint("mu", "nu"),
                SampleConstraint("komega", 0.0),
                ReferenceConstraint("psi", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
        "double_diffraction_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        "double_diffraction_horizontal": ConstraintSet(
            [
                SampleConstraint("komega", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
            extras={"h2": REQUIRED, "k2": REQUIRED, "l2": REQUIRED},
        ),
        "zone_vertical": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                DetectorConstraint("nu", 0.0),
            ],
            computed=["komega", "kappa", "kphi", "delta"],
            extras={"z0": REQUIRED, "z1": REQUIRED, "in_plane_residual": None},
        ),
        "zone_horizontal": ConstraintSet(
            [
                SampleConstraint("komega", 0.0),
                DetectorConstraint("delta", 0.0),
            ],
            computed=["mu", "kappa", "kphi", "nu"],
            extras={"z0": REQUIRED, "z1": REQUIRED, "in_plane_residual": None},
        ),
    }
    return AdHocDiffractometer(
        name="kappa6c",
        stages=stages,
        basis=basis,
        description=(
            f"Six-circle kappa diffractometer, psic-style outer axes "
            f"(transverse detector, synchrotron). "
            f"Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        kappa_pseudo_angle_convention=convention,
        modes=modes,
        default_mode="bisecting_vertical",
    )


def _reference_kappa4ch(
    alpha_deg: float = KAPPA_ALPHA_DEFAULT,
) -> AdHocDiffractometer:
    """Hand-built kappa4ch (the pre-#267 ``presets.kappa4ch`` body verbatim)."""
    basis = BASIS_BL
    VERTICAL = basis["vertical"]
    LONGITUDINAL = basis["longitudinal"]
    kax = kappa_axis_from_eulerian(+VERTICAL, +LONGITUDINAL, alpha_deg)
    stages = [
        Stage("komega", -VERTICAL, parent=None, role="sample"),
        Stage("kappa", kax, parent="komega", role="sample"),
        Stage("kphi", -VERTICAL, parent="kappa", role="sample"),
        Stage("ttheta", -VERTICAL, parent=None, role="detector"),
    ]
    convention = KappaPseudoAngleConvention(
        n_komega=-VERTICAL,
        n_kappa=kax,
        n_kphi=-VERTICAL,
        n_chi_eq=+LONGITUDINAL,
    )
    modes = {
        "bisecting": ConstraintSet(
            [VirtualBisectConstraint("omega", "ttheta")],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_kphi": ConstraintSet(
            [SampleConstraint("kphi", 0.0)],
            computed=["komega", "kappa", "ttheta"],
        ),
        "fixed_omega": ConstraintSet(
            [SampleConstraint("omega", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [SampleConstraint("chi", 90.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [SampleConstraint("phi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
        ),
        "fixed_psi": ConstraintSet(
            [ReferenceConstraint("psi", 0.0)],
            computed=["komega", "kappa", "kphi", "ttheta"],
            extras={"n_hat": REQUIRED, "psi": None},
        ),
    }
    return AdHocDiffractometer(
        name="kappa4ch",
        stages=stages,
        basis=basis,
        description=(
            f"Four-circle kappa diffractometer, horizontal scattering plane "
            f"(laboratory). Kappa alpha = {alpha_deg} deg."
        ),
        kappa_alpha_deg=alpha_deg,
        kappa_pseudo_angle_convention=convention,
        modes=modes,
        default_mode="bisecting",
    )


def _reference_sixc() -> AdHocDiffractometer:
    """Hand-built sixc (the pre-#267 ``presets.sixc`` body verbatim)."""
    basis = BASIS_YOU
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("omega", -TRANSVERSE, parent="alpha", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("delta", -TRANSVERSE, parent="alpha", role="detector"),
        Stage("gamma", +VERTICAL, parent="delta", role="detector"),
    ]
    modes = {
        "bisecting_4c": ConstraintSet(
            [
                SampleConstraint("alpha", 0.0),
                DetectorConstraint("gamma", 0.0),
                BisectConstraint("omega", "delta"),
            ],
            computed=["omega", "chi", "phi", "delta"],
        ),
        "fixed_gamma_5c": ConstraintSet(
            [
                DetectorConstraint("gamma", 0.0),
                SampleConstraint("alpha", 0.0),
                BisectConstraint("omega", "delta"),
            ],
            computed=["omega", "chi", "phi", "delta", "alpha"],
        ),
        "fixed_alpha_5c": ConstraintSet(
            [
                SampleConstraint("alpha", 0.0),
                BisectConstraint("omega", "delta"),
                DetectorConstraint("gamma", 0.0),
            ],
            computed=["omega", "chi", "phi", "delta", "gamma"],
        ),
        "fixed_alpha_zaxis": ConstraintSet(
            [
                SampleConstraint("alpha", 0.0),
                SampleConstraint("chi", 0.0),
                ReferenceConstraint("alpha_i", 0.0),
            ],
            computed=["omega", "delta", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "fixed_beta_zaxis": ConstraintSet(
            [
                DetectorConstraint("gamma", 0.0),
                SampleConstraint("chi", 0.0),
                ReferenceConstraint("beta_out", 0.0),
            ],
            computed=["omega", "delta", "alpha"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "alpha_eq_beta_zaxis": ConstraintSet(
            [
                SampleConstraint("chi", 0.0),
                SampleConstraint("phi", 0.0),
                ReferenceConstraint("a_eq_b", True),
            ],
            computed=["omega", "delta", "alpha", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
    }
    return AdHocDiffractometer(
        name="sixc",
        stages=stages,
        basis=basis,
        description=(
            "Lohmeier & Vlieg (1993) six-circle surface diffractometer "
            "(IUCr six-circle; Walko S(3D2)1). "
            "Sample and detector share the alpha (rotary table) base stage."
        ),
        modes=modes,
        default_mode="bisecting_4c",
    )


def _reference_zaxis() -> AdHocDiffractometer:
    """Hand-built zaxis (the pre-#267 ``presets.zaxis`` body verbatim)."""
    basis = BASIS_YOU
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("alpha", +VERTICAL, parent=None, role="sample"),
        Stage("Z", +LONGITUDINAL, parent="alpha", role="sample"),
        Stage("delta", -TRANSVERSE, parent="alpha", role="detector"),
        Stage("gamma", +VERTICAL, parent="delta", role="detector"),
    ]
    modes = {
        "zaxis": ConstraintSet(
            [ReferenceConstraint("alpha_i", 0.0)],
            computed=["Z", "delta", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
        "reflectivity": ConstraintSet(
            [ReferenceConstraint("a_eq_b", True)],
            computed=["Z", "delta", "alpha", "gamma"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
    }
    return AdHocDiffractometer(
        name="zaxis",
        stages=stages,
        basis=basis,
        description=(
            "Z-axis four-circle diffractometer (Bloch 1985; Walko 2016 (S1D2)1). "
            "Surface normal parallel to Z-axis. "
            "Sample and detector share the alpha (base) stage."
        ),
        modes=modes,
    )


def _reference_s2d2() -> AdHocDiffractometer:
    """Hand-built s2d2 (the pre-#267 ``presets.s2d2`` body verbatim)."""
    basis = BASIS_YOU
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("Z", +LONGITUDINAL, parent="mu", role="sample"),
        Stage("nu", +VERTICAL, parent=None, role="detector"),
        Stage("delta", -TRANSVERSE, parent="nu", role="detector"),
    ]
    modes = {
        "fixed_mu": ConstraintSet(
            [SampleConstraint("mu", 0.0)],
            computed=["Z", "nu", "delta"],
        ),
        "reflectivity": ConstraintSet(
            [ReferenceConstraint("a_eq_b", True)],
            computed=["mu", "Z", "nu", "delta"],
            extras={"n_hat": REQUIRED, "alpha_i": None, "beta_out": None},
        ),
    }
    return AdHocDiffractometer(
        name="s2d2",
        stages=stages,
        basis=basis,
        description=(
            "S2D2 four-circle diffractometer (Evans-Lutterodt & Tang 1995; "
            "Walko 2016 S2D2). "
            "Fully decoupled sample (mu, Z) and detector (nu, delta) axes."
        ),
        modes=modes,
    )


def _reference_fivec() -> AdHocDiffractometer:
    """Hand-built fivec (the pre-#267 ``presets.fivec`` body verbatim)."""
    basis = BASIS_YOU
    VERTICAL = basis["vertical"]
    TRANSVERSE = basis["transverse"]
    LONGITUDINAL = basis["longitudinal"]
    stages = [
        Stage("mu", +VERTICAL, parent=None, role="sample"),
        Stage("omega", -TRANSVERSE, parent="mu", role="sample"),
        Stage("chi", +LONGITUDINAL, parent="omega", role="sample"),
        Stage("phi", -TRANSVERSE, parent="chi", role="sample"),
        Stage("ttheta", -TRANSVERSE, parent="mu", role="detector"),
    ]
    modes = {
        "bisecting_4c": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                BisectConstraint("omega", "ttheta"),
            ],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_chi": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("chi", 90.0),
            ],
            computed=["omega", "phi", "ttheta"],
        ),
        "fixed_phi": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("phi", 0.0),
            ],
            computed=["omega", "chi", "ttheta"],
        ),
        "fixed_mu": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                BisectConstraint("omega", "ttheta"),
            ],
            computed=["omega", "chi", "phi", "ttheta"],
        ),
        "fixed_omega_noncoplanar": ConstraintSet(
            [
                SampleConstraint("mu", 0.0),
                SampleConstraint("omega", 0.0),
            ],
            computed=["mu", "chi", "phi", "ttheta"],
        ),
    }
    return AdHocDiffractometer(
        name="fivec",
        stages=stages,
        basis=basis,
        description=(
            "Five-circle diffractometer: fourcv on vertical mu base "
            "(Vlieg et al. 1987; Walko 2016 (S3D1)1). "
            "Sample and detector coupled through mu."
        ),
        modes=modes,
        default_mode="bisecting_4c",
    )


# ---------------------------------------------------------------------------
# Equivalence checker
# ---------------------------------------------------------------------------


def _assert_geometries_equivalent(
    declarative: AdHocDiffractometer,
    reference: AdHocDiffractometer,
) -> None:
    """Compare two geometries field-by-field for structural equivalence."""
    assert declarative.name == reference.name
    # basis
    for key in ("vertical", "longitudinal", "transverse"):
        np.testing.assert_array_equal(
            declarative.basis[key],
            reference.basis[key],
            err_msg=f"basis['{key}'] differs",
        )
    # stages
    decl_stages = list(declarative._stages.values())  # noqa: SLF001
    ref_stages = list(reference._stages.values())  # noqa: SLF001
    assert [s.name for s in decl_stages] == [s.name for s in ref_stages]
    for ds, rs in zip(decl_stages, ref_stages, strict=False):
        assert ds.name == rs.name
        np.testing.assert_array_equal(
            ds.axis,
            rs.axis,
            err_msg=f"stage {ds.name!r} axis differs",
        )
        assert ds.parent == rs.parent, f"stage {ds.name!r} parent differs"
        assert ds.role == rs.role, f"stage {ds.name!r} role differs"
    # modes
    assert list(declarative.modes) == list(reference.modes)
    for mname, decl_cs in declarative.modes.items():
        ref_cs = reference.modes[mname]
        assert decl_cs.constraints == ref_cs.constraints, (
            f"mode {mname!r} constraints differ: "
            f"{decl_cs.constraints!r} vs {ref_cs.constraints!r}"
        )
        assert decl_cs.computed == ref_cs.computed, f"mode {mname!r} computed differs"
        # extras: same keys, identical sentinel mapping for REQUIRED
        assert set(decl_cs.extras) == set(ref_cs.extras), (
            f"mode {mname!r} extras keys differ"
        )
        for ek in decl_cs.extras:
            assert decl_cs.extras[ek] is ref_cs.extras[ek] or (
                decl_cs.extras[ek] == ref_cs.extras[ek]
            ), f"mode {mname!r} extras[{ek!r}] differs"
    # default mode
    assert declarative.mode_name == reference.mode_name
    # Kappa attributes (None for non-kappa geometries; populated for kappa)
    assert declarative.kappa_alpha_deg == reference.kappa_alpha_deg, (
        "kappa_alpha_deg differs"
    )
    decl_conv = declarative.kappa_pseudo_angle_convention
    ref_conv = reference.kappa_pseudo_angle_convention
    if ref_conv is None:
        assert decl_conv is None, (
            "reference has no kappa_pseudo_angle_convention but declarative does"
        )
    else:
        assert decl_conv is not None, (
            "reference has a kappa_pseudo_angle_convention but declarative does not"
        )
        for axis_name in ("n_komega", "n_kappa", "n_kphi", "n_chi_eq"):
            np.testing.assert_allclose(
                getattr(decl_conv, axis_name),
                getattr(ref_conv, axis_name),
                atol=1e-12,
                err_msg=f"kappa_pseudo_angle_convention.{axis_name} differs",
            )


# ---------------------------------------------------------------------------
# Per-geometry equivalence tests
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name, reference_factory",
    [
        pytest.param("fourcv", _reference_fourcv, id="fourcv"),
        pytest.param("fourch", _reference_fourch, id="fourch"),
        pytest.param("psic", _reference_psic, id="psic"),
        pytest.param("kappa4cv", _reference_kappa4cv, id="kappa4cv"),
        pytest.param("kappa6c", _reference_kappa6c, id="kappa6c"),
        pytest.param("kappa4ch", _reference_kappa4ch, id="kappa4ch"),
        pytest.param("sixc", _reference_sixc, id="sixc"),
        pytest.param("zaxis", _reference_zaxis, id="zaxis"),
        pytest.param("s2d2", _reference_s2d2, id="s2d2"),
        pytest.param("fivec", _reference_fivec, id="fivec"),
    ],
)
def test_declarative_matches_reference(geom_name, reference_factory):
    """The declarative YAML build must equal the hand-built reference."""
    declarative = ahd.make_geometry(geom_name)
    reference = reference_factory()
    _assert_geometries_equivalent(declarative, reference)


@pytest.mark.parametrize(
    "geom_name",
    [
        pytest.param("fourcv", id="fourcv"),
        pytest.param("fourch", id="fourch"),
        pytest.param("psic", id="psic"),
        pytest.param("kappa4cv", id="kappa4cv"),
        pytest.param("kappa6c", id="kappa6c"),
        pytest.param("kappa4ch", id="kappa4ch"),
        pytest.param("sixc", id="sixc"),
        pytest.param("zaxis", id="zaxis"),
        pytest.param("s2d2", id="s2d2"),
        pytest.param("fivec", id="fivec"),
    ],
)
def test_declarative_basis_override_round_trips(geom_name):
    """Caller-supplied ``basis=`` reaches the constructed geometry."""
    g_default = ahd.make_geometry(geom_name)
    # Pick whichever of {BASIS_BL, BASIS_YOU} is *different* from the
    # file's declared default, so the override actually changes something.
    if all(
        np.array_equal(g_default.basis[k], BASIS_YOU[k])
        for k in ("vertical", "longitudinal", "transverse")
    ):
        other = BASIS_BL
    else:
        other = BASIS_YOU
    g_override = ahd.make_geometry(geom_name, basis=other)
    for k in ("vertical", "longitudinal", "transverse"):
        np.testing.assert_array_equal(g_override.basis[k], other[k])


# ---------------------------------------------------------------------------
# Kappa-specific overrides
# ---------------------------------------------------------------------------


@pytest.mark.parametrize(
    "geom_name, reference_factory",
    [
        pytest.param("kappa4cv", _reference_kappa4cv, id="kappa4cv"),
        pytest.param("kappa4ch", _reference_kappa4ch, id="kappa4ch"),
        pytest.param("kappa6c", _reference_kappa6c, id="kappa6c"),
    ],
)
@pytest.mark.parametrize(
    "alpha_deg",
    [
        pytest.param(30.0, id="alpha-30"),
        pytest.param(50.0, id="alpha-50-default"),
        pytest.param(67.5, id="alpha-67.5"),
    ],
)
def test_declarative_alpha_deg_override(geom_name, reference_factory, alpha_deg):
    """The kappa tilt angle override flows through the loader correctly."""
    declarative = ahd.make_geometry(geom_name, alpha_deg=alpha_deg)
    reference = reference_factory(alpha_deg=alpha_deg)
    _assert_geometries_equivalent(declarative, reference)
