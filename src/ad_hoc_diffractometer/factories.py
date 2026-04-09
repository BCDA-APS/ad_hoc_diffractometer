"""
factories.py — predefined diffractometer geometry factory functions.

Each factory returns a fully configured AdHocDiffractometer instance for
a well-known diffractometer geometry.  Add new geometries here as the
package grows.
"""

import numpy as np

from .constants import XHAT
from .constants import YHAT
from .constants import ZHAT
from .geometry import AdHocDiffractometer
from .stage import Stage


def geometry_psic() -> AdHocDiffractometer:
    """
    You (1999) '4S+2D' six-circle diffractometer (psic geometry).

    Basis: xHat=vertical, yHat=longitudinal, zHat=lateral.

    Sample stack (floor first):
        mu  (S2-1): vertical,     +x, right-handed
        eta (S2-2): lateral,      -z, left-handed
        chi (S2-3): longitudinal, +y, right-handed
        phi (S2-4): lateral,      -z, left-handed

    Detector stack (floor first):
        nu    (S1-1): vertical, +x, right-handed
        delta (S1-2): lateral,  -z, left-handed

    mu and nu share the same vertical rotation axis and are mechanically
    independent (same axis, not coupled).

    Reference: H. You, J. Appl. Cryst. 32, 614-623 (1999).
               DOI: 10.1107/S0021889899001223
    """
    stages = [
        Stage("mu", +XHAT, parent=None, role="sample"),
        Stage("eta", -ZHAT, parent="mu", role="sample"),
        Stage("chi", +YHAT, parent="eta", role="sample"),
        Stage("phi", -ZHAT, parent="chi", role="sample"),
        Stage("nu", +XHAT, parent=None, role="detector"),
        Stage("delta", -ZHAT, parent="nu", role="detector"),
    ]
    return AdHocDiffractometer(
        name="psic",
        stages=stages,
        description="You (1999) 4S+2D six-circle diffractometer",
    )


def geometry_fourc() -> AdHocDiffractometer:
    """
    Busing & Levy (1967) four-circle Eulerian diffractometer.

    Basis (Busing & Levy convention):
        x = lateral  (scattering vector at zero angles)
        y = longitudinal (along the beam)
        z = vertical
    Right-handed: lateral x longitudinal = vertical => +x x +y = +z.

    Sample stack (floor first):
        omega:     vertical, -z, left-handed
        chi:       lateral,  +x, right-handed  (NB: lateral = x in this basis)
        phi:       vertical, -z, left-handed

    Detector (floor, independent of sample stack):
        two_theta: vertical, -z, left-handed

    omega and two_theta share the same vertical axis and are mechanically
    independent (same relationship as S1-1/S2-1 in our equipment description).

    Reference: W.R. Busing & H.A. Levy, Acta Cryst. 22, 457-464 (1967).
    """
    ZHAT_BL = np.array([0.0, 0.0, 1.0])  # vertical in Busing & Levy frame
    XHAT_BL = np.array([1.0, 0.0, 0.0])  # lateral  in Busing & Levy frame
    basis = {
        "lateral": XHAT_BL,
        "longitudinal": np.array([0.0, 1.0, 0.0]),
        "vertical": ZHAT_BL,
    }
    stages = [
        Stage("omega", -ZHAT_BL, parent=None, role="sample"),
        Stage("chi", +XHAT_BL, parent="omega", role="sample"),
        Stage("phi", -ZHAT_BL, parent="chi", role="sample"),
        Stage("two_theta", -ZHAT_BL, parent=None, role="detector"),
    ]
    return AdHocDiffractometer(
        name="fourc",
        stages=stages,
        basis=basis,
        description="Busing & Levy (1967) four-circle Eulerian diffractometer",
    )


def geometry_sixc() -> AdHocDiffractometer:
    """
    Lohmeier & Vlieg (1993) six-circle surface diffractometer (sixc geometry).

    Basis: xHat=vertical, yHat=longitudinal (beam), zHat=lateral.

    The key structural difference from psic is that both sample and detector
    stacks share a common base stage (alpha, the rotary table), making this
    a coupled geometry rather than the decoupled S4D2 of psic.

    Stack (floor first):
        alpha (shared base): vertical, +x, right-handed  [rotary table]
          --> omega (sample): longitudinal, +y, right-handed
                --> chi:      longitudinal, +y, right-handed
                      --> phi: longitudinal, +y, right-handed
          --> delta (detector): lateral, -z, left-handed
                --> gamma:      vertical,  +x, right-handed

    Reference: M. Lohmeier & E. Vlieg, J. Appl. Cryst. 26, 706-716 (1993).
    """
    stages = [
        Stage("alpha", +XHAT, parent=None, role="sample"),
        Stage("omega", +YHAT, parent="alpha", role="sample"),
        Stage("chi", +YHAT, parent="omega", role="sample"),
        Stage("phi", +YHAT, parent="chi", role="sample"),
        Stage("delta", -ZHAT, parent="alpha", role="detector"),
        Stage("gamma", +XHAT, parent="delta", role="detector"),
    ]
    return AdHocDiffractometer(
        name="sixc",
        stages=stages,
        description=(
            "Lohmeier & Vlieg (1993) six-circle surface diffractometer. "
            "Sample and detector share the alpha (rotary table) base stage."
        ),
    )
