" horizontal tail "
import numpy as np
from gpkit import Model, Variable
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import AeroSurf
from gpkitmodels.GP.aircraft.wing.constant_taper_chord import c_bar
from gpkitmodels.GP.aircraft.wing.wing_interior import WingInterior
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin

#pylint: disable=invalid-name, too-many-locals, unused-variable
#pylint: disable=attribute-defined-outside-init

class HorizontalTail(AeroSurf):
    "horizontal tail model"
    flight_model = TailAero
    sparModel = None

    def setup(self, N=3, lam=0.8):
        self.ascs = AeroSurf.setup(self, N, lam)
        self.skin.substitutions.update({"\\rho_{CFRP}": 0.049})
        self.foam.substitutions.update({"\\bar{A}_{jh01}": 0.0548,
                                        "\\rho_{foam}": 0.024})
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        lh = Variable("l_h", "ft", "horizontal tail moment arm")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")

        return self.ascs, mh*(1+2.0/self.planform["AR"]) <= 2*np.pi
