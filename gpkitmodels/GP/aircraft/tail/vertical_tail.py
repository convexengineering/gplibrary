" vertical tail "
from gpkit import Model, Variable
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import AeroSurf
from gpkitmodels.GP.aircraft.wing.constant_taper_chord import c_bar
from gpkitmodels.GP.aircraft.wing.wing_interior import WingInterior
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin

#pylint: disable=invalid-name, too-many-locals, unused-variable
#pylint: disable=attribute-defined-outside-init

class VerticalTail(AeroSurf):
    "vertical tail model"
    flight_model = TailAero
    sparModel = None
    def setup(self, N=3, lam=0.8):

        self.ascs = AeroSurf.setup(self, N, lam)
        self.skin.substitutions.update({"\\rho_{CFRP}": 0.049})
        self.foam.substitutions.update({"\\bar{A}_{jh01}": 0.0548,
                                        "\\rho_{foam}": 0.024})

        Vv = Variable("V_v", "-", "vertical tail volume coefficient")
        lv = Variable("l_v", "ft", "vertical tail moment arm")

        return self.ascs, self.components
