" vertical tail "
from gpkit import Variable
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import Wing

#pylint: disable=attribute-defined-outside-init, unused-variable

class VerticalTail(Wing):
    "vertical tail model"
    flight_model = TailAero
    sparModel = None
    def setup(self, N=3, lam=0.8):

        self.ascs = Wing.setup(self, N, lam)
        self.planform.substitutions.update({"\\tau": 0.08})
        self.skin.substitutions.update({"\\rho_{CFRP}": 0.049})
        self.foam.substitutions.update({"\\bar{A}_{jh01}": 0.0548,
                                        "\\rho_{foam}": 0.024})

        Vv = Variable("V_v", "-", "vertical tail volume coefficient")
        lv = Variable("l_v", "ft", "vertical tail moment arm")

        return self.ascs, self.components

