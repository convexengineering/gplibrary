" vertical tail "
from gpkit import Variable
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import Wing
from gpkitmodels.GP.aircraft.wing.wing_core import WingCore

#pylint: disable=attribute-defined-outside-init, unused-variable, no-member

class VerticalTail(Wing):
    """vertical tail model

    Upper Unbounded
    ---------------
    lv, Vv, W

    Lower Unbounded
    ---------------
    lv, Vv, b

    """
    flight_model = TailAero
    fillModel = WingCore
    sparModel = None
    def setup(self, N=3):

        self.ascs = Wing.setup(self, N)
        self.planform.substitutions.update(
            {self.planform.tau: 0.08, self.planform.lam: 0.8})
        self.skin.substitutions.update({self.skin.rhocfrp: 0.049})
        self.foam.substitutions.update({self.foam.Abar: 0.0548,
                                        self.foam.rhocore: 0.024})

        Vv = self.Vv = Variable("V_v", "-", "vertical tail volume coefficient")
        lv = self.lv = Variable("l_v", "ft", "vertical tail moment arm")

        return self.ascs, self.components
