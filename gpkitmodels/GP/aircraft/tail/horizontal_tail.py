" horizontal tail "
import numpy as np
from gpkit import Variable
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import Wing
from gpkitmodels.GP.aircraft.wing.wing_core import WingCore

#pylint: disable=attribute-defined-outside-init, unused-variable, no-member

class HorizontalTail(Wing):
    "horizontal tail model"
    flight_model = TailAero
    fillModel = WingCore
    sparModel = None

    def setup(self, N=3):
        self.ascs = Wing.setup(self, N)
        self.planform.substitutions.update(
            {self.planform.AR: 4, self.planform.tau: 0.08,
             self.planform.lam: 0.8})
        self.skin.substitutions.update({self.skin.rhocfrp: 0.049})
        self.foam.substitutions.update({self.foam.Abar: 0.0548,
                                        self.foam.rhocore: 0.024})
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        lh = Variable("l_h", "ft", "horizontal tail moment arm")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")

        return self.ascs, mh*(1+2.0/self.planform["AR"]) <= 2*np.pi

