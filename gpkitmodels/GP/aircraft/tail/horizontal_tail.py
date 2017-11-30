" horizontal tail "
import numpy as np
from gpkit import parse_variables
from .tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import Wing
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin
from gpkitmodels.GP.aircraft.wing.wing_core import WingCore
from gpkitmodels.GP.materials.kevlar import Kevlar

#pylint: disable=attribute-defined-outside-init, no-member
#pylint: disable=exec-used, undefined-variable

class HorizontalTail(Wing):
    """ Horizontal Tail Model

    Variables
    ---------
    Vh                          [-]     horizontal tail volume coefficient
    lh                          [ft]    horizontal tail moment arm
    CLhmin              0.75    [-]     max downlift coefficient
    mh                          [-]     horizontal tail span effectiveness

    Upper Unbounded
    ---------------
    lh, Vh, W

    Lower Unbounded
    ---------------
    lh, Vh, b, mh

    LaTex Strings
    -------------
    Vh          V_{\\mathrm{h}}
    lh          l_{\\mathrm{h}}
    CLmin       C_{L_{\\mathrm{min}}}
    mh          m_{\\mathrm{h}}

    """
    flight_model = TailAero
    fillModel = WingCore
    sparModel = None

    def setup(self, N=3):
        exec parse_variables(HorizontalTail.__doc__)

        WingSkin.material = Kevlar()
        self.ascs = Wing.setup(self, N)
        self.planform.substitutions.update(
            {self.planform.AR: 4, self.planform.tau: 0.08,
             self.planform.lam: 0.8})
        if self.fillModel:
            self.foam.substitutions.update({self.foam.Abar: 0.0548,
                                            self.foam.material.rho: 0.024})

        return self.ascs, mh*(1+2.0/self.planform["AR"]) <= 2*np.pi
