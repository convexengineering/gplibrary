" vertical tail "
from gpkit import parse_variables
from .tail_aero import TailAero
from gplibrary.GP.aircraft.wing.wing import Wing
from gplibrary.GP.aircraft.wing.wing_core import WingCore
from gplibrary.GP.aircraft.wing.wing_skin import WingSkin

#pylint: disable=attribute-defined-outside-init, no-member, exec-used

class VerticalTail(Wing):
    """ Vertical Tail Model

    Variables
    ---------
    Vv                  [-]         vertical tail volume coefficient
    lv                  [ft]        vertical tail moment arm

    Upper Unbounded
    ---------------
    lv, Vv, W, planform.tau (if not sparModel)

    Lower Unbounded
    ---------------
    lv, Vv, planform.b, planform.tau (if not sparModel)
    spar.Sy (if sparModel), spar.J (if sparJ)

    LaTex Strings
    -------------
    Vv      V_{\\mathrm{v}}
    lv      l_{\\mathrm{v}}

    """

    flight_model = TailAero
    fillModel = WingCore
    sparModel = None

    def setup(self, N=3):
        exec parse_variables(VerticalTail.__doc__)

        self.ascs = Wing.setup(self, N)
        self.planform.substitutions.update(
            {self.planform.lam: 0.8, self.planform.AR: 4})
        if self.fillModel:
            self.foam.substitutions.update({self.foam.Abar: 0.0548,
                                            self.foam.material.rho: 0.024})

        return self.ascs
