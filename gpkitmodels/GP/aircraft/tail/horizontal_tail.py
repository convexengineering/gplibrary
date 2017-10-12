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

class HorizontalTail(Model):
    "horizontal tail model"
    def setup(self, N=3, lam=0.8):

        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        W = Variable("W", "lbf", "horizontal tail weight")
        lh = Variable("l_h", "ft", "horizontal tail moment arm")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")
        mfac = Variable("m_{fac}", 1.1, "-", "horizontal tail margin factor")

        cb, eta, deta, cbarmac = c_bar(lam, N)
        subdict = {"\\lambda": lam, "\\bar{c}": cb, "\\eta": eta, "AR": 5.0,
                   "\\bar{c}_{ave}": (cb[1:]+cb[:-1])/2, "\\tau": 0.08,
                   "\\bar{c}_{MAC}": cbarmac, "d\\eta": deta,
                   "C_{L_{max}}": 1.5}

        self.surf = AeroSurf(N=N)
        self.surf.substitutions.update(subdict)

        self.skin = WingSkin()
        self.skin.substitutions.update({"\\rho_{CFRP}": 0.049})
        self.foam = WingInterior()
        self.foam.substitutions.update({"\\bar{A}_{jh01}": 0.0548})
        self.foam.substitutions.update({"\\rho_{foam}": 0.024})

        self.components = [self.skin, self.foam]

        constraints = [
            W/mfac >= sum([c["W"] for c in self.components]),
            mh*(1+2.0/self.surf["AR"]) <= 2*np.pi,
            self.skin["W"] >= (self.skin["\\rho_{CFRP}"]*self.surf["S"]*2
                               * self.skin["t"]*self.skin["g"]),
            self.foam["W"] >= 2*(
                self.foam["g"]*self.foam["\\rho_{foam}"]
                *self.foam["\\bar{A}_{jh01}"]*self.surf["c_{ave}"]**2
                * (self.surf["b"]/2)*self.surf["d\\eta"]).sum()
            ]

        self.flight_model = TailAero

        return constraints, self.surf, self.components
