" vertical tail "
import numpy as np
from gpkit import Model, Variable
from tail_aero import TailAero
from gpkitmodels.GP.aircraft.wing.wing import AeroSurf
from gpkitmodels.GP.aircraft.wing.constant_taper_chord import c_bar
from gpkitmodels.GP.aircraft.wing.wing_interior import WingInterior
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin

class VerticalTail(Model):
    "vertical tail model"
    def setup(self, N=3, lam=0.8):

        Vv = Variable("V_v", "-", "vertical tail volume coefficient")
        W = Variable("W", "lbf", "vertical tail weight")
        lv = Variable("l_v", "ft", "vertical tail moment arm")
        mfac = Variable("m_{fac}", 1.1, "-", "vertical tail margin factor")

        cb, eta, deta, cbarmac = c_bar(lam, N)
        subdict = {"\\lambda": lam, "\\bar{c}": cb, "\\eta": eta,
                   "\\bar{c}_{ave}": (cb[1:]+cb[:-1])/2, "\\tau": 0.08,
                   "\\bar{c}_{MAC}": cbarmac, "d\\eta": deta, "C_{L_{max}}": 1.5}

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
            self.skin["W"] >= (self.skin["\\rho_{CFRP}"]*self.surf["S"]*2
                               * self.skin["t"]*self.skin["g"]),
            self.foam["W"] >= 2*(
                self.foam["g"]*self.foam["\\rho_{foam}"]
                *self.foam["\\bar{A}_{jh01}"]*self.surf["c_{ave}"]**2
                * (self.surf["b"]/2)*self.surf["d\\eta"]).sum()
            ]

        return constraints, self.surf, self.components

    def flight_model(self, state):
        return TailAero(self, state)
