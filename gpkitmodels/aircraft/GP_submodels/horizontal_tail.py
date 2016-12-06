" horizontal tail "
import numpy as np
from gpkit import Model, Variable
from tail_aero import TailAero

class HorizontalTail(Model):
    "horizontal tail model"
    def setup(self, lam=0.8):
        Sh = Variable("S", "ft**2", "horizontal tail area")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        ARh = Variable("AR_h", "-", "horizontal tail aspect ratio")
        Abar = Variable("\\bar{A}_{NACA0008}", 0.0548, "-",
                        "cross sectional area of NACA 0008")
        rhofoam = Variable("\\rho_{foam}", 1.5, "lbf/ft^3",
                           "Density of formular 250")
        rhoskin = Variable("\\rho_{skin}", 0.1, "g/cm**2",
                           "horizontal tail skin density")
        bh = Variable("b", "ft", "horizontal tail span")
        W = Variable("W", "lbf", "horizontal tail weight")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        lh = Variable("l_h", "ft", "horizontal tail moment arm")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")
        cth = Variable("c_{t_h}", "ft", "horizontal tail tip chord")
        crh = Variable("c_{r_h}", "ft", "horizontal tail root chord")
        lamh = Variable("\\lambda", lam, "-", "horizontal tail taper ratio")
        lamhfac = Variable("\\lambda_h/(\\lambda_h+1)", lam/(lam+1), "-",
                           "horizontal tail taper ratio factor")
        CLhtmax = Variable("C_{L_{max}}", "-", "maximum CL of horizontal tail")
        mfac = Variable("m_{fac}", 1.1, "-", "horizontal tail margin factor")
        tau = Variable("\\tau", 0.08, "-", "horizontal tail thickness ratio")

        constraints = [
            bh**2 == ARh*Sh,
            mh*(1+2/ARh) <= 2*np.pi,
            W/mfac >= g*rhoskin*Sh + rhofoam*Sh**2/bh*Abar,
            cth == 2*Sh/bh*lamhfac,
            crh == cth/lam,
            lamh == lamh
            ]

        return constraints

    def flight_model(self, state):
        return TailAero(self, state)
