" horizontal tail "
import numpy as np
from gpkit import Model, Variable

class HorizontalTail(Model):
    "horizontal tail model"
    def __init__(self, lam=0.8, **kwargs):
        Sh = Variable("S", "ft**2", "horizontal tail area")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        ARh = Variable("AR_h", "-", "horizontal tail aspect ratio")
        Abar = Variable("\\bar{A}_{NACA0008}", 0.0548, "-",
                        "cross sectional area of NACA 0008")
        rhofoam = Variable("\\rho_{foam}", 1.5, "lbf/ft^3",
                           "Density of formular 250")
        rhoskin = Variable("\\rho_{skin}", 0.1, "g/cm**2",
                           "horizontal tail skin density")
        bh = Variable("b_h", "ft", "horizontal tail span")
        W = Variable("W", "lbf", "horizontal tail weight")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        lh = Variable("l_h", "ft", "horizontal tail moment arm")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")
        cth = Variable("c_{t_h}", "ft", "horizontal tail tip chord")
        crh = Variable("c_{r_h}", "ft", "horizontal tail root chord")
        lamhfac = Variable("\\lambda_h/(\\lambda_h+1)", lam/(lam+1), "-",
                           "horizontal tail taper ratio factor")
        CLhtmax = Variable("C_{L_{max}}", "-", "maximum CL of horizontal tail")
        mfac = Variable("m_{fac}", 1.1, "-", "horizontal tail margin factor")
        tau = Variable("\\tau", 0.08, "-", "horizontal tail thickness ratio")

        self.flight_model = HorizontalTailAero

        constraints = [
            bh**2 == ARh*Sh,
            mh*(1+2/ARh) <= 2*np.pi,
            W/mfac >= g*rhoskin*Sh + rhofoam*Sh**2/bh*Abar,
            cth == 2*Sh/bh*lamhfac,
            crh == cth/lam,
            lh == lh,
            CLhmin == CLhmin,
            CLhtmax == CLhtmax,
            Vh == Vh,
            tau == tau,
            ]

        Model.__init__(self, None, constraints, **kwargs)

class HorizontalTailAero(Model):
    "horizontal tail aero model"
    def __init__(self, static, state, **kwargs):

        Re = Variable("Re", "-", "horizontal tail reynolds number")
        Cd = Variable("C_d", "-", "horizontal tail drag coefficient")

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["S"]/static["b_h"]
                   / state["\\mu"]),
            Cd**70.5599 >= (7.42688e-90*(Re/1000)**-33.0637
                            * (static["\\tau"]*100)**18.0419
                            + 5.02826e-163*(Re/1000)**-18.7959
                            * (static["\\tau"]*100)**53.1879
                            + 4.22901e-77*(Re/1000)**-41.1704
                            * (static["\\tau"]*100)**28.4609)
            ]

        Model.__init__(self, None, constraints, **kwargs)
