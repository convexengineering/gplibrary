" vertical tail "
from gpkit import Model, Variable
from tail_aero import TailAero

class VerticalTail(Model):
    "vertical tail model"
    def setup(self, lam=0.7):

        W = Variable("W", "lbf", "one vertical tail weight")
        Sv = Variable("S", "ft**2", "total vertical tail surface area")
        Vv = Variable("V_v", 0.05, "-", "vertical tail volume coefficient")
        ARv = Variable("AR_v", "-", "vertical tail aspect ratio")
        bv = Variable("b", "ft", "one vertical tail span")
        rhofoam = Variable("\\rho_{foam}", 1.5, "lbf/ft^3",
                           "Density of formular 250")
        rhoskin = Variable("\\rho_{skin}", 0.1, "g/cm**2",
                           "vertical tail skin density")
        Abar = Variable("\\bar{A}_{NACA0008}", 0.0548, "-",
                        "cross sectional area of NACA 0008")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        lv = Variable("l_v", "ft", "horizontal tail moment arm")
        ctv = Variable("c_{t_v}", "ft", "vertical tail tip chord")
        crv = Variable("c_{r_v}", "ft", "vertical tail root chord")
        lamv = Variable("\\lambda", lam, "-", "vertical tail taper ratio")
        lamvfac = Variable("\\lambda_v/(\\lambda_v+1)", lam/(lam+1), "-",
                           "vertical tail taper ratio factor")
        CLvtmax = Variable("C_{L_{max}}", 1.1, "-",
                           "maximum CL of vertical tail")
        mfac = Variable("m_{fac}", 1.1, "-", "vertical tail margin factor")
        tau = Variable("\\tau", 0.08, "-", "vertical tail thickness ratio")

        constraints = [bv**2 == ARv*Sv,
                       W/mfac >= rhofoam*Sv**2/bv*Abar + g*rhoskin*Sv,
                       ctv == 2*Sv/bv*lamvfac,
                       crv == ctv/lam,
                       lamv == lamv
                      ]

        return constraints

    def flight_model(self, state):
        return TailAero(self, state)
