" wing skin "
from gpkit import Model, Variable

class WingSkin(Model):
    "wing skin model"
    def __init__(self, S, croot, b, **kwargs):

        rhocfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        W = Variable("W", "lbf", "wing skin weight")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")
        t = Variable("t", "in", "wing skin thickness")
        tmin = Variable("t_{min}", 0.012, "in",
                        "minimum gague wing skin thickness")
        Jtbar = Variable("\\bar{J/t}", 0.01114, "1/mm",
                         "torsional moment of inertia")

        self.loading = WingSkinL

        constraints = [W >= rhocfrp*S*2*t*g,
                       t >= tmin,
                       Jtbar == Jtbar,
                       b == b,
                       croot == croot]

        Model.__init__(self, None, constraints, **kwargs)

class WingSkinL(Model):
    "wing skin loading model for torsional loads in skin"
    def __init__(self, static, **kwargs):

        taucfrp = Variable("\\tau_{CFRP}", 570, "MPa", "torsional stress limit")
        Cmw = Variable("C_{m_w}", 0.121, "-", "negative wing moment coefficent")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 45, "m/s", "never exceed vehicle speed")

        constraints = [
            taucfrp >= (1/static["\\bar{J/t}"]/(static["c_{root}"])**2
                        / static["t"]*Cmw*static["S"]*rhosl*Vne**2)]

        Model.__init__(self, None, constraints, **kwargs)
