" wing skin "
from gpkit import Model, Variable

class WingSkin(Model):
    "wing skin model"
    def setup(self):

        rhocfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        W = Variable("W", "lbf", "wing skin weight")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")
        t = Variable("t", "in", "wing skin thickness")
        tmin = Variable("t_{min}", 0.012, "in", "wing skin min gauge")
        Jtbar = Variable("\\bar{J/t}", 0.01114, "1/mm",
                         "torsional moment of inertia")

        constraints = [t >= tmin]

        self.loading = WingSkinL

        return constraints

class WingSkinL(Model):
    "wing skin loading model for torsional loads in skin"
    def setup(self, static):

        taucfrp = Variable("\\tau_{CFRP}", 570, "MPa", "torsional stress limit")
        Cmw = Variable("C_{m_w}", 0.121, "-", "negative wing moment coefficent")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 45, "m/s", "never exceed vehicle speed")

        constraints = [
            taucfrp >= (1/static["\\bar{J/t}"]/(static["c_{root}"])**2
                        / static.skin["t"]*Cmw*static["S"]*rhosl*Vne**2)]

        return constraints
