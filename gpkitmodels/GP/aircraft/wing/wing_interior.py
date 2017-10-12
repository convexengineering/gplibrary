" wing interior "
from gpkit import Model, Variable

class WingInterior(Model):
    "wing interior model"
    def setup(self):

        W = Variable("W", "lbf", "interior mass of wing")
        dmdy = Variable("(dm/dy)", "kg/m", "foam mass per length")
        rhofoam = Variable("\\rho_{foam}", 0.036, "g/cm^3", "foam density")
        Abar = Variable("\\bar{A}_{jh01}", 0.0753449, "-",
                        "jh01 non dimensional area")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [Abar == Abar]

        return constraints
