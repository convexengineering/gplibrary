" wing interior "
from gpkit import Model, Variable

class WingInterior(Model):
    "wing interior model"
    def setup(self, surface):

        W = Variable("W", "lbf", "interior mass of wing")
        rhofoam = Variable("\\rho_{foam}", 0.036, "g/cm^3", "foam density")
        Abar = Variable("\\bar{A}_{jh01}", 0.0753449, "-",
                        "jh01 non dimensional area")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [W >= 2*(g*rhofoam*Abar*surface.cave**2
                               * surface.b/2*surface.deta).sum()]

        return constraints
