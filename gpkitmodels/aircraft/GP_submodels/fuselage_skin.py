" fuselage skin "
import numpy as np
from gpkit import Model, Variable

class FuselageSkin(Model):
    "fuselage skin model"
    def __init__(self, S, d, l):

        W = Variable("W", "lbf", "fuselage skin weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        rhokevlar = Variable("\\rho_{kevlar}", 1.3629, "g/cm**3",
                             "kevlar density")
        t = Variable("t", "in", "skin thickness")
        tmin = Variable("t_{min}", 0.03, "in", "minimum skin thickness")
        I = Variable("I", "m**4", "wing skin moment of inertia")

        self.loading = FuselageSkinL

        constraints = [W >= S*rhokevlar*t*g,
                       t >= tmin,
                       I <= np.pi*(d/2)**3*t,
                       l == l]

        Model.__init__(self, None, constraints)

class FuselageSkinL(Model):
    "fuselage skin loading"
    def __init__(self, static, Wcent):

        Mh = Variable("M_h", "N*m", "horizontal axis center fuselage moment")
        Nmax = Variable("N_{max}", 5, "-", "max loading")
        sigmakevlar = Variable("\\sigma_{Kevlar}", 190, "MPa",
                               "stress strength of Kevlar")

        constraints = [Mh >= Nmax*Wcent/4*static["l"],
                       sigmakevlar >= Mh*static["d"]/2/static["I"]]

        Model.__init__(self, None, constraints)
