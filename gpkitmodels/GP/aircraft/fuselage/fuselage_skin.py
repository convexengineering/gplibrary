" fuselage skin "
import numpy as np
from gpkit import Model, Variable

class FuselageSkin(Model):
    "fuselage skin model"
    def setup(self, S, R, l):

        W = Variable("W", "lbf", "fuselage skin weight")
        m = Variable("m", "kg", "fuselage skin mass")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        rhokevlar = Variable("\\rho_{kevlar}", 1.3629, "g/cm**3",
                             "kevlar density")
        t = Variable("t", "in", "skin thickness")
        tmin = Variable("t_{min}", 0.03, "in", "minimum skin thickness")
        I = Variable("I", "m**4", "wing skin moment of inertia")
        #Ig = Variable("I_G", "kg*m**2", "mass moment of inertia")
        E = Variable("E", 30, "GPa", "Young's Modulus of Kevlar")

        self.l = l

        constraints = [m >= S*rhokevlar*t,
                       W >= m*g,
                       t >= tmin,
                       I <= np.pi*R**3*t,
                       #Ig >= m*(4*R**2 + 4*R*t + t**2),
                       l == l,
                       E == E]

        return constraints

    def loading(self, Wcent):
        return FuselageSkinL(self, Wcent)

    def landing(self, Wcent):
        return FuselageLanding(self, Wcent)

class FuselageSkinL(Model):
    "fuselage skin loading"
    def setup(self, static, Wcent):

        Mh = Variable("M_h", "N*m", "horizontal axis center fuselage moment")
        Nmax = Variable("N_{max}", 5, "-", "max loading")
        sigmakevlar = Variable("\\sigma_{Kevlar}", 190, "MPa",
                               "stress strength of Kevlar")
        q = Variable("q", "N/m", "distributed load")
        kappa = Variable("\\kappa", 0.05, "-", "maximum tip deflection ratio")



        constraints = [
            Mh >= Nmax*Wcent/4*static.l,
            sigmakevlar >= Mh*static["R"]/static["I"],
            q >= Wcent*Nmax/static.l,
            kappa*static.l/2 >= (
                q*(static.l/2)**4/(8*static["E"]*static["I"]))
            ]

        return constraints

class FuselageLanding(Model):
    "fuselage loading case"
    def setup(self, static, Wcent):

        F = Variable("F", "lbf", "maximum landing force")
        Nmax = Variable("N_{max}", 5, "-", "maximum landing load factor")
        a = Variable("a", "m/s**2", "landing vertical acceleration")
        omegadot = Variable("\\dot{\\omega}", "1/s**2",
                            "angular acceleration about rear fuselage")
        Mg = Variable("M_G", "N*m", "landing moment about center of mass")
        sigmakevlar = Variable("\\sigma_{Kevlar}", 190, "MPa",
                               "stress strength of Kevlar")

        constraints = [F >= Wcent*Nmax,
                       a >= F/static["m"],
                       omegadot >= a/(static["l_{body}"]/2),
                       Mg >= static["I_G"]*omegadot,
                       sigmakevlar >= Mg*static["R"]/static["I"]
                      ]

        return constraints
