" fuselage.py "
import numpy as np
from gpkit import Variable, Model

class FuelTank(Model):
    """
    Returns the weight of the fuel tank.  Assumes a cylinder shape with some
    fineness ratio
    """
    def __init__(self, Wfueltot, **kwargs):

        W = Variable("W", "lbf", "fuel tank weight")
        f = Variable("f", 0.03, "-", "fraction fuel tank weight to fuel weight")
        mfac = Variable("m_{fac}", 1.1, "-", "fuel volume margin factor")
        rhofuel = Variable("\\rho_{fuel}", 6.01, "lbf/gallon",
                           "density of 100LL")
        Vol = Variable("\\mathcal{V}", "ft^3", "fuel tank volume")

        constraints = [W >= f*Wfueltot,
                       Vol/mfac >= Wfueltot/rhofuel,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def __init__(self, Wfueltot, **kwargs):

        d = Variable("d", "ft", "fuselage diameter")
        l = Variable("l", "ft", "fuselage length")
        S = Variable("S", "ft^2", "Fuselage surface area")
        Volavn = Variable("\\mathcal{V}_{avn}", 0.125, "ft^3",
                          "Avionics volume")
        W = Variable("W", "lbf", "Fuselage weight")
        mfac = Variable("m_{fac}", 2.1, "-", "Fuselage weight margin factor")
        hengine = Variable("h_{engine}", 6, "in", "engine height")
        phi = Variable("\\phi", 6, "-", "fuselage fineness ratio")

        self.fueltank = FuelTank(Wfueltot)
        self.skin = FuselageSkin(S, d, l)
        self.components = [self.fueltank, self.skin]
        self.flight_model = FuselageAero
        self.loading = FuselageLoading

        constraints = [
            phi == l/d,
            S >= np.pi*d*l + np.pi*d**2,
            np.pi*(d/2)**2*l >= self.fueltank["\\mathcal{V}"] + Volavn,
            d >= hengine,
            W/mfac >= self.fueltank["W"] + self.skin["W"],
            ]

        Model.__init__(self, None, [self.components, constraints], **kwargs)

class FuselageLoading(Model):
    "fuselage loading cases"
    def __init__(self, fuselage, Wcent):

        skinloading = fuselage.skin.loading(fuselage.skin, Wcent)

        Model.__init__(self, None, skinloading)

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

class FuselageAero(Model):
    "fuselage drag model"
    def __init__(self, static, state, **kwargs):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == state["V"]*state["\\rho"]*static["l"]/state["\\mu"],
            Cf >= 0.455/Re**0.3
            ]

        Model.__init__(self, None, constraints, **kwargs)

