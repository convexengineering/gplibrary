" fuselage.py "
import numpy as np
from gpkit import Variable, Model
from fuel_tank import FuelTank
from fuselage_skin import FuselageSkin

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
        phi = Variable("\\phi", 12, "-", "fuselage fineness ratio")

        self.fueltank = FuelTank(Wfueltot)
        self.skin = FuselageSkin(S, d, l)
        self.components = [self.fueltank, self.skin]
        self.flight_model = FuselageAero
        self.loading = FuselageLoading

        constraints = [
            phi == l/(d/2),
            S >= np.pi*d*l + np.pi*d**2,
            np.pi*(d/2)**2*l >= self.fueltank["\\mathcal{V}"] + Volavn,
            d >= hengine,
            W/mfac >= self.fueltank["W"] + self.skin["W"],
            ]

        Model.__init__(self, None, [self.components, constraints], **kwargs)

class FuselageLoading(Model):
    "fuselage loading cases"
    def __init__(self, fuselage, Wcent):

        loading = [fuselage.skin.loading(fuselage.skin, Wcent)]
        loading.append(fuselage.skin.landing(fuselage.skin, Wcent))

        Model.__init__(self, None, loading)

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

