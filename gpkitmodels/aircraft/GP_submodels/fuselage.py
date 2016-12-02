" fuselage.py "
import numpy as np
from gpkit import Variable, Model
from fuel_tank import FuelTank
from fuselage_skin import FuselageSkin

class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def setup(self, Wfueltot):

        d = Variable("d", "ft", "fuselage diameter")
        l = Variable("l", "ft", "fuselage length")
        S = Variable("S", "ft^2", "fuselage cross sectional area")
        Volavn = Variable("\\mathcal{V}_{avn}", 0.125, "ft^3",
                          "Avionics volume")
        W = Variable("W", "lbf", "Fuselage weight")
        mfac = Variable("m_{fac}", 2.1, "-", "Fuselage weight margin factor")
        hengine = Variable("h_{engine}", 12, "in", "engine height")
        kbody = Variable("k_{body}", "-", "fuselage body fineness ratio")
        knose = Variable("k_{nose}", 2, "-", "fuselage nose finess ratio")
        ktail = Variable("k_{tail}", 4, "-", "fuselage tail finess ratio")
        Swet = Variable("S_{wet}", "ft**2", "fuselage wetted area")

        self.fueltank = FuelTank(Wfueltot)
        self.skin = FuselageSkin(Swet, d, l)
        self.components = [self.fueltank, self.skin]

        constraints = [
            kbody == l/(d/2),
            knose == knose,
            ktail == ktail,
            Swet >= np.pi*d*l*1.5,
            S >= np.pi*(d/2)**2,
            np.pi*(d/2)**2*l >= self.fueltank["\\mathcal{V}"] + Volavn,
            d >= hengine,
            W/mfac >= self.fueltank["W"] + self.skin["W"],
            ]

        return self.components, constraints

    def loading(self, Wcent):
        return FuselageLoading(self, Wcent)

    def flight_model(self, state):
        return FuselageAero(self, state)

class FuselageLoading(Model):
    "fuselage loading cases"
    def __init__(self, fuselage, Wcent):

        loading = [fuselage.skin.loading(Wcent)]
        loading.append(fuselage.skin.landing(Wcent))

        Model.__init__(self, None, loading)

class FuselageAero(Model):
    "fuselage drag model"
    def __init__(self, static, state, **kwargs):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")
        Reref = Variable("Re_{ref}", 1e6, "-", "reference Reynolds number")
        Cfref = Variable("C_{r_{ref}}", "-",
                         "reference skin friction coefficient")
        Cd = Variable("C_d", "-", "fuselage drag coefficient")

        constraints = [
            Re == state["V"]*state["\\rho"]*static["l"]/state["\\mu"],
            Cf >= 0.455/Re**0.3,
            Cfref == 0.455/Reref**0.3,
            Cd**0.996232 >= Cf/Cfref*(
                0.00243049*static["k_{body}"]**0.033607
                * static["k_{nose}"]**1.21682 * static["k_{tail}"]**0.306251
                + 0.00255095*static["k_{body}"]**-0.0316887
                * static["k_{nose}"]**-0.585489 * static["k_{tail}"]**1.15394
                + 0.0436011 * static["k_{body}"]**0.0545722
                * static["k_{nose}"]**0.258228 * static["k_{tail}"]**-1.42664
                + 0.00970479 * static["k_{body}"]**0.8661
                * static["k_{nose}"]**-0.209136 * static["k_{tail}"]**-0.156166)
            ]

        Model.__init__(self, None, constraints, **kwargs)

