" fuselage.py "
import numpy as np
from gpkit import Variable, Model
from fuel_tank import FuelTank
from fuselage_skin import FuselageSkin

class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def setup(self, Wfueltot):

        R = Variable("R", "ft", "fuselage radius")
        l = Variable("l", "ft", "fuselage length")
        S = Variable("S", "ft^2", "fuselage cross sectional area")
        Volavn = Variable("\\mathcal{V}_{avn}", 0.125, "ft^3",
                          "avionics volume")
        Volpay = Variable("\\mathcal{V}_{pay}", 1.0, "ft^3", "payload volume")
        W = Variable("W", "lbf", "Fuselage weight")
        mfac = Variable("m_{fac}", 2.1, "-", "Fuselage weight margin factor")
        lbody = Variable("l_{body}", "ft", "center body length")
        kbody = Variable("k_{body}", "-", "fuselage body fineness ratio")
        knose = Variable("k_{nose}", "-", "fuselage nose finess ratio")
        kbulk = Variable("k_{bulk}", "-", "fuselage bulk finess ratio")
        Swet = Variable("S_{wet}", "ft**2", "fuselage wetted area")
        Sbody = Variable("S_{body}", "ft**2", "wetted surface area of body")
        Snose = Variable("S_{nose}", "ft**2", "wetted surface area of nose")
        Sbulk = Variable("S_{bulk}", "ft**2", "wetted surface area of bulk")
        Volbody = Variable("\\mathcal{V}_{body}", "ft**3", "volume of body")
        Volnose = Variable("\\mathcal{V}_{nose}", "ft**3", "volume of nose")
        Volbulk = Variable("\\mathcal{V}_{bulk}", "ft**3", "volume of bulk")
        Voleng = Variable("\\mathcal{V}_engine", 0.5, "ft**3", "fuslage volume")

        self.fueltank = FuelTank(Wfueltot)
        self.skin = FuselageSkin(Swet, R, lbody)
        self.components = [self.fueltank, self.skin]

        constraints = [
            kbody == lbody/R,
            knose == knose,
            kbulk == kbulk,
            Swet >= Sbody + Snose + Sbulk,
            Sbody >= 2*np.pi*R*lbody,
            Snose**(8./5.) >= (
                (2*np.pi*R**2)**(8./5.)*(1./3. + 2./3.*(knose)**(8./5.))),
            Sbulk >= R**2*(0.012322*kbulk**2 + 1.524925*kbulk + 0.502498),
            Volbody <= np.pi*R**2*lbody,
            Volnose <= 4./6.*np.pi*knose*R**3,
            # Volbulk >= R**2*(0.060402*kbulk**2 + 3.44103*kbulk + 2.5531569),
            Volnose >= Volpay,
            l <= 3*R*(kbody*knose*kbulk)**(1./3),
            S >= np.pi*R**2,
            Volbody >= self.fueltank["\\mathcal{V}"] + Volavn,
            W/mfac >= self.fueltank["W"] + self.skin["W"],
            ]

        return self.components, constraints

    def loading(self, Wcent):
        return FuselageLoading(self, Wcent)

    def flight_model(self, state):
        return FuselageAero(self, state)

class FuselageLoading(Model):
    "fuselage loading cases"
    def setup(self, fuselage, Wcent):

        loading = [fuselage.skin.loading(Wcent)]
        loading.append(fuselage.skin.landing(Wcent))

        return loading

class FuselageAero(Model):
    "fuselage drag model"
    def setup(self, static, state):

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
                * static["k_{nose}"]**1.21682 * static["k_{bulk}"]**0.306251
                + 0.00255095*static["k_{body}"]**-0.0316887
                * static["k_{nose}"]**-0.585489 * static["k_{bulk}"]**1.15394
                + 0.0436011 * static["k_{body}"]**0.0545722
                * static["k_{nose}"]**0.258228 * static["k_{bulk}"]**-1.42664
                + 0.00970479 * static["k_{body}"]**0.8661
                * static["k_{nose}"]**-0.209136 * static["k_{bulk}"]**-0.156166)
            ]

        return constraints

