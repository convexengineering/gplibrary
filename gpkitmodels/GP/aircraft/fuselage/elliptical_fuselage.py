" elliptical fuselage.py "
import numpy as np
from gpkit import Variable, Model
from fuel_tank import FuelTank
from fuselage_skin import FuselageSkin

class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def setup(self):

        R = Variable("R", "ft", "fuselage radius")
        l = Variable("l", "ft", "fuselage length")
        S = Variable("S", "ft^2", "wetted fuselage area")
        W = Variable("W", "lbf", "Fuselage weight")
        mfac = Variable("m_{fac}", 2.0, "-", "Fuselage weight margin factor")
        f = Variable("f", "-", "fineness ratio of lenth to diameter")
        k = Variable("k", "-", "fuselage form factor")
        Vol = Variable("\\mathcal{V}", "ft**3", "fuselae volume")
        rhofuel = Variable("\\rho_{fuel}", 6.01, "lbf/gallon",
                           "density of 100LL")
        rhocfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        t = Variable("t", 0.024, "in",
                     "minimum gague fuselage skin thickness")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [
            f == l/R/2,
            k >= 1 + 60/f**3 + f/400,
            3*(S/np.pi)**1.6075 >= 2*(l*R*2)**1.6075 + (2*R)**(2*1.6075),
            Vol <= 4*np.pi/3*(l/2)*R**2,
            W/mfac >= S*rhocfrp*t*g,
            ]

        return constraints

    def flight_model(self, state):
        return FuselageAero(self, state)

class FuselageAero(Model):
    "fuselage drag model"
    def setup(self, static, state):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")
        Cd = Variable("C_d", "-", "fuselage drag coefficient")

        constraints = [
            Re == state["V"]*state["\\rho"]*static["l"]/state["\\mu"],
            Cf >= 0.455/Re**0.3,
            Cd >= Cf*static["k"]
            ]

        return constraints

