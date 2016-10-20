from gpkit import Variable, Model, units, VectorVariable
from gpkit.constraints.set import ConstraintSet
from gpkit import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
import numpy as np
import matplotlib.pyplot as plt
from solar_irradiance import get_Eirr
from gasmale import SteadyLevelFlight, Aerodynamics, Atmosphere, Fuel
from gasmale import BreguetEndurance, GasMALE
from wind_speeds import get_windspeed

class Power(Model):
    def __init__(self, N, latitude, day, **kwargs):

        # day must be in julian day Jan 1st = 1
        # latitude is in degrees
        esirr, td, tn = get_Eirr(latitude, day)

        ESirr = Variable("(E/S)_{irr}", esirr, "W*hr/m^2",
                         "Average daytime solar energy")
        Poper = VectorVariable(N, "P_{oper}", "W", "Aircraft operating power")
        Pacc = Variable("P_{acc}", 0.0, "W", "Accessory power draw")
        eta_solar = Variable("\\eta_{solar}", 0.2, "-",
                             "Solar cell efficiency")
        eta_charge = Variable("\\eta_{charge}", 0.98, "-",
                              "Battery charging efficiency")
        eta_discharge = Variable("\\eta_{discharge}", 0.98, "-",
                                 "Battery discharging efficiency")
        tday = Variable("t_{day}", td, "hr", "Daylight span")
        tnight = Variable("t_{night}", tn, "hr", "Night span")
        Ssolar = Variable("S_{solar}", "ft**2", "solar cell area")
        Ebatt = Variable("E_{batt}", "J", "total battery energy")
        P_shaft = VectorVariable(N, "P_{shaft}", "hp", "Shaft power")
        S = Variable("S", "ft**2", "wing area")

        constraints = [
            ESirr*eta_solar*Ssolar >= Poper*tday + Ebatt/eta_charge,
            Ssolar <= S,
            Poper >= P_shaft + Pacc,
            Ebatt >= Poper*tnight/eta_discharge]
        Model.__init__(self, None, constraints, **kwargs)

class SolarSimple(Model):
    def __init__(self, latitude=45, avail=90, day=355, altitude=16000, N=1, **kwargs):
        # http://sky-sailor.ethz.ch/docs/Conceptual_Design_of_Solar_Powered_Airplanes_for_continuous_flight2.pdf

        slf = SteadyLevelFlight(N, [0.8]*N)
        aero = Aerodynamics(N, jh01=False)
        atm = Atmosphere(N, [altitude]*N)
        power = Power(N, latitude, day)
        aero.substitutions.update({"AR": 25})
        wind = get_windspeed(latitude, avail, altitude)

        fstructures = Variable("f_{structures}", 0.35, "-",
                               "fractional structural weight")
        Wstructures = Variable("W_{structures}", "lbf", "structural weight")
        mtow = Variable("MTOW", "lbf", "max take off weight")
        Vmin = Variable("V_{min}", wind, "m/s", "minimum velocity")
        Wpay = Variable("W_{pay}", 1, "lbf", "payload")
        Wsolar = Variable("W_{solar}", "lbf", "solar cell weight")
        rhosolar = Variable("\\rho_{solar}", 0.3, "kg/m^2",
                            "solar cell area density")
        g = Variable("g", 9.81, "m/s**2", "gravitational constant")
        hbatt = Variable("h_{batt}", 350, "W*hr/kg", "battery energy density")
        Wbatt = Variable("W_{batt}", "lbf", "battery weight")

        constraints = [mtow >= Wstructures + Wpay + Wsolar + Wbatt,
                       Wstructures >= mtow*fstructures,
                       Wsolar >= rhosolar*power["S_{solar}"]*g,
                       Wbatt >= power["E_{batt}"]/hbatt*g,
                       slf["V"] >= Vmin,
                       slf["W_{N+1}"] == mtow,
                       slf["W_{N}"] == mtow,
                      ]

        cost = mtow

        lc = LinkedConstraintSet([constraints, slf, aero, atm, power])

        Model.__init__(self, cost, lc, **kwargs)

if __name__ == "__main__":
    M = SolarSimple(latitude=35, avail=80, altitude=70000)
    sol = M.solve("mosek")
