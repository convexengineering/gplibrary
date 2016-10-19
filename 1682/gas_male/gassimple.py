import numpy as np
from gpkit import Model, Variable, LinkedConstraintSet
from gasmale import SteadyLevelFlight, Aerodynamics, Atmosphere, Fuel
from gasmale import BreguetEndurance, GasMALE
from plotting import plot_sweep
from gen_tex import gen_tex_fig
import matplotlib.pyplot as plt
from wind_speeds import get_windspeed

class GasSimple(Model):
    def __init__(self, latitude=45, altitude=16000, N=5, avail=90, **kwargs):

        slf = SteadyLevelFlight(N, [0.8]*N)
        fuel = Fuel(N)
        aero = Aerodynamics(N)
        atm = Atmosphere(N, [altitude]*5)
        be = BreguetEndurance(N)
        aero.substitutions.update({"AR": 25})
        be.substitutions.update({"BSFC": np.linspace(0.6, 0.6, N)})
        wind = get_windspeed(latitude, avail, altitude)

        etaengine = Variable("\\eta_{engine}", 0.95, "-", "motor efficiency")
        fstructures = Variable("f_{structures}", 0.35, "-",
                               "fractional structural weight")
        Wstructures = Variable("W_{structures}", "lbf", "structural weight")
        mtow = Variable("MTOW", 200, "lbf", "max take off weight")
        Wfueltot = Variable("W_{fuel-tot}", "lbf", "total fuel weight")
        tflight = Variable("t_{flight}", "days", "flight time")
        Vmin = Variable("V_{min}", wind, "m/s", "minimum velocity")
        Wpay = Variable("W_{pay}", 1, "lbf", "payload")

        constraints = [mtow >= fuel["W_{start}"],
                       mtow >= Wstructures + Wfueltot + Wpay,
                       Wstructures >= mtow*fstructures,
                       Wstructures + Wpay <= fuel["W_{end}"],
                       Wfueltot >= fuel["W_{fuel-fs}"],
                       slf["P_{shaft}"] == be["P_{shaft-tot}"]*etaengine,
                       be["t"] >= tflight/N,
                       slf["V"] >= Vmin,
                      ]

        cost = 1/tflight

        lc = LinkedConstraintSet([constraints, slf, fuel, aero, atm, be])

        Model.__init__(self, cost, lc, **kwargs)

if __name__ == "__main__":
    M = GasSimple()
    sol = M.solve("mosek")
