from propeller import Propeller
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState

def eta_test():

    fs = FlightState()
    p = Propeller(fs)
    p.substitutions[p.T] = 100
    p.cost = 1/p.eta
    sol = p.solve()
    print sol.table()

def test():
    "tests"
    eta_test()

if __name__ == "__main__":
    test()

