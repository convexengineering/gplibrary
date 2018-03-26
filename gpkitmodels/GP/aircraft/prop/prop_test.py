" propelle tests "
from propeller import Propeller
#from qprop import QProp
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState

def eta_test():

    fs = FlightState()
    p = Propeller(fs)
    p.substitutions[p.T] = 100
    p.cost = 1/p.eta
    sol = p.solve()
    print sol.table()

def qprop_test():

    fs = FlightState()
    p = QProp(fs)
    p.substitutions[p.T] = 100
    p.cost = 1/p.eta
    sol = p.solve()
    print sol.table()

def test():
    "tests"
    eta_test()
    #qprop_test()

if __name__ == "__main__":
    test()

