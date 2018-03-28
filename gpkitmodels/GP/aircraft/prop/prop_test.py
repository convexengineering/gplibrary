" propelle tests "
from propeller import Propeller
from gpkit import units
#from qprop import QProp
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState

def eta_test():

    fs = FlightState()
    p = Propeller()
    pp = p.flight_model(fs)
    pp.substitutions[pp.T] = 100
    pp.cost = 1/pp.eta + pp.Q/(100.*units("N*m"))
    sol = pp.solve()
    print sol.table()

#def qprop_test():
#
#    fs = FlightState()
#    p = QProp(fs)
#    p.substitutions[p.T] = 100
#    p.cost = 1/p.eta
#    sol = p.solve()
#    print sol.table()

def test():
    "tests"
    eta_test()
    #qprop_test()

if __name__ == "__main__":
    test()

