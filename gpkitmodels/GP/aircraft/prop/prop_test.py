" propeller tests "
from gpkitmodels.GP.aircraft.prop.propeller import Propeller
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
from gpkit import Model

def simpleprop_test():
    " test simple propeller model "
    fs = FlightState()
    p = Propeller()
    pp = p.flight_model(p, fs)
    m = Model(1/pp.eta, [fs, p, pp])
    m.substitutions.update({"rho": 1.225, "V": 50, "T": 100})
    m.solve()

def test():
    "tests"
    simpleprop_test()

if __name__ == "__main__":
    test()

