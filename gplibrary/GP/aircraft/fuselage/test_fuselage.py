from gpkit import Model
from gplibrary.GP.aircraft.fuselage.elliptical_fuselage import Fuselage
from gplibrary.GP.aircraft.wing.wing_test import FlightState

def test_ellp():
    "elliptical fuselage test"
    f = Fuselage()
    fs = FlightState()
    faero = f.flight_model(f, fs)
    f.substitutions[f.Vol] = 1.33

    m = Model(f.W*faero.Cd, [f, fs, faero])
    m.solve()

def test():
    "tests"
    test_ellp()

if __name__ == "__main__":
    test()

