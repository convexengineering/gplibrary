" propeller tests "
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp


from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
from gpkit import units, Model

def simpleprop_test():
    " test simple propeller model "
    fs = FlightState()
    Propeller.flight_model = ActuatorProp
    p = Propeller()
    pp = p.flight_model(p, fs)
    m = Model(1/pp.eta  + p.W/(100.*units("lbf"))+ pp.Q/(100.*units("N*m")),
              [fs, p, pp])
    m.substitutions.update({"rho": 1.225, "V": 50, "T": 100, "omega":1000})
    sol = m.solve()
    #print sol.table()


def test():
    "tests"
    simpleprop_test()
if __name__ == "__main__":
    test()

