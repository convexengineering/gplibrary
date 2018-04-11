" propeller tests "
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp
from gpkitmodels.GP.aircraft.prop.propeller import SimpleQProp
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState
from gpkit import units, Model

def simpleprop_test():
    " test simple propeller model "
    fs = FlightState()
    Propeller.flight_model = ActuatorProp
    p = Propeller()
    pp = p.flight_model(p, fs)
    m = Model(1/pp.eta + pp.Q/(100.*units("N*m")) + p.W/(100.*units("lbf")),
              [fs, p, pp])
    m.substitutions.update({"rho": 1.225, "V": 50, "T": 100, "T_m": 40})
    m.solve()

def OE_eta_test():
    " simple qprop test "

    fs = FlightState()
    Propeller.flight_model = SimpleQProp
    p = Propeller()
    pp = p.flight_model(p, fs)
    pp.substitutions[pp.T] = 1000
    pp.substitutions[pp.AR_b] = 12
    m = Model(1/pp.eta + pp.Q/(10.*units("N*m")) + p.W/(100.*units("lbf")),
              [pp, p])
    m.localsolve()


def test():
    "tests"
    simpleprop_test()
    OE_eta_test()

if __name__ == "__main__":
    test()

