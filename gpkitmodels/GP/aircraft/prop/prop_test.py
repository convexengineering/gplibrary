" propeller tests "
from propeller import Actuator_Propeller
from propeller import One_Element_Propeller
from gpkit import units
#from qprop import QProp
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState

def simpleprop_test():
    " test simple propeller model "
    fs = FlightState()
    p = Actuator_Propeller()
    pp = p.flight_model(p, fs)
    m = Model(1/pp.eta+ pp.Q/(100.*units("N*m")), [fs, p, pp])
    m.substitutions.update({"rho": 1.225, "V": 50, "T": 100})
    m.solve()

def OE_eta_test():

    fs = FlightState()
    p = One_Element_Propeller()
    pp = p.flight_model(p,fs)
    pp.substitutions[pp.T] = 1000
    #pp.substitutions[pp.omega] = 500
    pp.substitutions[pp.AR_b] = 12
    pp.cost = 1/pp.eta + pp.Q/(10.*units("N*m"))
    sol = pp.localsolve()
    #sol = pp.debug()
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
    simpleprop_test()
    OE_eta_test()
    #qprop_test()

if __name__ == "__main__":
    test()

