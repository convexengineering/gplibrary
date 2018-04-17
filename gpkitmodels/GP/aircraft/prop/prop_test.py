" propeller tests "
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp
from gpkitmodels.GP.aircraft.prop.propeller import SimpleQProp, Multi_Element_Propeller

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

def OE_eta_test():
    " simple qprop test "

    fs = FlightState()
    Propeller.flight_model = SimpleQProp
    p = Propeller()
    pp = p.flight_model(p, fs)
    pp.substitutions[pp.T] = 100
    pp.substitutions[pp.AR_b] = 12
    m = Model(1/pp.eta + pp.Q/(10.*units("N*m"))+ p.W/(100.*units("lbf")),
              [pp, p])
    #m.debug()
    sol = m.localsolve()
    #print sol.table()

def ME_eta_test():

    fs  = FlightState()
    p   = Multi_Element_Propeller()
    p.substitutions[p.T_m]  = 100
    #p.substitutions[p.R]    = 2
    pp = p.flight_model(p,fs)
    pp.substitutions[pp.T]  = 100
   
    #pp.substitutions[pp.omega] = 500
    
    #pp.substitutions[pp.omega] = 1000
    pp.cost = 1./pp.eta + pp.Q/(100.*units("N*m"))
    #sol = pp.debug()
    sol = pp.localsolve(iteration_limit = 400)
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
    OE_eta_test()
    simpleprop_test()
    ME_eta_test()
if __name__ == "__main__":
    test()

