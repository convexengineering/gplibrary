" propeller tests "
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp
from gpkitmodels.SP.aircraft.prop.propeller import BladeElementProp

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
    m.solve()

def ME_eta_test():

    fs  = FlightState()
    Propeller.flight_model = BladeElementProp
    p   = Propeller()
    #pp = p.flight_model(p,fs, onDesign = False)
    pp = p.flight_model(p,fs)
    pp.substitutions[pp.T]  = 100
    pp.cost = 1./pp.eta + pp.Q/(1000.*units("N*m")) + p.T_m/(1000*units('N'))
    sol = pp.localsolve(iteration_limit = 400)
    #print sol.table()



class MultiPointProp(Model):
    def setup(self):
        fs1  = FlightState()
        fs2  = FlightState()
        fs2.substitutions.update({"V":25})
        Propeller.flight_model = BladeElementProp
        p   = Propeller()
        #pp = p.flight_model(p,fs, onDesign = False)
        pp1 = p.flight_model(p,fs1, MIL = True)
        pp1.substitutions[pp1.T]  = 100
        pp2 = p.flight_model(p,fs2)
        pp2.substitutions[pp2.T]  = 150

        self.cost = (.5*(1./pp1.eta )#+ pp1.Q/(1000.*units("N*m")))
                    +.5*(1./pp2.eta )#+ pp2.Q/(1000.*units("N*m"))) 
                    + p.T_m/(1000*units('N')) 
                    + p.W/(100*units('lbf')))
    
        return p,pp1, fs1, pp2, fs2

def ME_multiobjectve_test():
    mpp = MultiPointProp()
    sol = mpp.localsolve(iteration_limit = 400)
    print sol.table()

def test():
    "tests"
    simpleprop_test()
    ME_eta_test()
    ME_multiobjectve_test()
if __name__ == "__main__":
    #test()
    mpp = MultiPointProp()
    sol = mpp.localsolve(iteration_limit = 400)
    print sol.table()

