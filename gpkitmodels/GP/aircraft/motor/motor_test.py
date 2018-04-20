from gpkit import Model, parse_variables, SignomialsEnabled, SignomialEquality, units
from motor import Propulsor, Motor, MotorPerf
from gpkitmodels.GP.aircraft.prop.propeller import Propeller, ActuatorProp
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState

class Propulsor_Test(Model):
    """Propulsor Test Model
    """

    def setup(self):
        fs = FlightState()
        p = Propulsor()
        pp = p.flight_model(p,fs)
        pp.substitutions[pp.prop.T] = 100
        self.cost = 1./pp.motor.etam + p.W/(1000*units('lbf')) + 1./pp.prop.eta

        return fs,p,pp

class Actuator_Propulsor_Test(Model):
    """Propulsor Test Model
    """

    def setup(self):
        fs = FlightState()
        Propulsor.prop_flight_model = ActuatorProp
        p = Propulsor()
        pp = p.flight_model(p,fs)
        pp.substitutions[pp.prop.T] = 100
        self.cost = pp.motor.Pelec/(1000*units('W')) + p.W/(1000*units('lbf'))

        return fs,p,pp


def actuator_propulsor_test():

    test = Actuator_Propulsor_Test()
    test.solve()
    
<<<<<<< HEAD


=======
def ME_propulsor_test():
    test = MultiElement_Propulsor_Test()
    sol = test.localsolve()
    print sol.table()
>>>>>>> parent of 46e07bc... updates to tests

def propulsor_test():

    test = Propulsor_Test()
    test.solve()

class Motor_P_Test(Model):
    def setup(self):
        fs = FlightState()
        m  = Motor()
        mp = MotorPerf(m,fs)
        self.mp = mp
        mp.substitutions[m.Qmax] = 100
        mp.substitutions[mp.Q]    = 10
        self.cost = 1./mp.etam + m.W/(100.*units('lbf'))
        return self.mp, fs, m

class speed_280_motor(Model):
    def setup(self):
        fs = FlightState()
        m  = Motor()
        mp = MotorPerf(m,fs)
        self.mp = mp
        mp.substitutions[m.Qmax] = 100
        mp.substitutions[mp.R]   = .7
        mp.substitutions[mp.i0]  = .16
        mp.substitutions[mp.Kv]  = 3800
        mp.substitutions[mp.v]   = 6
        self.cost = 1./mp.etam
        return self.mp, fs
class hacker_q150_45_motor(Model):
    def setup(self):
        fs = FlightState()
        m  = Motor()
        mp = MotorPerf(m,fs)
        self.mp = mp
        mp.substitutions[m.Qmax] = 10000
        mp.substitutions[mp.R]   = .033
        mp.substitutions[mp.i0]  = 4.5
        mp.substitutions[mp.Kv]  = 29
        self.cost = 1./mp.etam
        return self.mp, fs

def motor_test():
    test = Motor_P_Test()
    test.solve()

def test():
    motor_test()
    actuator_propulsor_test()
    propulsor_test()

    
if __name__ == "__main__":
    test()
