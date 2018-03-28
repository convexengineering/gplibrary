from gpkit import Model, parse_variables, SignomialsEnabled, SignomialEquality, units
from motor import Propulsor, ElecMotor, ElecMotor_Performance
from gpkitmodels.GP.aircraft.wing.wing_test import FlightState

class Propulsor_Test(Model):
    """Propulsor Test Model
    """

    def setup(self):
        fs = FlightState()
        p = Propulsor()
        pp = p.flight_model(fs)
        pp.substitutions[pp.prop.T] = 50
        self.cost = 1./pp.motor.etam + p.W/(10000000*units('lbf')) + 1./pp.prop.eta

        return fs,p,pp

def propulsor_test():

    test = Propulsor_Test()
    #sol = test.debug()
    sol = test.solve()
    #print sol.table()

class Motor_P_Test(Model):
    def setup(self):
        fs = FlightState()
        m  = ElecMotor()
        mp = ElecMotor_Performance(m,fs)
        self.mp = mp
        mp.substitutions[m.Qmax] = 100
        mp.substitutions[mp.Q]    = 10
        self.cost = 1./mp.etam
        return self.mp, fs

def motor_test():
    test = Motor_P_Test()
    sol = test.solve()
    #sol = test.debug()

    #print sol.table()
    
def motor_eta_speed():
    test = Motor_P_Test()
    omega = [.01, 10, 100, 200, 400, 600, 800, 1000, 1200, 1400, 1600, 1800, 2000, 3000, 4000, 10000, 100000]
    eta = []
    for o in omega:
        test = Motor_P_Test()
        test.substitutions[test.mp.omega] = o
        sol = test.solve()
        eta.append(sol["freevariables"]["etam"])

    print omega
    print eta
    #sol = test.debug()
    plt.plot(omega, eta)
    plot.show()
    print sol.table()

def test():
    motor_test()
    propulsor_test()
    
if __name__ == "__main__":
    test()