from gpkit import Model, Variable,Vectorize,parse_variables, SignomialsEnabled, SignomialEquality, units
from gpkit.constraints.tight import Tight as TCS
from propeller import BladeElementProp, Propeller
from gpkitmodels import g

class FlightState(Model):
    """ Flight State

    Variables
    ---------
    rho         1.225       [kg/m**3]       air density
    mu          1.789e-5    [N*s/m^2]       air viscosity
    V                       [kts]           speed
    qne                     [kg/s^2/m]      never exceed dynamic pressure
    Vne         160         [kts]           never exceed speed
    a           295         [m/s]           speed of sound
    """
    def setup(self):
        exec parse_variables(FlightState.__doc__)
        return [qne == 0.5*rho*Vne**2]

class TCSTest(Model):
    """ Test case for BE propeller model

    Variables
    ---------
    Mtip        .5          [-]         Max tip mach number
    omega_max   10000       [-]       maximum rotation rate
    m_total                  [-]        fds
    """
    def setup(self):
        exec parse_variables(TCSTest.__doc__)

        return [TCS([Mtip <= omega_max]),
                    m_total >= Mtip]
class PropTest(Model):
    """ Test case for BE propeller model

    Variables
    ---------
    Mtip                    [-]         Max tip mach number
    omega_max               [rpm]       maximum rotation rate
    T                       [N]       total thrust
    t_flight    100         [min]       nominal flight time
    Estar       200         [W*hr/kg]   battery specific energy
    P                       [W]         prop shaft power
    m_batt                  [kg]        battery weight
    V_cruise                [m/s]       cruise speed
    m_total                 [kg]        total system weight
    Q                       [N*m]       prop torque
    eta                     [-]       prop efficiency
    """
    
    def setup(self, N=5, MIL = False, DragModel = 1, f= .05, cl_spec = -1):
        exec parse_variables(PropTest.__doc__)
        self.prop = prop = Propeller(N, f)
        self.state = state = FlightState()
        prop.flight_model = BladeElementProp
        self.prop_perf = prop_perf = prop.flight_model(prop, state, N, MIL = MIL, DragModel = DragModel, cl_spec = cl_spec)


        Qprop = prop_perf.Q
        RPM_prop = prop_perf.omega
        Tprop = prop_perf.T

        constraints = [Tprop == T,
                        Qprop*RPM_prop == P,
                        eta == prop_perf.eta,
                        Q == Qprop,
                        P*t_flight/Estar == m_batt,
                        state.V == V_cruise,
                        m_total >= m_batt + prop.W/g]

        return constraints, prop, state, prop_perf


class MOPropTest(Model):
    """ Multi-objective propeller test case

    Variables
    ---------
    Mtip                    [-]         Max tip mach number
    omega_max               [rpm]       maximum rotation rate
    L_D_c       15          [-]         cruise L/D
    T_W_TO      .3          [-]         Takeoff Thrust-to-weight
    W           2700        [lbf]       vehicle MTOW      
    t_flight    5          [min]       nominal flight time
    t_TO        5           [min]       nominal takeoff time
    Estar       200         [W*hr/kg]   battery specific energy
    P_cruise                [W]         prop shaft power
    P_TO                    [W]         prop shaft power
    m_batt                  [kg]        battery weight
    V_cruise    120         [m/s]       cruise speed
    V_TO        10          [m/s]       cruise speed
    m_total                 [kg]        total system weight
    Rmax        1.5         [m]         Propeller radius
    """
    
    def setup(self, N=5, MIL = False, DragModel = 1, f= .05, cl_spec = -1):
        exec parse_variables(MOPropTest.__doc__)
        self.prop = prop = Propeller(N, f)
        self.TOstate = TOstate = FlightState()
        self.state = state = FlightState()
        prop.flight_model = BladeElementProp
        self.prop_perf = prop_perf = prop.flight_model(prop, state, N, MIL = False, DragModel = DragModel)
        self.TO_perf = TO_perf = prop.flight_model(prop, TOstate, N, MIL = False, DragModel = DragModel)

        Qprop = prop_perf.Q
        RPM_prop = prop_perf.omega
        Qprop_TO = TO_perf.Q
        RPM_prop_TO = TO_perf.omega
        Tprop_c = prop_perf.T
        Tprop_TO = TO_perf.T

        constraints = [Tprop_c == W/L_D_c,
                        Tprop_TO == W*T_W_TO,
                        Qprop*RPM_prop == P_cruise,
                        Qprop_TO*RPM_prop_TO== P_TO,
                        P_cruise*t_flight/Estar + P_TO*t_TO/Estar <= m_batt,
                        state.V == V_cruise,
                        TOstate.V == V_TO,
                        m_total >= m_batt + prop.W/g]
        return constraints, prop, TOstate, state, prop_perf, TO_perf

if __name__ == '__main__':
    #M = PropTest(N = 19, MIL = True, DragModel = 0)
    #M.cost = M.m_total
    #sol = M.localsolve()
    #print sol.table()
    M = MOPropTest()
    M.cost = M.m_total
    sol = M.localsolve()
    print sol.summary()



