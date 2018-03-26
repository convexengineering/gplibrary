" propeller model "
from numpy import pi
from gpkit import Model, parse_variables, SignomialsEnabled, SignomialEquality

class Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    R          10           [m]         prop radius
    W          .01          [lbf]        prop weight

    """
    def setup(self):
        exec parse_variables(Propeller.__doc__)

    def flight_model(self,state):
        return Propeller_Performance(self, state)

class Propeller_Performance(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [N]         thrust
    Tc                      [-]         coefficient of thrust
    etaadd     0.7          [-]         swirl and nonuniformity losses
    etav       0.85         [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                     [-]         overall efficiency
    z1         self.helper  [-]         efficiency helper 1
    z2                      [-]         efficiency helper 2
    lam                     [-]         advance ratio
    CT                      [-]         thrust coefficient
    CP                      [-]         power coefficient
    Q                       [N*m]       torque
    omega                   [RPM]       propeller rotation rate 

    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self,parent,  state):
        exec parse_variables(Propeller_Performance.__doc__)

        V       = state.V
        rho     = state.rho
        R       = parent.R        

        constraints = [eta <= etav*etai,
                Tc == T/(0.5*rho*V**2*pi*R**2),
                z2 >= Tc + 1,
                etai*(z1 + z2**0.5/etaadd) <= 2,
                lam == V/(omega*R),
                CT == Tc*lam**2,
                #CP == Q*omega/(.5*rho*(omega*R)**3*math.pi*R**2),
                #eta == CT*lam/CP, #Is this the same eta? check 
                ]
        return constraints