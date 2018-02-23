" propeller model "
from numpy import pi
from gpkit import Model, parse_variables

class Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [N]         thrust
    Tc                      [-]         coefficient of thrust
    R          0.3          [m]         prop radius
    etaadd     0.7          [-]         swirl and nonuniformity losses
    etav       0.85         [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                     [-]         overall efficiency
    z1         self.helper  [-]         efficiency helper 1
    z2                      [-]         efficiency helper 2

    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self, state):
        exec parse_variables(Propeller.__doc__)

        V = state.V
        rho = state.rho

        return [eta <= etav*etai,
                Tc == T/(0.5*rho*V**2*pi*R**2),
                z2 >= 1 + Tc,
                etai*(z1 + z2**0.5/etaadd) <= 2]
