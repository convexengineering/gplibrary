" propeller model "
from numpy import pi
from gpkit import Model, Variable,Vectorize,parse_variables, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from gpfit.fit_constraintset import XfoilFit
import os
import pandas as pd

class ActuatorProp(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [lbf]       thrust
    Tc                      [-]         coefficient of thrust
    etaadd     .7           [-]         swirl and nonuniformity losses
    etav       .85          [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                     [-]         overall efficiency
    z1         self.helper  [-]         efficiency helper 1
    z2                      [-]         efficiency helper 2
    lam                     [-]         advance ratio
    CT                      [-]         thrust coefficient
    CP                      [-]         power coefficient
    Q                       [N*m]       torque
    omega                   [rpm]       propeller rotation rate
    omega_max  10000        [rpm]       max rotation rate
    P_shaft                 [kW]        shaft power
    M_tip      .5           [-]         Tip mach number
    a          295          [m/s]       Speed of sound at altitude
    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self, static, state = None, rho = None, V = None):
        exec parse_variables(ActuatorProp.__doc__)

        V = state.V
        rho = state.rho
        R = static.R

        constraints = [eta <= etav*etai,
                       Tc >= T/(0.5*rho*V**2*pi*R**2),
                       z2 >= Tc + 1,
                       etai*(z1 + z2**0.5/etaadd) <= 2,
                       lam >= V/(omega*R),
                       CT >= Tc*lam**2,
                       CP <= Q*omega/(.5*rho*(omega*R)**3*pi*R**2),
                       eta >= CT*lam/CP,
                       omega <= omega_max,
                       P_shaft == Q*omega,
                       (M_tip*a)**2 >= (omega*R)**2 + V**2,
                       static.T_m >= T
                      ]
        return constraints



class Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    R                               [ft]            prop radius
    W                               [lbf]           prop weight
    K           4e-4                [1/ft^2]        prop weight scaling factor
    T_m                             [lbf]           prop max static thrust

    Variables of length N
    ---------------------
    c                               [ft]            prop chord 
    """

    flight_model = ActuatorProp

    def setup(self, N = 1):
        exec parse_variables(Propeller.__doc__)

        return [W >= K*T_m*R**2]




