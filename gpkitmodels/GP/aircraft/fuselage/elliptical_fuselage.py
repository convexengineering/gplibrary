" elliptical fuselage.py "
import numpy as np
from gpkit import Variable, Model, parse_variables
from gpkitmodels.GP.materials import cfrpfabric
from gpkitmodels import g

class FuselageAero(Model):
    """ Fuselage Aerodyanmic Model

    Variables
    ---------
    Cf              [-]         fuselage skin friction coefficient
    Re              [-]         fuselage reynolds number
    Cd              [-]         fuselage drag coefficient
    mfac    1.0     [-]         fuselage drag margin

    """
    def setup(self, static, state):
        exec parse_variables(FuselageAero.__doc__)

        V = state.V
        rho = state.rho
        l = static.l
        mu = state.mu
        k = static.k

        constraints = [
            Re == V*rho*l/mu,
            Cf >= 0.455/Re**0.3,
            Cd/mfac >= Cf*k
            ]

        return constraints

class Fuselage(Model):
    """ Fuselage Model

    Variables
    ---------
    R                   [ft]            fuselage radius
    l                   [ft]            fuselage length
    S                   [ft^2]          wetted fuselage area
    W                   [lbf]           fuselage weight
    mfac    2.0         [-]             fuselage weight margin factor
    f                   [-]             fineness ratio of lenth to diameter
    k                   [-]             fuselage form factor
    Vol                 [ft^3]          fuselae volume
    rhofuel 6.01        [lbf/gallon]    density of 100LL
    rhocfrp 1.6         [g/cm^3]        density of CFRP
    t                   [in]            fuselage skin thickness
    nply    2           [-]             number of plys

    """
    material = cfrpfabric
    flight_model = FuselageAero

    def setup(self):
        exec parse_variables(Fuselage.__doc__)

        rhocfrp = self.material.rho
        tmin = self.material.tmin

        constraints = [
            f == l/R/2,
            k >= 1 + 60/f**3 + f/400,
            3*(S/np.pi)**1.6075 >= 2*(l*R*2)**1.6075 + (2*R)**(2*1.6075),
            Vol == 4*np.pi/3*(l/2)*R**2,
            W/mfac >= S*rhocfrp*t*g,
            t >= nply*tmin,
            ]

        return constraints


