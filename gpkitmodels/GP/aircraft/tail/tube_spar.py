" tube spar "
from numpy import pi
from gpkitmodels.GP.materials import cfrpfabric
from gpkitmodels import g
from gpkit import Model, parse_variables

class TubeSpar(Model):
    """ Tail Boom Model

    Variables
    ---------
    mfac        1.0             [-]         weight margin factor
    k           0.8             [-]         taper index
    kfac        self.minusk2    [-]         (1-k/2)
    W                           [lbf]       spar weight

    Variables of length N-1
    -----------------------
    I                           [m^4]       moment of inertia
    d                           [in]        diameter
    t                           [in]        thickness
    dm                          [kg]        segment mass
    Sy                          [m^3]       section modulus

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    J, l, I0

    LaTex Strings
    -------------
    kfac        (1-k/2)
    mfac        m_{\\mathrm{fac}}

    """

    minusk2 = lambda self, c: 1-c(self.k)/2.
    material = cfrpfabric

    @parse_variables(__doc__, globals())
    def setup(self, N, surface):
        deta = surface.deta
        tmin = self.material.tmin
        rho = self.material.rho
        l = surface.l

        self.weight = W/mfac >= g*dm.sum()

        return [I <= pi*t*d**3/8.0,
                Sy <= 2*I/d,
                dm >= pi*rho*d*deta*t*kfac*l,
                self.weight,
                t >= tmin]
