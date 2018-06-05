" wing interior "
from gpkit import Model, parse_variables
from gplibrary.GP.materials import foamhd
from gplibrary import g

#pylint: disable=exec-used, no-member, undefined-variable

class WingCore(Model):
    """ Wing Core Model

    Variables
    ---------
    W                           [lbf]       wing core weight
    Abar            0.0753449   [-]         normalized cross section area

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    cave, b, surface.deta

    LaTex Strings
    -------------
    rhocore                 \\rho_{\\mathrm{core}}
    Abar                    \\bar{A}

    """
    material = foamhd

    def setup(self, surface):
        self.surface = surface
        exec parse_variables(WingCore.__doc__)

        cave = self.cave = surface.cave
        b = self.b = surface.b
        deta = surface.deta
        rho = self.material.rho

        return [W >= 2*(g*rho*Abar*cave**2*b/2*deta).sum()]
