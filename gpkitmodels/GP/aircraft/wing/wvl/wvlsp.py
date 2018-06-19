import numpy as np
from os.path import abspath, dirname
from os import sep
from wvlfun import wvl
from gpkit import Variable, Model, VectorVariable, SignomialsEnabled, parse_variables
from gpkit.constraints.tight import Tight as TCS
from gpkitmodels.GP.aircraft.wing.wing import WingAero
from gpkitmodels.GP.aircraft.wing import wing

class WVL(WingAero):
    """ weissenger vortex lifting line

    Variables of length Na
    ---------------------
    dy                      [m]         spanwise section length
    G                       [m]         Gamma/Vinf vortex filament strength

    Variables of length [Na,Na]
    ---------------------------
    B       self.Bmatrix    [-]         B matrix

    """

    Bmatrix = wvl

    def setup(self, static, state,
              fitdata=dirname(abspath(wing.__file__)) + sep + "jho_fitdata.csv"):
        self.wingaero = super(WVL, self).setup(static, state, fitdata)
        Na = (static.N-1)*2
        exec parse_variables(WVL.__doc__)

        S = self.static.planform.S
        dylist = np.array(list(np.flip(static.planform.deta, 0)) + list(static.planform.deta))

        for dd in dylist:
            print type(dd)

        with SignomialsEnabled():
            self.Di = TCS([self.CDi >= 2*np.dot(G, np.dot(B, G))/S])
            self.constraints.extend([
                dy == dylist*static.planform.b/2,
                self.Di,
                TCS([self.CL <= sum(2*G*dy/S)])
                ])

        return self.wingaero, self.constraints
