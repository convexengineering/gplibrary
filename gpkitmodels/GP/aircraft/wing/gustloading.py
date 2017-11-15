" spar loading for gust case "
import os
from numpy import pi, hstack, cos
import pandas as pd
from gpkit import parse_variables
from gpfit.fit_constraintset import FitCS
from .sparloading import SparLoading

#pylint: disable=invalid-name, no-member, arguments-differ, exec-used
#pylint: disable=attribute-defined-outside-init, undefined-variable

class GustL(SparLoading):
    """ Gust Loading Model

    Variables
    ---------
    vgust       10      [m/s]       gust velocity
    Ww                  [lbf]       wing weight
    v                   [m/s]       vehicle speed
    cl                  [-]         wing lift coefficient

    Variables of length wing.N
    --------------------------
    agust                           [-]         gust angle of attack
    cosminus1   self.return_cosm1   [-]         1 minus cosine factor

    Upper Unbounded
    ---------------
    v, cl, I, tshear, Sy, cave

    Lower Unbounded
    ---------------
    Ww, b

    LaTex Strings
    -------------
    vgust               V_{\\mathrm{gust}}
    Ww                  W_{\\mathrm{w}}
    cl                  c_l
    agust               \\alpha_{\\mathrm{gust}}
    cosminus1           (cos(x)-1)

    """
    new_qbarFun = None
    new_SbarFun = None

    return_cosm1 = lambda self, c: hstack(
        [1e-3, 1-cos(c[self.wing.planform.eta][1:]*pi/2)])

    def setup(self, wing):
        exec parse_variables(GustL.__doc__)
        self.load = SparLoading.setup(self, wing)

        self.b = self.wing.planform.b
        self.I = self.wing.spar.I
        self.Sy = self.wing.spar.Sy
        self.cave = self.wing.planform.cave
        self.tshear = self.wing.spar.tshear

        path = os.path.dirname(os.path.abspath(__file__))
        df = pd.read_csv(path + os.sep + "arctan_fit.csv").to_dict(
            orient="records")[0]

        constraints = [
            # fit for arctan from 0 to 1, RMS = 0.044
            FitCS(df, agust, [cosminus1*vgust/v]),
            self.beam["\\bar{q}"] >= self.wing.planform.cbar*(
                1 + 2*pi*agust/cl*(1+Ww/self.W)),
            ]

        return self.load, constraints
