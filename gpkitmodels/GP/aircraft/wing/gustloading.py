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
    v, cl, wing.spar.I, wing.spar.tshear, wing.spar.Sy
    wing.spar.J (if wingSparJ)
    theta (if not wingSparJ), M (if not wingSparJ)

    Lower Unbounded
    ---------------
    Ww, wing.planform.b, wing.planform.cave
    state.qne (if wingSparJ)
    theta (if not wingSparJ), M (if not wingSparJ)

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
        [1e-10, 1-cos(c[self.wing.planform.eta][1:]*pi/2)])

    def setup(self, wing, state):
        exec parse_variables(GustL.__doc__)
        self.load = SparLoading.setup(self, wing, state)

        cbar = self.wing.planform.cbar
        W = self.W  # from SparLoading

        path = os.path.dirname(os.path.abspath(__file__))
        df = pd.read_csv(path + os.sep + "arctan_fit.csv").to_dict(
            orient="records")[0]

        constraints = [
            # fit for arctan from 0 to 1, RMS = 0.044
            FitCS(df, agust, [cosminus1*vgust/v]),
            self.beam["\\bar{q}"] >= cbar*(1 + 2*pi*agust/cl*(1+Ww/W)),
            ]

        return self.load, constraints
