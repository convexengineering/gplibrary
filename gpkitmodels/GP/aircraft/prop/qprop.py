" qprop model "
from numpy import pi
from gpkit import Model, parse_variables
from gpkit.constraints.tight import Tight as TCS
from gpfit.fit_constraintset import FitCS
import pandas as pd
from os.path import abspath

class QProp(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [N]         thrust
    B          2            [-]         number of blades
    r          0.5          [m]         local radius
    R          1            [m]         propeller radius
    dr         1            [m]         local delta radius
    Q                       [N*m]       torque
    G                       [m**2/s]    prop circulation
    Omega     1000          [rpm]     rotation rate
    eta                     [-]         prop efficiency
    W                       [m/s]       local total velocity
    Wa                      [m/s]       axial local total velocity
    Wt                      [m/s]       tangential local total velocity
    va                      [m/s]       axial local rotor-induced velocity
    vt                      [m/s]       tangential local rotor-induced velocity
    ep         0.1          [-]         inverse local lift to drag
    cl         1.0          [-]         local coefficient of lift
    c          0.1          [m]         local blade chord
    F          0.5          [-]         Prandtl factor
    lamw                    [-]         local wake advance ratio
    f                       [-]         Prandtl exponent

    """

    def setup(self, state):
        exec parse_variables(QProp.__doc__)

        V = state.V
        rho = state.rho

        df = pd.read_csv(abspath("arccosfit.csv"))
        fd = df.to_dict(orient="records")[0]

        return [eta <= V*T/Omega/Q,
                Q >= rho*B*G*(Wa + ep*Wt)*r*dr,
                G >= 0.5*W*c*cl,
                W**2 >= Wa**2 + Wt**2,
                Wa >= V + va,
                va >= vt*Wt/Wa,
                Wt >= Omega*r,
                # vt**2*(1 + (4*lamw*R/pi/B/r)**2) >= (B*G/4/pi/r/F)**2,
                vt >= (B*G/4/pi/r/F),
                # lamw >= r/R*(vt/va),
                # f >= B/2*(0.5)*R/r*(va/vt),
                # FitCS(fd, F, [f])
               ]
