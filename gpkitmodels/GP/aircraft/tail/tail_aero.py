" tail aerodynamics "
import os
import pandas as pd
from gpkit import Model, parse_variables
from gpfit.fit_constraintset import XfoilFit

#pylint: disable=exec-used, attribute-defined-outside-init, undefined-variable
#pylint: disable=no-member

class TailAero(Model):
    """Tail Aero Model

    Variables
    ---------
    Re          [-]     Reynolds number
    Cd          [-]     drag coefficient

    Upper Unbounded
    ---------------
    Cd, Re
    static.planform.S, static.planform.b
    state.V, state.rho (if not rhovalue)

    Lower Unbounded
    ---------------
    static.planform.S, static.planform.tau, static.planform.b
    state.V, state.rho (if not rhovalue)

    LaTex Strings
    -------------
    Cd      C_d

    """
    def setup(self, static, state):
        self.state = state
        self.static = static
        exec parse_variables(TailAero.__doc__)

        cmac = static.planform.cmac
        b = static.planform.b
        S = static.planform.S
        tau = static.planform.tau
        rho = state.rho
        self.rhovalue = rho.key.value
        V = state.V
        mu = state.mu
        path = os.path.dirname(__file__)
        fd = pd.read_csv(path + os.sep + "tail_dragfit.csv").to_dict(
            orient="records")[0]

        constraints = [
            Re == V*rho*S/b/mu,
            # XfoilFit(fd, Cd, [Re, static["\\tau"]],
            #          err_margin="RMS", airfoil="naca 0008")
            XfoilFit(fd, Cd, [Re, tau], err_margin="RMS")
            ]

        return constraints
