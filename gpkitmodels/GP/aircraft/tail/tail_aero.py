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

    Upper Bounded by state
    ----------------------
    \\rho, \\mu

    Lower Bounded by state
    ----------------------
    \\rho, \\mu, V

    Upper Unbounded
    ---------------
    Cd, Re, S, b

    Lower Unbounded
    ---------------
    S, tau

    LaTex Strings
    -------------
    Cd      C_d

    """
    def setup(self, static, state):
        exec parse_variables(TailAero.__doc__)

        cmac = self.cmac = static.planform.cmac
        b = self.b = static.planform.b
        S = self.S = static.planform.S
        tau = self.tau = static.planform.tau
        self.state = state
        rho = self.rho = state.rho
        V = self.V = state.V
        mu = self.mu = state.mu
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
