" wing.py "
import numpy as np
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkit import parse_variables, SignomialsEnabled

#pylint: disable=attribute-defined-outside-init, invalid-name

class Wing(WingGP):
    """ SP wing model

    Variables
    ---------
    mw          [-]     span wise effectiveness

    """
    def setup(self, N=5):
        exec(parse_variables(Wing.__doc__))

        self.wing = WingGP.setup(self, N=N)
        with SignomialsEnabled():
            constraints = [mw*(1 + 2/self.planform["AR"]) >= 2*np.pi]

        return self.wing, constraints
