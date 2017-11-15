" wing.py "
import numpy as np
from gpkit.constraints.array import ArrayConstraint
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkit import Variable, SignomialsEnabled

#pylint: disable=attribute-defined-outside-init, invalid-name

class Wing(WingGP):

    def setup(self, N=5):

        self.wing = WingGP.setup(self, N=N)
        mw = Variable("m_w", "-", "span wise effectiveness")

        with SignomialsEnabled():
            constraints = [mw*(1 + 2/self.planform["AR"]) >= 2*np.pi]

        return self.wing, constraints
