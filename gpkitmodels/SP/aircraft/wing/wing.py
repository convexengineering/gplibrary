" wing.py "
import numpy as np
from gpkit.constraints.array import ArrayConstraint
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkit import Variable, SignomialsEnabled

#pylint: disable=attribute-defined-outside-init, invalid-name

def Wing(N=5, lam=0.5, spar="CapSpar", hollow=False):

    wing = WingGP(N=N, lam=lam, spar=spar, hollow=hollow)
    mw = Variable("m_w", "-", "span wise effectiveness")

    with SignomialsEnabled():
        constraints = [mw*(1 + 2/wing["AR"]) >= 2*np.pi]

    wing.append(constraints[0])
    wing.substitutions.update(constraints[0].substitutions)
    return wing
