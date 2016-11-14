" tail boom flexibility "
import numpy as np
from gpkit import Model, Variable

class TailBoomFlexibility(Model):
    "tail boom flexibility model"
    def __init__(self, htail, tailboom, wing, state, **kwargs):

        Fne = Variable("F_{NE}", "-", "tail boom flexibility factor")
        deda = Variable("d\\epsilon/d\\alpha", "-", "wing downwash derivative")
        SMcorr = Variable("SM_{corr}", 0.35, "-", "corrected static margin")

        # signomial helper variables
        sph1 = Variable("sph1", "-", "first term involving $V_h$")
        sph2 = Variable("sph2", "-", "second term involving $V_h$")

        constraints = [
            Fne >= (1 + htail["m_h"]*0.5*state["V_{NE}"]**2*state["\\rho_{sl}"]
                    * htail["S"]*tailboom["l"]**2/tailboom["E"]
                    / tailboom["I_0"]*tailboom["(1-k/2)"]),
            sph1*(wing["m_w"]*Fne/htail["m_h"]/htail["V_h"]) + deda <= 1,
            sph2 <= htail["V_h"]*htail["(C_{L_h})_{min}"]/wing["C_{L_{max}}"],
            (sph1 + sph2).mono_lower_bound({"sph1": .48, "sph2": .52}) >= (
                SMcorr + wing["C_M"]/wing["C_{L_{max}}"]),
            deda >= wing["m_w"]*wing["S"]/wing["b"]/4/np.pi/htail["l_h"]]

        Model.__init__(self, None, constraints, **kwargs)
