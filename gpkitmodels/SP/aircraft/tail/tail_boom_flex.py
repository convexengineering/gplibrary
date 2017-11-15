" tail boom flexibility "
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled

class TailBoomFlexibility(Model):
    "tail boom flexibility model"
    def setup(self, htail, tailboom, wing, state):

        Fne = Variable("F_{NE}", "-", "tail boom flexibility factor")
        deda = Variable("d\\epsilon/d\\alpha", "-", "wing downwash derivative")
        SMcorr = Variable("SM_{corr}", 0.55, "-", "corrected static margin")

        # signomial helper variables
        sph1 = Variable("sph1", "-", "first term involving $V_h$")
        sph2 = Variable("sph2", "-", "second term involving $V_h$")

        constraints = [
            Fne >= (1 + htail.mh*0.5*state.Vne**2*state.rhosl
                    * htail["S"]*htail.lh**2/tailboom.E
                    / tailboom.I0*tailboom.kfac),
            sph1*(wing["m_w"]*Fne/htail.mh/htail.Vh) + deda <= 1,
            sph2 <= htail.Vh*htail.CLhmin/wing.planform.CLmax,
            # (sph1 + sph2).mono_lower_bound({"sph1": .48, "sph2": .52}) >= (
            #     SMcorr + wing["C_M"]/wing["C_{L_{max}}"]),
            deda >= wing["m_w"]*wing["S"]/wing["b"]/4/np.pi/htail.lh]

        with SignomialsEnabled():
            constraints.extend([sph1 + sph2 >= SMcorr + wing.planform.CM/wing.planform.CLmax])

        return constraints
