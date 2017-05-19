" wing.py "
import numpy as np
from gpkit import Variable, Model, SignomialsEnabled, ConstraintSet
from gpkit.constraints.array import ArrayConstraint
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP

#pylint: disable=attribute-defined-outside-init, invalid-name

def clean_model(model):
    for vk in model.varkeys:
        if len(vk.models) > 1:
            delete_model(vk)
    for cs in model:
        if isinstance(cs, ArrayConstraint):
            vvars = [cs.right, cs.left]
            for v in vvars:
                if "descr" in v.__dict__.keys():
                    delete_model(v)
        elif isinstance(cs, ConstraintSet):
            for c in cs:
                if isinstance(c, ArrayConstraint):
                    vvars = [c.right, c.left]
                    for v in vvars:
                        if "descr" in v.__dict__.keys():
                            delete_model(v)
                elif isinstance(c, Model):
                    delete_model(c)
                    clean_model(c)
        elif isinstance(cs, Model):
            delete_model(cs)
            clean_model(cs)

def delete_model(vk):
    if isinstance(vk, Model):
        parentm = vk.name[0]
        if parentm in vk.name[1:]:
            ind = vk.name[1:].index(parentm)
            del vk.name[ind+1]
            del vk.num[ind+1]
    else:
        parentm = vk.models[0]
        if parentm in vk.models[1:]:
            ind = vk.models[1:].index(parentm)
            del vk.models[ind+1]
            del vk.modelnums[ind+1]

class Wing(Model):
    "The thing that creates the lift"
    def setup(self, N=5, lam=0.5, spar="CapSpar", hollow=False):

        self.gpwing = WingGP(N=N, lam=lam, spar=spar, hollow=hollow)
        clean_model(self.gpwing)
        del self.gpwing.name[0]
        del self.gpwing.num[0]
        mw = Variable("m_w", "-", "span wise effectiveness")

        with SignomialsEnabled():
            constraints = [mw*(1 + 2/self.gpwing["AR"]) >= 2*np.pi]

        self.spar = self.gpwing.spar
        self.wingskin = self.gpwing.wingskin
        self.components = self.gpwing.components

        if not hollow:
            self.winginterior = self.gpwing.winginterior

        return constraints, self.gpwing

    def flight_model(self, state):
        " what happens during flight "
        return self.gpwing.flight_model(state)

    def loading(self, Wcent, Wwing=None, V=None, CL=None):
        " loading cases "
        return self.gpwing.loading(Wcent, Wwing, V, CL)

