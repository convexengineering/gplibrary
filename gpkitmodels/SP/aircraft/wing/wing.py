" wing.py "
import numpy as np
from gpkit import Variable, Model, SignomialsEnabled, ConstraintSet
from gpkit.constraints.array import ArrayConstraint
from gpkitmodels.GP.aircraft.wing.wing import Wing as WingGP
from gpkit.varkey import VarKey
from gpkit.keydict import KeySet
from gpkit import Variable, Model, Vectorize, SignomialsEnabled
from gpkitmodels.GP.aircraft.wing.wing_interior import WingInterior
from gpkitmodels.GP.aircraft.wing.wing_skin import WingSkin
from gpkitmodels.GP.aircraft.wing.capspar import CapSpar
from gpkitmodels.GP.aircraft.wing.tube_spar import TubeSpar
from gpkitmodels.GP.aircraft.wing.constant_taper_chord import c_bar
from gpkitmodels.GP.aircraft.wing.wing import WingAero, WingLoading

#pylint: disable=attribute-defined-outside-init, invalid-name

def clean_model(model, parentn):
    for vk in model.varkeys:
        if len(vk.models) > 1:
            delete_model(vk, parentn)
    for cs in model:
        if isinstance(cs, ArrayConstraint):
            vvars = [cs.right, cs.left]
            for v in vvars:
                if "descr" in v.__dict__.keys():
                    delete_model(v, parentn)
        elif isinstance(cs, ConstraintSet):
            for c in cs:
                if isinstance(c, ArrayConstraint):
                    vvars = [c.right, c.left]
                    for v in vvars:
                        if "descr" in v.__dict__.keys():
                            delete_model(v, parentn)
                elif isinstance(c, Model):
                    delete_model(c, parentn)
                    clean_model(c, parentn)
        elif isinstance(cs, Model):
            delete_model(cs, parentn)
            clean_model(cs, parentn)

def delete_model(vk, parentn):
    if isinstance(vk, Model):
        ind = vk.name.index(parentn)
        if parentn in vk.name[ind+1:]:
            del vk.name[ind+1]
            del vk.num[ind+1]
    else:
        ind = vk.models.index(parentn)
        if parentn in vk.models[ind+1:]:
            del vk.models[ind+1]
            del vk.modelnums[ind+1]

# class Wing(Model):
#     "The thing that creates the lift"
#     def __init__(self, **kwargs):
#         super(Wing, self).__init__(**kwargs)
#         clean_model(self.gpwing, "Wing")
#         self.gpwing.name = []
#         self.gpwing.num = []
#
#         newvks = KeySet()
#         for vk in self.gpwing.varkeys:
#             newvk = VarKey(**vk.descr)
#             newvks[newvk] = vk
#
#         newvks.update(self.unique_varkeys)
#         self.varkeys = newvks
#
#     def setup(self, N=5, lam=0.5, spar="CapSpar", hollow=False):
#
#         self.gpwing = WingGP(N=N, lam=lam, spar=spar, hollow=hollow)
#
#         mw = Variable("m_w", "-", "span wise effectiveness")
#
#         with SignomialsEnabled():
#             constraints = [mw*(1 + 2/self.gpwing["AR"]) >= 2*np.pi]
#
#         self.spar = self.gpwing.spar
#         self.wingskin = self.gpwing.wingskin
#         self.components = self.gpwing.components
#
#         if not hollow:
#             self.winginterior = self.gpwing.winginterior
#
#         return constraints, self.gpwing
#
#     def flight_model(self, state):
#         " what happens during flight "
#         return self.gpwing.flight_model(state)
#
#     def loading(self, Wcent, Wwing=None, V=None, CL=None):
#         " loading cases "
#         return self.gpwing.loading(Wcent, Wwing, V, CL)

class Wing(Model):
    "The thing that creates the lift"
    def setup(self, N=5, lam=0.5, spar="CapSpar", hollow=False):

        W = Variable("W", "lbf", "weight")
        mfac = Variable("m_{fac}", 1.2, "-", "wing weight margin factor")
        S = Variable("S", "ft^2", "surface area")
        AR = Variable("AR", "-", "aspect ratio")
        b = Variable("b", "ft", "wing span")
        tau = Variable("\\tau", 0.115, "-", "airfoil thickness ratio")
        CLmax = Variable("C_{L_{max}}", 1.39, "-", "maximum CL of JHO1")
        CM = Variable("C_M", 0.14, "-", "wing moment coefficient")
        croot = Variable("c_{root}", "ft", "root chord")
        cmac = Variable("c_{MAC}", "ft", "mean aerodynamic chord")
        lamw = Variable("\\lambda", lam, "-", "wing taper ratio")
        mw = Variable("m_w", "-", "span wise effectiveness")

        cb, _ = c_bar(lam, N)
        with Vectorize(N):
            cbar = Variable("\\bar{c}", cb, "-",
                            "normalized chord at mid element")
        with Vectorize(N-1):
            cbave = Variable("\\bar{c}_{ave}", (cb[1:]+cb[:-1])/2, "-",
                             "normalized mid section chord")
            cave = Variable("c_{ave}", "ft", "mid section chord")

        constraints = [b**2 == S*AR,
                       cave == cbave*S/b,
                       croot == S/b*cb[0],
                       cmac == S/b]

        with SignomialsEnabled():
            constraints.extend([mw*(1 + 2/AR) >= 2*np.pi])

        if spar == "CapSpar":
            self.spar = CapSpar(b, cave, tau, N)
        elif spar == "TubeSpar":
            self.spar = TubeSpar(b, cave, tau, N)
        self.wingskin = WingSkin(S, croot, b)
        self.components = [self.spar, self.wingskin]

        if not hollow:
            self.winginterior = WingInterior(cave, b, N)
            self.components.extend([self.winginterior])

        constraints.extend([W/mfac >= sum(c["W"] for c in self.components)])

        return self.components, constraints

    def flight_model(self, state):
        return WingAero(self, state)

    def loading(self, Wcent, Wwing=None, V=None, CL=None):
        return WingLoading(self, Wcent, Wwing, V, CL)
