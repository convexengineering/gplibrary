" wing.py "
import numpy as np
from gpkit import Variable, Model, Vectorize, SignomialsEnabled
from wing_interior import WingInterior
from wing_skin import WingSkin
from capspar import CapSpar
from tube_spar import TubeSpar
from constant_taper_chord import c_bar
from gpkitmodels.tools.fit_constraintset import FitCS
from gpkit.constraints.tight import Tight as TCS
import pandas as pd
import os

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
        mw = Variable("m_w", 2.0*np.pi/(1+2.0/23), "-",
                      "assumed span wise effectiveness")
        croot = Variable("c_{root}", "ft", "root chord")
        cmac = Variable("c_{MAC}", "ft", "mean aerodynamic chord")
        lamw = Variable("\\lambda", lam, "-", "wing taper ratio")
        cb, _ = c_bar(lam, N)
        with Vectorize(N):
            cbar = Variable("\\bar{c}", cb, "-",
                            "normalized chord at mid element")
        with Vectorize(N-1):
            cbave = Variable("\\bar{c}_{ave}", (cb[1:]+cb[:-1])/2, "-",
                             "normalized mid section chord")
            cave = Variable("c_{ave}", "ft", "mid section chord")

        constraints = [b**2 == S*AR,
                       lamw == lamw,
                       cbar == cbar,
                       cave == cbave*S/b,
                       croot == S/b*cb[0],
                       cmac == S/b]

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

class WingLoading(Model):
    "wing loading cases"
    def setup(self, wing, Wcent, Wwing=None, V=None, CL=None):

        loading = [wing.wingskin.loading()]
        loading.append(wing.spar.loading(Wcent))
        if Wwing:
            loading.append(wing.spar.gustloading(Wcent, Wwing, V, CL))

        return loading

class WingAero(Model):
    "wing aerodynamic model with profile and induced drag"
    def setup(self, static, state):
        "wing drag model"
        Cd = Variable("C_d", "-", "wing drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        e = Variable("e", 0.9, "-", "span efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        cdp = Variable("c_{dp}", "-", "wing profile drag coeff")

        path = os.path.abspath(__file__).replace(os.path.basename(__file__), "")
        df = pd.read_csv(path + os.sep + "jho_fitdata.csv")

        constraints = [
            Cd >= cdp + CL**2/np.pi/static["AR"]/e,
            Re == state["\\rho"]*state["V"]*static["c_{MAC}"]/state["\\mu"],
            FitCS(df, cdp, [CL, Re], airfoil="jho1.dat")
            # FitCS(df, cdp, [CL, Re])
            ]

        return constraints
