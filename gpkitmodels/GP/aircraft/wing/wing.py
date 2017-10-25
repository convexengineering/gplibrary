" wing.py "
import os
import numpy as np
import pandas as pd
from gpkit import Variable, Model, Vectorize
from .wing_interior import WingInterior
from .wing_skin import WingSkin
from .capspar import CapSpar
from gpfit.fit_constraintset import XfoilFit

#pylint: disable=invalid-name, attribute-defined-outside-init, unused-variable
#pylint: disable=too-many-instance-attributes, too-many-locals, no-member

class Planform(Model):
    "The thing that creates the lift"
    def setup(self, N):

        S = Variable("S", "ft^2", "surface area")
        AR = Variable("AR", "-", "aspect ratio")
        b = Variable("b", "ft", "wing span")
        tau = Variable("\\tau", 0.115, "-", "airfoil thickness ratio")
        CLmax = Variable("C_{L_{max}}", 1.39, "-", "maximum CL of JHO1")
        CM = Variable("C_M", 0.14, "-", "wing moment coefficient")
        croot = Variable("c_{root}", "ft", "root chord")
        cmac = Variable("c_{MAC}", "ft", "mean aerodynamic chord")
        lam = Variable("\\lambda", 0.5, "-", "wing taper ratio")
        return_cmac = lambda c: 2./3*(1+c[lam]+c[lam]**2)/(1+c[lam])
        cbarmac = Variable("\\bar{c}_{MAC}", return_cmac, "-", "non-dim MAC")
        with Vectorize(N):
            eta = Variable("\\eta", np.linspace(0, 1, N), "-", "(2y/b)")
            return_c = lambda c: np.array([2./(1+c[lam])*(1+(c[lam]-1)*e)
                                           for e in c[eta]])
            cbar = Variable("\\bar{c}", return_c, "-", "non-dim chord at nodes")

        with Vectorize(N-1):
            cave = Variable("c_{ave}", "ft", "mid section chord")
            return_avg = lambda c: (return_c(c)[:-1] + return_c(c)[1:])/2.
            cbave = Variable("\\bar{c}_{ave}", return_avg, "-",
                             "non-dim mid section chord")
            return_deta = lambda c: np.diff(c[eta])
            deta = Variable("d\\eta", return_deta, "-", "\\Delta (2y/b)")

        return [b**2 == S*AR,
                cave == cbave*S/b,
                croot == S/b*cbar[0],
                cmac == croot*cbarmac]

class WingLoading(Model):
    "wing loading cases"
    def setup(self, wing, Wcent, Wwing=None, V=None, CL=None):

        loading = [wing.skin.loading(wing)]
        loading.append(wing.spar.loading(wing, Wcent))
        if Wwing:
            loading.append(wing.spar.gustloading(wing, Wcent, Wwing, V, CL))

        return loading

class WingAero(Model):
    "wing aerodynamic model with profile and induced drag"
    def setup(self, static, state):
        "wing drag model"
        Cd = Variable("C_d", "-", "wing drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        CLstall = Variable("C_{L_{stall}}", 1.3, "-", "stall CL")
        e = Variable("e", 0.9, "-", "span efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        cdp = Variable("c_{dp}", "-", "wing profile drag coeff")

        path = os.path.dirname(__file__)
        df = pd.read_csv(path + os.sep + "jho_fitdata.csv")
        fd = df.to_dict(orient="records")[0]

        constraints = [
            Cd >= cdp + CL**2/np.pi/static["AR"]/e,
            Re == state["\\rho"]*state["V"]*static["c_{MAC}"]/state["\\mu"],
            # XfoilFit(fd, cdp, [CL, Re], airfoil="jho1.dat"),
            XfoilFit(fd, cdp, [CL, Re]),
            CL <= CLstall
            ]

        return constraints

class Wing(Model):
    """
    Aicraft wing model for constant tapered wing
    INPUTS
    ------
    N : int             number of sections
    lam : float         taper ratio
    hollow: boolean     True if wing is not hollow (filled with foam)
    """

    sparModel = CapSpar
    fillModel = WingInterior
    flight_model = WingAero
    loading = WingLoading

    def setup(self, N=5):

        self.N = N

        W = Variable("W", "lbf", "wing weight")
        mfac = Variable("m_{fac}", 1.2, "-", "wing weight margin factor")

        self.planform = Planform(N)
        self.skin = WingSkin(self.planform)
        self.components = [self.skin]

        if self.sparModel:
            self.spar = self.sparModel(N, self.planform)
            self.components.extend([self.spar])
        if self.fillModel:
            self.foam = self.fillModel(self.planform)
            self.components.extend([self.foam])

        constraints = [W/mfac >= sum(c["W"] for c in self.components)]

        return constraints, self.planform, self.components
