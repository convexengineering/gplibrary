" wing.py "
import os
import numpy as np
import pandas as pd
from gpkit import Variable, Model, Vectorize
from .wing_interior import WingInterior
from .wing_skin import WingSkin
from .capspar import CapSpar
from .constant_taper_chord import c_bar
from gpfit.fit_constraintset import XfoilFit

#pylint: disable=invalid-name, attribute-defined-outside-init, unused-variable
#pylint: disable=too-many-instance-attributes

class Wing(Model):
    """
    Aicraft wing model for constant tapered wing
    INPUTS
    ------
    N : int             number of sections
    lam : float         taper ratio
    hollow: boolean     True if wing is not hollow (filled with foam)
    """
    def setup(self, N=5, lam=0.5, hollow=False):

        self.N = N

        W = Variable("W", "lbf", "wing weight")
        mfac = Variable("m_{fac}", 1.2, "-", "wing weight margin factor")

        cb, eta, deta, cbarmac = c_bar(lam, N)
        subdict = {"\\lambda": lam, "\\bar{c}": cb, "\\eta": eta,
                   "\\bar{c}_{ave}": (cb[1:]+cb[:-1])/2,
                   "\\bar{c}_{MAC}": cbarmac, "d\\eta": deta}

        self.surf = AeroSurf(N)
        self.surf.substitutions.update(subdict)
        self.spar = CapSpar(N, self.surf)
        self.skin = WingSkin()
        self.components = [self.spar, self.skin]


        constraints = [
            W/mfac >= sum(c["W"] for c in self.components),
            self.skin["W"] >= (self.skin["\\rho_{CFRP}"]*self.surf["S"]*2
                               * self.skin["t"]*self.skin["g"]),
            ]

        if not hollow:
            self.foam = WingInterior()
            self.components.extend([self.foam])
            constraints.extend([
                self.foam["W"] >= 2*(
                    self.foam["g"]*self.foam["\\rho_{foam}"]
                    * self.foam["\\bar{A}_{jh01}"]*self.surf["c_{ave}"]**2
                    * (self.surf["b"]/2)*self.surf["d\\eta"]).sum()])

        self.flight_model = WingAero
        self.loading = WingLoading

        return constraints, self.surf, self.components

class AeroSurf(Model):
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
        lamw = Variable("\\lambda", "-", "wing taper ratio")
        cbarmac = Variable("\\bar{c}_{MAC}", "-", "non-dim MAC")
        with Vectorize(N):
            cbar = Variable("\\bar{c}", "-",
                            "normalized chord at mid element")
            eta = Variable("\\eta", "-", "(2y/b)")
        with Vectorize(N-1):
            cbave = Variable("\\bar{c}_{ave}", "-",
                             "normalized mid section chord")
            cave = Variable("c_{ave}", "ft", "mid section chord")
            deta = Variable("d\\eta", "-", "\\Detla (2y/b)")

        constraints = [b**2 == S*AR,
                       cave == cbave*S/b,
                       croot == S/b*cbar[0],
                       cmac == croot*cbarmac,
                       cbar == cbar,
                       eta == eta,
                       deta == deta,
                       lamw == lamw]

        return constraints

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
