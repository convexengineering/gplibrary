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
    def return_c(self, c):
        " return normalized chord distribution "
        return np.array([2./(1+c[self.lam])*(1+(c[self.lam]-1)*e) for e
                         in c[self.eta]])

    def return_cmac(self, c):
        " return normalized MAC "
        cbar = self.return_c(c)
        lam = cbar[1:]/cbar[:-1]
        maci = 2./3*cbar[:-1]*(1 + lam + lam**2)/(1 + lam)
        deta = np.diff(c[self.eta])
        num = sum([(cbar[i] + cbar[i+1])/2*maci[i]*deta[i] for i
                   in range(len(deta))])
        den = sum([(cbar[i] + cbar[i+1])/2*deta[i] for i in range(len(deta))])
        return num/den/cbar[0]

    def setup(self, N):

        S = Variable("S", "ft^2", "surface area")
        AR = Variable("AR", "-", "aspect ratio")
        b = Variable("b", "ft", "wing span")
        tau = Variable("\\tau", 0.115, "-", "airfoil thickness ratio")
        CLmax = Variable("C_{L_{max}}", 1.39, "-", "maximum CL of JHO1")
        CM = Variable("C_M", 0.14, "-", "wing moment coefficient")
        croot = Variable("c_{root}", "ft", "root chord")
        cmac = Variable("c_{MAC}", "ft", "mean aerodynamic chord")
        self.lam = Variable("\\lambda", 0.5, "-", "wing taper ratio")
        with Vectorize(N):
            self.eta = Variable("\\eta", np.linspace(0, 1, N), "-", "(2y/b)")
            cbar = Variable("\\bar{c}", self.return_c, "-", "non-dim chord at nodes")

        with Vectorize(N-1):
            cave = Variable("c_{ave}", "ft", "mid section chord")
            return_avg = lambda c: (self.return_c(c)[:-1] + self.return_c(c)[1:])/2.
            cbave = Variable("\\bar{c}_{ave}", return_avg, "-",
                             "non-dim mid section chord")
            return_deta = lambda c: np.diff(c[self.eta])
            deta = Variable("d\\eta", return_deta, "-", "\\Delta (2y/b)")
        cbarmac = Variable("\\bar{c}_{MAC}", self.return_cmac, "-",
                           "non-dim MAC")

        return [b**2 == S*AR,
                cave == cbave*S/b,
                croot == S/b*cbar[0],
                cmac == croot*cbarmac]

class WingAero(Model):
    "wing aerodynamic model with profile and induced drag"
    def setup(self, static, state, fitdata="jho_fitdata.csv"):
        "wing drag model"
        Cd = Variable("C_d", "-", "wing drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        CLstall = Variable("C_{L_{stall}}", 1.3, "-", "stall CL")
        e = Variable("e", 0.9, "-", "span efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        cdp = Variable("c_{dp}", "-", "wing profile drag coeff")

        path = os.path.dirname(os.path.abspath(fitdata))
        df = pd.read_csv(path + os.sep + fitdata)
        fd = df.to_dict(orient="records")[0]

        if fd["d"] == 2:
            independentvars = [CL, Re]
        elif fd["d"] == 3:
            independentvars = [CL, Re, static["\\tau"]]

        constraints = [
            Cd >= cdp + CL**2/np.pi/static["AR"]/e,
            Re == state["\\rho"]*state["V"]*static["c_{MAC}"]/state["\\mu"],
            # XfoilFit(fd, cdp, [CL, Re], airfoil="jho1.dat"),
            XfoilFit(fd, cdp, independentvars),
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

    def setup(self, N=5, lam=0.5):
        # TODO: phase out lam in later version

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
