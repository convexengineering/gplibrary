" wing.py "
from os import sep
from os.path import abspath, dirname
import numpy as np
import pandas as pd
from gpkit import Variable, Model, parse_variables
from .wing_interior import WingInterior
from .wing_skin import WingSkin
from .capspar import CapSpar
from gpfit.fit_constraintset import XfoilFit

#pylint: disable=no-member, invalid-name, unused-argument, exec-used
#pylint: disable=undefined-variable

class Planform(Model):
    """ Planform Area Definition

    Scalar Variables
    ---------
    S                                   [ft^2]  surface area
    AR                                  [-]     aspect ratio
    b                                   [ft]    span
    tau             0.115               [-]     airfoil thickness ratio
    CLmax           1.39                [-]     maximum lift coefficient
    CM              0.14                [-]     wing moment coefficient
    croot                               [ft]    root chord
    cmac                                [ft]    mean aerodynamic chord
    lam             0.5                 [-]     taper ratio
    cbarmac         self.return_cmac    [-]     non-dim MAC

    Variables of length N
    ---------------------
    eta         np.linspace(0,1,N)      [-]     (2y/b)
    cbar        self.return_c           [-]     non-dim chord at nodes

    Variables of length N-1
    -----------------------
    cave                                [ft]    mid section chord
    cbave       self.return_avg         [-]     non-dim mid section chord
    deta        self.return_deta        [-]     \\Delta (2y/b)

    Upper Unbounded
    ---------------
    S, AR, b, croot, cmac

    Lower Unbounded
    ---------------
    S, AR, b, croot, cmac

    LaTex Strings
    -------------
    tau         \\tau
    CLmax       C_{L_{\\mathrm{max}}}
    CM          C_M
    croot       c_{\\mathrm{root}}
    cmac        c_{\\mathrm{MAC}}
    lam         \\lambda
    cbarmac     \\bar{c}_{\\mathrm{MAC}}

    """
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

    return_avg = lambda self, c: (self.return_c(c)[:-1]
                                  + self.return_c(c)[1:])/2.
    return_deta = lambda self, c: np.diff(c[self.eta])

    def setup(self, N):
        exec parse_variables(Planform.__doc__)

        return [b**2 == S*AR,
                cave == cbave*S/b,
                croot == S/b*cbar[0],
                cmac == croot*cbarmac]

class WingAero(Model):
    """
    Wing Aero Model

    Variables
    ---------
    Cd                      [-]     wing drag coefficient
    CL                      [-]     lift coefficient
    CLstall         1.3     [-]     stall CL
    e               0.9     [-]     span efficiency
    Re                      [-]     reynolds number
    cdp                     [-]     wing profile drag coefficient

    Upper Unbounded
    ---------------
    Cd, Re

    Lower Unbounded
    ---------------
    CL

    LaTex Strings
    -------------
    Cd              C_d
    CL              C_L
    CLstall         C_{L_{\\mathrm{stall}}}
    cdp             c_{d_p}

    """
    def setup(self, static, state,
              fitdata=dirname(abspath(__file__)) + sep + "jho_fitdata.csv"):
        exec parse_variables(WingAero.__doc__)

        df = pd.read_csv(fitdata)
        fd = df.to_dict(orient="records")[0]

        if fd["d"] == 2:
            independentvars = [self.CL, self.Re]
        elif fd["d"] == 3:
            independentvars = [self.CL, self.Re, static["\\tau"]]

        return [Cd >= cdp + CL**2/np.pi/static.planform.AR/e,
                Re == (state["\\rho"]*state["V"]*static.planform.cmac
                       / state["\\mu"]),
                # XfoilFit(fd, cdp, [CL, Re], airfoil="jho1.dat"),
                XfoilFit(fd, cdp, independentvars),
                CL <= CLstall
               ]

class Wing(Model):
    """
    Wing Model

    Variables
    ---------
    W                   [lbf]       wing weight
    mfac        1.2     [-]         wing weight margin factor
    """

    sparModel = CapSpar
    fillModel = WingInterior
    flight_model = WingAero

    def setup(self, N=5):
        exec parse_variables(Wing.__doc__)

        self.N = N

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
