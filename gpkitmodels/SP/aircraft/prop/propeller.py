" propeller model "
from numpy import pi
from gpkit import Model, Variable,Vectorize,parse_variables, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from gpfit.fit_constraintset import XfoilFit, FitCS
import os
import pandas as pd



class BladeElementPerf(Model):
    """ Single element of a propeller blade

    Variables
    ---------
    dT                      [lbf]       thrust
    eta_i                   [-]         local induced efficiency
    dQ                      [N*m]       torque
    omega                   [rpm]       propeller rotation rate 
    Wa                      [m/s]       Axial total relative velocity
    Wt                      [m/s]       Tangential total relative velocity
    Wr                      [m/s]       Total relative velocity
    va                      [m/s]       Axial induced velocity
    vt                      [m/s]       Tangential induced velocity
    G                       [m^2/s]     Circulation
    cl                      [-]         local lift coefficient
    cd                      [-]         local drag coefficient
    B           2           [-]         number of blades
    r                       [m]         local radius
    lam_w                   [-]         advance ratio
    eps                     [-]         blade efficiency
    AR_b                    [-]         blade aspect ratio
    AR_b_max    50          [-]         max blade aspect ratio
    Re                      [-]         blade reynolds number
    f                       [-]         intermediate tip loss variable
    F                       [-]         Prandtl tip loss factor
    cl_max      .7          [-]         max airfoil cl
    dr                      [m]         length of blade element
    M                       [-]         Mach number
    a           295         [m/s]       Speed of sound at altitude
    alpha                   [-]         local angle of attack
    beta_max    1.05        [-]         max twist angle

    """

    def setup(self,static,  state):
        exec parse_variables(BladeElementPerf.__doc__)

        V       = state.V
        rho     = state.rho
        R       = static.R
        mu      = state.mu
        path = os.path.dirname(__file__)
        fd_cl = pd.read_csv(path + os.sep + "dae51_fitdata_cl_a.csv").to_dict(
            orient="records")[0]
        fd_cd = pd.read_csv(path + os.sep + "dae51_fitdata.csv").to_dict(
            orient="records")[0]

        c = static.c
        beta = static.beta
        constraints = [TCS([Wa>=V + va]),
                        TCS([Wt + vt<=omega*r]),
                        TCS([G == (1./2.)*Wr*c*cl]),
                        F == (2./pi)*(1.01116*f**.0379556)**(10), #This is the GP fit of arccos(exp(-f))
                        M == Wr/a,
                        lam_w == (r/R)*(Wa/Wt),
                        va == vt*(Wt/Wa),
                        eps == cd/cl,
                        TCS([dQ >= rho*B*G*(Wa+eps*Wt)*r*dr]),
                        AR_b == R/c,
                        AR_b <= AR_b_max,
                        Re == Wr*c*rho/mu,
                        eta_i == (V/(omega*r))*(Wt/Wa),
                        TCS([f+(r/R)*B/(2*lam_w) <= (B/2.)*(1./lam_w)]),                   
                        #TCS(XfoilFit(fd_cl, cl, [alpha,Re], name="clpolar")),
                        #beta <= beta_max
                        XfoilFit(fd_cd, cd, [cl,Re], name="cdpolar"),
                        #TCS([beta >= alpha + Wa/Wt]), 
                        #TCS([cl**0.135559 >= 0.250583 * (alpha*180./pi)**0.214389 * (Re)**-0.0391487
                        #    + 0.251724 * (alpha*180./pi)**0.190369 * (Re)**-0.032165
                        #    + 0.351677 * (alpha*180./pi)**-0.051975 * (Re)**0.0428966]),
                        cl <= cl_max
                    ]
        with SignomialsEnabled():
            constraints += [SignomialEquality(Wr**2,(Wa**2+Wt**2)),
                            SignomialEquality(beta, alpha + (0.946041 * (Wa/Wt)**0.996025)),
                            SignomialEquality(cl**0.135559,0.250583 * (alpha*180./pi)**0.214389 * (Re)**-0.0391487
                            + 0.251724 * (alpha*180./pi)**0.190369 * (Re)**-0.032165
                            + 0.351677 * (alpha*180./pi)**-0.051975 * (Re)**0.0428966),
                            TCS([dT <= rho*B*G*(Wt-eps*Wa)*dr]),
                            TCS([vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2) >= (B*G/(4.*pi*r))**2]),

            ]
        return constraints, state



class BladeElementProp(Model):
    """ Performance for a propeller with multiple elements

    Variables
    ---------
    Mtip        .5          [-]         Max tip mach number
    omega_max   10000       [rpm]       maximum rotation rate
    eta                     [-]         overall efficiency
    omega                   [rpm]       rotation rate
    T                       [lbf]       total thrust
    Q                       [N*m]       total torque
    """
    def setup(self,static,  state, N = 5):
        exec parse_variables(BladeElementProp.__doc__)
        
        with Vectorize(N):
            blade = BladeElementPerf(static, state)


        constraints = [blade.dr == static.R/(N),
                        blade.omega == omega,
                        blade.r[0] == static.R/(2.*N)]

        for n in range(1,N):
            constraints += [TCS([blade.r[n] >= blade.r[n-1] + static.R/N]),
                            blade.eta_i[n] == blade.eta_i[n-1],
                            ]

        constraints += [TCS([Q >= sum(blade.dQ)]),
                        eta == state.V*T/(omega*Q),
                        blade.M[-1] <= Mtip,
                        static.T_m >= T,
                        omega <= omega_max
                        ]

        with SignomialsEnabled():
            constraints += [TCS([T <= sum(blade.dT)])] 

        return constraints, blade 

