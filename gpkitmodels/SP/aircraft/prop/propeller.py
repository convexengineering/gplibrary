" propeller model "
from numpy import pi
from gpkit import Model, Variable,Vectorize,parse_variables, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from gpfit.fit_constraintset import XfoilFit, FitCS
import os
import pandas as pd
from scipy.interpolate import interp1d

class QPROPDrag(Model):
    """QPROP DRAG MODEL

    Variables
    ---------
    cl                      [-]         lift coefficient
    cd                      [-]         local drag coefficient
    cd0          .01        [-]         base drag coefficient
    cd2           .008      [-]         drag rise with CL
    REref         150000    [-]         reference reynolds number
    re                      [-]         local reynolds number
    """

    def setup(self):
        exec parse_variables(QPROPDrag.__doc__)

        constraints = [cd >= (cd0 + cd2*cl**2)*(re/REref)**-.5]

        return constraints

class BladeElementPerf(Model):
    """ Single element of a propeller blade

    Three viscous drag models in use:
    0 - No viscous drag
    1 - pieceswise fit (QPROP methodology)
    2 - Xfoil airfoil fit

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
    r                       [m]         local radius
    lam_w                   [-]         advance ratio
    eps                     [-]         blade efficiency
    Re                      [-]         blade reynolds number
    f                       [-]         intermediate tip loss variable
    F                       [-]         Prandtl tip loss factor
    cl_max      5.0         [-]         max airfoil cl
    dr                      [m]         length of blade element
    M                       [-]         Mach number
    alpha                   [-]         local angle of attack
    beta_max    90          [-]         max twist angle
    alpha_max   20          [-]         stall AoA
    cl0         .0          [-]         zero AoA lift
    dclda       .1097       [-]         lift slope (per degree)
    h_ati                   [-]         atan^-1 helper variable
    h_aci                   [-]         acos(exp(-x))^-1 helper variable
    

    """

    def setup(self,static,  state, DragModel = 0):
        exec parse_variables(BladeElementPerf.__doc__)

        V       = state.V
        rho     = state.rho
        R       = static.R
        mu      = state.mu
        a       = state.a
        path = os.path.dirname(__file__)
        fd_cl = pd.read_csv(path + os.sep + "dae51_fitdata_cl_a.csv").to_dict(
            orient="records")[0]
        fd_cd = pd.read_csv(path + os.sep + "dae51_fitdata_cd_a.csv").to_dict(
            orient="records")[0]
        B = static.B
        c = static.c
        beta = static.beta

        constraints = [TCS([Wa>=V + va]),
                        TCS([Wt + vt<=omega*r]),
                        #Wt == omega*r,
                        TCS([G == (1./2.)*Wr*c*cl]),
                        h_aci**2.97951 >= (0.409975 * (f)**-1.45436
                                    + 0.173256 * (f)**0.157267),
                        #F == (2./pi)*(1.01116*f**.0379556)**(10), #This is the GP fit of arccos(exp(-f))
                        F == (2./pi)*(1/h_aci),
                        M == Wr/a,
                        lam_w == (r/R)*(Wa/Wt),
                        va == vt*(Wt/Wa),
                        TCS([dQ >= rho*B*G*(Wa+eps*Wt)*r*dr]),
                        Re == Wr*c*rho/mu,
                        eta_i == (V/(omega*r))*(Wt/Wa),
                        #eta_i*(1.+va/V) + 1 <= (vt/(omega*r)),
                        TCS([f+(r/R)*B/(2*lam_w) <= (B/2.)*(1./lam_w)]),                   
                        #TCS(XfoilFit(fd_cl, cl, [alpha,Re], name="clpolar")),
                        #beta <= beta_max,
                        #XfoilFit(fd_cd, cd, [cl,Re], name="cdpolar"),
                        #cd == 1e-30,
                        alpha <= alpha_max,
                        alpha >= .000001,
                        TCS([cl <= cl_max]),
                        
                    ]

        with SignomialsEnabled():
            constraints += [SignomialEquality(Wr**2,(Wa**2+Wt**2)),
                            #SignomialEquality(Wt + vt,omega*r),
                            #SignomialEquality(Wa,V+va),
                            #SignomialEquality(beta*pi/180., alpha*pi/180. + (0.946041 * (Wa/Wt)**0.996025)),
                            SignomialEquality(beta*pi/180., alpha*pi/180. + 1/h_ati),
                            
                            h_ati**1.83442 <= 0.966692 * (Wa/Wt)**-1.84391
                                        + 0.596688 * (Wa/Wt)**-0.0973781,
                            #SignomialEquality(cl**0.163122,0.237871 * (alpha)**-0.0466595 * (Re)**0.0255029
                            #                             + 0.351074 * (alpha)**0.255199 * (Re)**-0.021581
                            #                             + 0.224209 * (alpha)**-0.0427746 * (Re)**0.0258791),
                            SignomialEquality(cl, cl0 + alpha*dclda),
                            #cl <= cl0 + alpha*dclda,
                            #TCS(cl <=cl0 + alpha*dclda),
                            TCS([dT <= rho*B*G*(Wt-eps*Wa)*dr]),
                            TCS([vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2) >= (B*G/(4.*pi*r))**2]),
                            

            ]

        if DragModel == 0:
            constraints += [eps == 1e-30/cl]
            drag = []
        elif DragModel == 1:
            drag = self.drag = QPROPDrag()
            constraints += [drag.cl == cl,
                            drag.re == Re,
                            eps == drag.cd/cl]
        elif DragModel == 2:
            constraints += [XfoilFit(fd_cd, cd, [cl,Re], name="cdpolar"),
                            eps == drag.cd/cl]

        return constraints, drag

        



class BladeElementProp(Model):
    """ Performance for a propeller with multiple elements

    Variables
    ---------
    Mtip        10.5        [-]         Max tip mach number
    omega_max   100000      [rpm]       maximum rotation rate
    eta                     [-]         overall efficiency
    omega                   [rpm]       rotation rate
    T                       [lbf]       total thrust
    Q                       [N*m]       total torque
    """

    
    def setup(self,static,  state, N = 5, MIL = False, DragModel = 1, cl_spec = -1):
        exec parse_variables(BladeElementProp.__doc__)
        dr = static.dr
        r  = static.r
        with Vectorize(N):
            blade = self.blade = BladeElementPerf(static, state, DragModel)


        constraints = [blade.dr == dr*static.R,
                        blade.omega == omega,]

        #for n in range(1,N):
        #    constraints += [TCS([blade.r[n] >= blade.r[n-1] + static.R/N])]

        #    if MIL:    #Enforce MIL constraint at design condition
        #        constraints += [blade.eta_i[n] == blade.eta_i[n-1]]

        constraints += [TCS([Q >= sum(blade.dQ)]),
                        eta == state.V*T/(omega*Q),
                        blade.M[-1] <= Mtip,
                        static.T_m >= T,
                        omega <= omega_max
                        ]

        #Specify Design cls where desired
        if cl_spec != -1:
            XIdes = cl_spec[0]
            CLdes = cl_spec[1]
            fcl = interp1d(XIdes, CLdes)
            for i in range(N):
                cl_l = fcl(static.r_pc[i])
                constraints += [blade.cl[i] == float(cl_l)]

        with SignomialsEnabled():
            constraints += [TCS([T <= sum(blade.dT)])] 
            for n in range(N):
                constraints += [blade.r[n] == r[n]*static.R
                            #---TCS([blade.r[n] <= blade.r[n-1] + static.R/N])
                            ]
            for n in range(1,N):
                if MIL:    #Enforce MIL constraint at design condition
                    constraints += [blade.eta_i[n] == blade.eta_i[n-1]]

        return constraints, blade
class Propeller(Model):
    """ Propeller Model

    *Note: R_hub should be specified using f_hub (hub radius fraction) so that the control point distribution can be precomputed
    *       r/R is precomputed based on 1-f_hub and N
    Variables
    ---------
    R                               [m]             prop tip radius
    R_hub                           [m]             prop hub radius
    f_hub                           [-]             prop hub fraction
    W                               [lbf]           prop weight
    K           4e-4                [1/ft^2]        prop weight scaling factor
    T_m                             [lbf]           prop max static thrust
    B           2                   [-]             number of blades


    Variables of length N
    ---------------------
    c                               [m]             prop chord 
    beta                            [-]             local prop angle
    """

    flight_model = BladeElementProp

    def setup(self, N = 5, f = .05):
        exec parse_variables(Propeller.__doc__)
        self.N = N
        
        dR = (1-f)/N;
        self.r_pc = r_pc = [f+dR/2]
        for i in range(1,N):
            r_pc.append(r_pc[i-1]+dR)
        self.r = VectorVariable(N, "r", r_pc, "-", "r/R radial locations")
        self.dr = Variable("dr", dR, "-", "delta_r/R")

        constraints = [W >= K*T_m*R**2,
                        f_hub == R_hub/R,
                        f_hub == f]
        return constraints

