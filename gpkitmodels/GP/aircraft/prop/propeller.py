" propeller model "
from numpy import pi
from gpkit import Model, Variable,Vectorize,parse_variables, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from gpfit.fit_constraintset import XfoilFit
import os
import pandas as pd

class ActuatorProp(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [lbf]       thrust
    Tc                      [-]         coefficient of thrust
    etaadd     .7           [-]         swirl and nonuniformity losses
    etav       .85          [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                     [-]         overall efficiency
    z1         self.helper  [-]         efficiency helper 1
    z2                      [-]         efficiency helper 2
    lam                     [-]         advance ratio
    CT                      [-]         thrust coefficient
    CP                      [-]         power coefficient
    Q                       [N*m]       torque
    omega                   [rpm]       propeller rotation rate
    omega_max  10000        [rpm]       max rotation rate
    P_shaft                 [kW]        shaft power
    M_tip      .5           [-]         Tip mach number
    a          295          [m/s]       Speed of sound at altitude
    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self, static, state):
        exec parse_variables(ActuatorProp.__doc__)

        V = state.V
        rho = state.rho
        R = static.R

        constraints = [eta <= etav*etai,
                       Tc >= T/(0.5*rho*V**2*pi*R**2),
                       z2 >= Tc + 1,
                       etai*(z1 + z2**0.5/etaadd) <= 2,
                       lam >= V/(omega*R),
                       CT >= Tc*lam**2,
                       CP <= Q*omega/(.5*rho*(omega*R)**3*pi*R**2),
                       eta >= CT*lam/CP,
                       omega <= omega_max,
                       P_shaft == Q*omega,
                       (M_tip*a)**2 >= (omega*R)**2 + V**2,
                       static.T_m >= T
                      ]
        return constraints

class Blade_Element_Performance(Model):
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
    cl                     [-]          local lift coefficient
    cd                      [-]         local drag coefficient
    B           2           [-]         number of blades
    r                       [m]         local radius
    lam_w                   [-]         advance ratio
    eps                     [-]         blade efficiency
    AR_b                   [-]          blade aspect ratio
    AR_b_max    50          [-]         max blade aspect ratio
    Re                      [-]         blade reynolds number
    f                       [-]         intermediate tip loss variable
    F                       [-]         Prandtl tip loss factor
    cl_max      .6         [-]          max airfoil cl
    dr                      [m]         length of blade element
    M                       [-]         Mach number
    a           295         [m/s]       Speed of sound at altitude

    """

    def setup(self,static,  state):
        exec parse_variables(Blade_Element_Performance.__doc__)

        V       = state.V
        rho     = state.rho
        R       = static.R
        mu      = state.mu
        path = os.path.dirname(__file__)
        fd = pd.read_csv(path + os.sep + "dae51_fitdata.csv").to_dict(
            orient="records")[0]
        c = static.c
        constraints = [TCS([Wa>=V + va]),
                        TCS([Wt + vt<=omega*r]),
                        TCS([G == (1./2.)*Wr*c*cl]),
                        F == (2./pi)*(1.01116*f**.0379556)**(1./.1), #This is the GP fit of arccos(exp(-f))
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
                        
                        XfoilFit(fd, cd, [cl,Re], name="polar"),
                        
                        cl <= cl_max
                    ]
        with SignomialsEnabled():
            constraints += [SignomialEquality(Wr**2,(Wa**2+Wt**2)),
                            
                            TCS([dT <= rho*B*G*(Wt-eps*Wa)*dr]),
                            #SignomialEquality(vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2),(B*G/(4.*pi*r))**2),
                            TCS([vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2)>=(B*G/(4.*pi*r))**2]),
                            #SignomialEquality(f+(r/R)*B/(2*lam_w),(B/2.)*(1./lam_w))
                            
                    
            ]
        return constraints, state



class MultiElementProp(Model):
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
        exec parse_variables(MultiElementProp.__doc__)
        
        with Vectorize(N):
            blade = Blade_Element_Performance(static, state)


        constraints = [blade.dr == static.R/(N),
                        blade.omega == omega,
                        blade.r[0] == static.R/(2.*N)]

        
        with SignomialsEnabled():
            for n in range(1,N):
                constraints += [TCS([blade.r[n] >= blade.r[n-1] + static.R/N]),
                                blade.eta_i[n] == blade.eta_i[n-1],
                                ]
            constraints += [TCS([T <= sum(blade.dT)]),
                            TCS([Q >= sum(blade.dQ)]),
                            eta == state.V*T/(omega*Q),
                            blade.M[-1] <= Mtip,
                            static.T_m>=T,
                            omega <= omega_max
                            ]


        return constraints, blade 



class Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    R                               [ft]            prop radius
    W                               [lbf]           prop weight
    K           4e-4                [1/ft^2]        prop weight scaling factor
    T_m                             [lbf]           prop max static thrust

    Variables of length N
    ---------------------
    c                               [ft]            prop chord 
    """

    flight_model = ActuatorProp

    def setup(self, N = 5):
        exec parse_variables(Propeller.__doc__)
        self.N = N
        return [W >= K*T_m*R**2]
