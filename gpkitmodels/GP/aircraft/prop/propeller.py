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

    def setup(self, static, state = None, rho = None, V = None):
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
    eta_b                   [-]         local efficiency
    dQ                      [N*m]       torque
    omega                   [rpm]       propeller rotation rate 
    Wa                      [m/s]       Axial total relative velocity
    Wt                      [m/s]       Tangential total relative velocity
    Wr                      [m/s]       Total relative velocity
    va                      [m/s]       Axial induced velocity
    vt                      [m/s]       Tangential induced velocity
    G                       [m^2/s]     Circulation
    c                       [m]         local chord
    cl                     [-]         local lift coefficient
    cd                      [-]         local drag coefficient
    cdp                     [-]         local drag coefficient
    B           3           [-]         number of blades
    r                       [m]         local radius
    lam_w                   [-]         advance ratio
    eps                     [-]         blade efficiency
    AR_b                   [-]         blade aspect ratio
    AR_b_max    50          [-]         max blade aspect ratio
    Ut                      [m/s]         tangential freestreem
    Re                      [-]         blade reynolds number
    f                       [-]         intermediate tip loss variable
    F                       [-]         Prandtl tip loss factor
    cl_max      .9         [-]         max airfoil cl
    dr                      [m]         length of blade element
    FF                      [-]         form factor
    t_c         .1        [-]         airfoil t/c
    C_f                     [-]         skin friction coefficient
    M                       [-]         Mach number
    a           295         [m/s]       Speed of sound at altitude

    """

    def setup(self,parent,  state):
        exec parse_variables(Blade_Element_Performance.__doc__)

        V       = state.V
        rho     = state.rho
        R       = parent.R
        mu      = state.mu
        path = os.path.dirname(__file__)
        fd = pd.read_csv(path + os.sep + "dae51_fitdata.csv").to_dict(
            orient="records")[0]

        constraints = [#TCS([Wr**2 >= (Wa**2+Wt**2)]),
                        #Ut == omega*r,
                        TCS([G == (1./2.)*Wr*c*cl]),
                        #TCS([f+(r/R)*B/(2*lam_w) <= (B/2.)*(1./lam_w)]), #.2 = 1-.8 = 1-r/R
                        F == (2./pi)*(1.02232*f**.0458729)**(1./.1), #This is the GP fit of arccos(exp(-f))
                        M == Wr/a,
                        lam_w == (r/R)*(Wa/Wt),
                        va == vt*(Wt/Wa),
                        eps == cd/cl,
                        #TCS([1 >= etai*(1+va/V)+vt/Ut]),
                        #TCS([1 >= etap*(1+eps*Wt/Wa)+eps*Wt/Wa]),
                        TCS([dQ >= rho*B*G*(Wa+eps*Wt)*r*dr]),
                        #Wa == V,
                        #Wt == omega*r,
                        AR_b == R/c,
                        AR_b <= AR_b_max,
                        Re == Wr*c*rho/mu,
                        #eta == T*V/(Q*omega),
                        eta_b == (V/(omega*r))*(Wt/Wa),

                        #TCS([cd >= cl**2/(pi*AR_b*.85) + cdp]),
                        #TCS([FF >= 1. + .3*t_c + (100*t_c)**.4]),
                        #C_f == .027/(Re**(1./7.)),
                        #cdp >= FF*C_f*(2.+.4*t_c),
                        
                        XfoilFit(fd, cd, [cl,Re], name="polar"),
                        
                        cl <= cl_max
                    ]
        with SignomialsEnabled():
            constraints += [TCS([Wa>=V + va]),
                            #SignomialEquality(Wa, V+va),
                            SignomialEquality(Wr**2,(Wa**2+Wt**2)),
                            #TCS([Wr**2 >= (Wa**2+Wt**2)]),
                            SignomialEquality(Wt,omega*r-vt),
                            #TCS([Wt<=omega*r-vt]),
                            TCS([dT <= rho*B*G*(Wt-eps*Wa)*dr]),
                            SignomialEquality(vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2),(B*G/(4.*pi*r))**2),
                            #TCS([vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2)<=(B*G/(4.*pi*r))**2]),
                            #TCS([f+(r/R)*B/(2*lam_w) <= (B/2.)*(1./lam_w)]),
                            SignomialEquality(f+(r/R)*B/(2*lam_w),(B/2.)*(1./lam_w))
                            #SignomialEquality(1,etap*(1+eps*Wt/Wa)+eps*Wa/Wt),
                            #SignomialEquality(1,etai*(1+va/V)+vt/Ut),
                            #SignomialEquality(dQ,rho*B*G*(Wa+eps*Wt)*r*dr),
                            #SignomialEquality(dT,rho*B*G*(Wt-eps*Wa)*dr),
                    
            ]
        return constraints, state



class Multi_Element_Propeller_Performance(Model):
    """ Performance for a propeller with multiple elements

    Variables
    ---------
    Mtip        .5          [-]         tip mach number
    """
    def setup(self,parent,  state, N = 4):
        #exec parse_variables(Multi_Element_Propeller_Performance.__doc__)
        T = self.T = Variable('T', 'lbf', 'Overall thrust')
        Q = self.Q = Variable('Q', 'N*m', 'Overall torque')
        eta = self.eta = Variable('eta', '-', 'Overall efficiency')
        omega = self.omega =  Variable('omega', 'rpm', 'rotation rate')
        omega_max  =  Variable('omega_max',10000, 'rpm', 'rotation rate')
        Mtip       =  Variable('Mtip', .5, '-', 'Tip Mach')
        with Vectorize(N):
            blade = Blade_Element_Performance(parent, state)


        constraints = [blade.dr == parent.R/(N),
                        blade.omega == omega,
                        blade["r"][0] == parent.R/(2.*N)]

        
        with SignomialsEnabled():
            for n in range(1,N):
                constraints += [TCS([blade["r"][n] >= blade["r"][n-1] + parent.R/N]),
                                #SignomialEquality(blade["r"][n],blade["r"][n-1] + parent.R/N),
                                blade["eta_b"][n] == blade["eta_b"][n-1]
                                ]
            constraints += [TCS([T <= sum(blade["dT"][n] for n in range(N))]),
                            TCS([Q >= sum(blade["dQ"][n] for n in range(N))]),
                            eta == state.V*T/(omega*Q),
                            blade["M"][-1] <= Mtip,
                            #parent.T_m>=T
                            omega <= omega_max
                            ]


        return constraints, blade
"hi"
class SimpleQProp(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [lbf]       thrust
    etaadd     0.7          [-]         swirl and nonuniformity losses
    etap                    [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                      [-]         overall efficiency
    Q                       [N*m]       torque
    omega                   [rpm]     propeller rotation rate
    omega_max    10000      [rpm]     max rotation rate
    P_shaft                 [kW]        shaft power
    M_tip       .5          [-]         Tip mach number
    a           295         [m/s]       Speed of sound at altitude
    Wa                      [m/s]       Axial total relative velocity
    Wt                      [m/s]       Tangential total relative velocity
    Wr                       [m/s]       Total relative velocity
    va                      [m/s]       Axial induced velocity
    vt                      [m/s]       Tangential induced velocity
    G                       [m^2/s]     Circulation
    c                       [m]         local chord
    cl         .4             [-]         local lift coefficient
    cd                      [-]         local drag coefficient
    B           2           [-]         number of blades
    r                       [m]         local radius
    lam_w                   [-]         advance ratio
    eps                     [-]         blade efficiency
    AR_b                    [-]         blade aspect ratio
    AR_b_max    18          [-]         max blade aspect ratio
    Ut                      [m/s]         tangential freestreem
    Re                      [-]         blade reynolds number
    f                       [-]         intermediate tip loss variable
    F                       [-]         Prandtl tip loss factor
    cl_max      1.3         [-]         max airfoil cl
    eta_check              [-]         efficiency check
    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self, static, state):
        exec parse_variables(SimpleQProp.__doc__)

        V = state.V
        rho = state.rho
        R = static.R
        mu = state.mu
        path = os.path.dirname(__file__)
        fd = pd.read_csv(path + os.sep + "dae51_fitdata.csv").to_dict(
            orient="records")[0]

        constraints = [TCS([eta <= etap*etai]),
                        omega <= omega_max,
                        P_shaft == Q*omega,
                        Ut == omega*r,
                        TCS([G == (1./2.)*Wr*c*cl]),
                        r == .8*R,      #Assume 80% chord is representative
                        f == (B/2.)*.2*(1./lam_w), #.2 = 1-.8 = 1-r/R
                        F == (2./pi)*(1.02232*f**.0458729)**(1./.1), #This is the GP fit of arccos(exp(-f))
                        lam_w == (r/R)*(Wa/Wt),
                        va == vt*(Wt/Wa),
                        eps == cd/cl,
                        TCS([1 >= etai*(1+va/V)+vt/Ut]),
                        TCS([1 >= etap*(1+eps*Wt/Wa)+eps*Wt/Wa]),
                        TCS([Q >= rho*B*G*(Wa+eps*Wt)*R**2]),
                        
                        TCS([(M_tip*a)**2 >= Wa**2 + Wt**2]),
                        AR_b == R/c,
                        AR_b <= AR_b_max,
                        Re == Wr*c*rho/mu,
                        eta == T*V/(Q*omega),
                        XfoilFit(fd, cd, [cl,Re], name="polar"),
                        static.T_m >= T,
                        cl <= cl_max
                    ]
        with SignomialsEnabled():
            constraints += [TCS([Wa>=V + va]),
                            SignomialEquality(Wr**2,(Wa**2+Wt**2)),
                            #SignomialEquality(Wa,V+va),
                            SignomialEquality(Wt,omega*r-vt),
                            #TCS([Wt>=omega*r-vt]),
                            TCS([T <= rho*B*G*(Wt-eps*Wa)*R]),
                            SignomialEquality(vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2),(B*G/(4.*pi*r))**2),
                            #TCS([vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2)>=(B*G/(4.*pi*r))**2]),
                            #SignomialEquality(1,etap*(1+eps*Wt/Wa)+eps*Wa/Wt),
                            #SignomialEquality(1,etai*(1+va/V)+vt/Ut),
                            #SignomialEquality(Q,rho*B*G*(Wa+eps*Wt)*R**2),
                            #SignomialEquality(T,rho*B*G*(Wt-eps*Wa)*R),
            ]
        return constraints, state

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

    flight_model = SimpleQProp

    def setup(self, N = 1):
        exec parse_variables(Propeller.__doc__)

        return [W >= K*T_m*R**2]



class Multi_Element_Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    R                               [ft]            prop radius
    W                               [lbf]           prop weight
    K           4e-4                [1/ft^2]        prop weight scaling factor
    T_m                             [lbf]           prop max static thrust
    """

    flight_model = Multi_Element_Propeller_Performance
    def setup(self):
        exec parse_variables(Multi_Element_Propeller.__doc__)

        constraints = [W >= K*T_m*R**2]

        return constraints
