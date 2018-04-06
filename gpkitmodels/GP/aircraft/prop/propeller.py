" propeller model "
from numpy import pi
from gpkit import Model, parse_variables, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight as TCS
from gpfit.fit_constraintset import XfoilFit
import os
import pandas as pd


class Actuator_Propeller_Performance(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [lbf]       thrust
    Tc                      [-]         coefficient of thrust
    etaadd     0.7          [-]         swirl and nonuniformity losses
    etav       0.85         [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                     [-]         overall efficiency
    z1         self.helper  [-]         efficiency helper 1
    z2                      [-]         efficiency helper 2
    lam                     [-]         advance ratio
    CT                      [-]         thrust coefficient
    CP                      [-]         power coefficient
    Q                       [N*m]       torque
    omega                   [rpm]     propeller rotation rate 
    omega_max    10000      [rpm]     max rotation rate
    P_shaft                 [kW]        shaft power
    M_tip       .5          [-]         Tip mach number
    a           295         [m/s]       Speed of sound at altitude
    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self,parent,  state):
        exec parse_variables(Actuator_Propeller_Performance.__doc__)

        V       = state.V
        rho     = state.rho
        R       = parent.R


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
                (M_tip*a)**2 >= (omega*R)**2 + V**2

                ]
        return constraints, state

class Actuator_Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    R                               [ft]            prop radius
    W                               [lbf]           prop weight
    K           4e-4                [1/ft^2]        prop weight scaling factor
    T_m_static   40.                [lbf]           prop max static thrust
    """

    flight_model = Actuator_Propeller_Performance
    def setup(self):
        exec parse_variables(Actuator_Propeller.__doc__)

        constraints = [W >= K*T_m_static*R**2]

        return constraints

class One_Element_Propeller_Performance(Model):
    """ Propeller Model

    Variables
    ---------
    T                       [lbf]       thrust
    etaadd     0.7          [-]         swirl and nonuniformity losses
    etap                    [-]         viscous losses
    etai                    [-]         inviscid losses
    eta                     [-]         overall efficiency
    Q                       [N*m]       torque
    omega                   [rpm]     propeller rotation rate 
    omega_max    10000      [rpm]     max rotation rate
    P_shaft                 [kW]        shaft power
    M_tip       .5          [-]         Tip mach number
    a           295         [m/s]       Speed of sound at altitude
    Wa                      [m/s]       Axial total relative velocity
    Wt                      [m/s]       Tangential total relative velocity
    W                       [m/s]       Total relative velocity
    va                      [m/s]       Axial induced velocity
    vt                      [m/s]       Tangential induced velocity
    G                       [m^2/s]     Circulation
    c                       [m]         local chord
    cl                      [-]         local lift coefficient
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
    """

    def helper(self, c):
        return 2. - 1./c[self.etaadd]

    def setup(self,parent,  state):
        exec parse_variables(One_Element_Propeller_Performance.__doc__)

        V       = state.V
        rho     = state.rho
        R       = parent.R
        mu      = state.mu
        path = os.path.dirname(__file__)
        fd = pd.read_csv(path + os.sep + "dae51_fitdata.csv").to_dict(
            orient="records")[0]

        constraints = [TCS([eta == etap*etai]),
                        omega <= omega_max,
                        P_shaft == Q*omega,
                        TCS([W**2 >= (Wa**2+Wt**2)]),
                        Ut == omega*r,
                        TCS([G == (1./2.)*W*c*cl]),
                        r == .8*R,      #Assume 80% chord is representative
                        f == (B/2.)*.2*(1./lam_w), #.2 = 1-.8 = 1-r/R
                        F == (2./pi)*(1.02232*f**.0458729)**(1./.1),
                        #vt**2. >= (B*G/(4.*pi*r))**2*(1./(F**2.*(1.+(4.*lam_w*R/(pi*B*r))**2.))),
                        #TCS([vt == (B*G/(4.*pi*r*F))]),
                        lam_w == (r/R)*(Wa/Wt),
                        va == vt*(Wt/Wa),
                        eps == cd/cl,
                        TCS([1 >= etai*(1+va/V)+vt/Ut]),
                        TCS([1 >= etap*(1+eps*Wt/Wa)+eps*Wt/Wa]),
                        TCS([Q >= rho*B*G*(Wa+eps*Wt)*R**2]),
                        
                        #TCS([(M_tip*a)**2 >= Wa**2 + Wt**2]),
                        AR_b == R/c,
                        AR_b <= AR_b_max,
                        Re == V*c*rho/mu,
                        XfoilFit(fd, cd, [cl,Re], name="polar"),
                        parent.T_m >= T

                    ]
        with SignomialsEnabled():
            constraints += [TCS([Wa<=V + va]),
                            #SignomialEquality(Wt,omega*r-vt),
                            TCS([Wt>=omega*r-vt]),
                            TCS([T <= rho*B*G*(Wt-eps*Wa)*R]),
                            #SignomialEquality(vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2),(B*G/(4.*pi*r))**2),
                            TCS([vt**2*F**2*(1.+(4.*lam_w*R/(pi*B*r))**2)>=(B*G/(4.*pi*r))**2]),

                            #SignomialEquality(W**2,(Wa**2+Wt**2)),
                            #SignomialEquality(1,etav*(1+eps*Wt/Wa)+eps*Wa/Wt),
                            #SignomialEquality(1,etai*(1+va/V)+vt/Ut),
                            #SignomialEquality(Q,rho*B*G*(Wa+eps*Wt)*R**2),
                            #SignomialEquality(T,rho*B*G*(Wt-eps*Wa)*R),
                    
            ]
        return constraints, state

class One_Element_Propeller(Model):
    """ Propeller Model

    Variables
    ---------
    R                               [ft]            prop radius
    W                               [lbf]           prop weight
    K           4e-4                [1/ft^2]        prop weight scaling factor
    T_m                             [lbf]           prop max static thrust
    """

    flight_model = One_Element_Propeller_Performance
    def setup(self):
        exec parse_variables(One_Element_Propeller.__doc__)

        constraints = [W >= K*T_m*R**2]

        return constraints