" cylindrical fuselage.py "
import numpy as np
from gpkit import Variable, Model, parse_variables, SignomialsEnabled, SignomialEquality
from fuselage_skin import FuselageSkin
from gpkit.constraints.tight import Tight as TCS


class FuselageLoading(Model):
    "fuselage loading cases"
    def setup(self, fuselage, Wcent):

        loading = [fuselage.skin.loading(Wcent)]

        return loading

class FuselageAero(Model):
    """fuselage drag model
   
    Variables
    ---------
    Cf                      [-]             fuselage skin friction coefficient
    Re                      [-]             fuselage Reynolds number
    Cd                      [-]             fuselage drag coefficient
    mfac            1.1     [-]             fuselage drag scaling factor

    """
    #Reref                   [-]             reference Reynolds number
    #Cfref                   [-]             reference skin friction coefficient
    def setup(self, static, state):
        exec parse_variables(FuselageAero.__doc__)

        rho = state.rho
        l = static.l
        mu = state.mu
        fbody = static.fbody
        fbulk = static.fbulk
        fnose = static.fnose
        V = state.V
        k =  static.k

        constraints = [
            Re == V*rho*l/mu,
            Cf >= 0.455/Re**0.3,
            #Cfref == 0.455/Reref**0.3,
            Cf >= 0.455/Re**0.3,
            Cd/mfac >= Cf*k
            #(Cd/mfac)**0.996232 >= Cf/Cfref*(
            #    0.00243049*fbody**0.033607
            #    * fnose**1.21682 * fbulk**0.306251
            #    + 0.00255095*fbody**-0.0316887
            #    * fnose**-0.585489 * fbulk**1.15394
            #    + 0.0436011 * fbody**0.0545722
            #    * fnose**0.258228 * fbulk**-1.42664
            #    + 0.00970479 * fbody**0.8661
            #    * fnose**-0.209136 * fbulk**-0.156166)
            ]

        return constraints
class Fuselage(Model):
    """The thing that carries the motor, and payload

    Variables
    ---------
    R                           [in]            fuselage radius
    l               16          [in]            fuselage length
    Sw                          [ft^2]          total wetted fuselage area
    Sw_c                        [ft^2]          center wetted fuselage area
    Sw_n                        [ft^2]          nose wetted fuselage area
    Sw_t                        [ft^2]          tail wetted fuselage area
    W                           [gf]           fuselage weight
    mfac            2.0         [-]             fuselage weight margin factor
    fbody                       [-]             fuselage body fineness ratio
    fnose                       [-]             fuselage nose fineness ratio
    fbulk           1.0         [-]             fuselage bulk fineness ratio
    f                           [-]             overall fineness ratio
    k                           [-]             fuselage form factor
    Vol                         [ft^3]          fuselage volume
    Vol_nose                    [in^3]          nose volume
    Vol_tail                    [in^3]          tail volume
    Vol_center                  [in^3]          tail volume
    rhocfrp 1.6                 [g/cm^3]        density of CFRP
    t                           [in]            fuselage skin thickness
    nply            2           [-]             number of plys
    lbody                       [in]            center body length
    lnose                       [in]            nose cone length
    ltail                       [in]            tail cone length
    Vol_payload                 [in^3]          Req'd payload volume
    Rm                          [in]            Motor radius
    t_ins           5           [mm]            Insulation thickness
    t_TEG           6           [mm]            TEG thickness
    t_margin        5           [mm]            Thickness margin
    fnose_max       2           [-]             max nose fineness ratio
    """

    flight_model = FuselageAero
    loading = FuselageLoading

    def setup(self):
        exec parse_variables(Fuselage.__doc__)

        self.skin = FuselageSkin(Sw, R, lbody)
        self.components = [self.skin]

        constraints = [
            fbody == lbody/R,
            fnose == lnose/R,
            fbulk == ltail/R,
            Sw >= Sw_c + Sw_n + Sw_t,
            Sw_c >= 2*np.pi*R*lbody,
            Sw_n**(8./5.) >= (
                (2*np.pi*R**2)**(8./5.)*(1./3. + 2./3.*(fnose)**(8./5.))),

            Sw_t >= R**2*(0.012322*fbulk**2 + 1.524925*fbulk + 0.502498),
            
            Vol_center == np.pi*R**2*lbody,
            Vol_nose == 4./3.*np.pi*R**2.0*lnose,
            #Vol_tail == 4./3.*np.pi*R**3.0,
            f == l/R/2.,
            k >= 1 + 60/f**3 + f/400,
            #l <= 3.*R*(fbody*fnose*fbulk)**(1./3),
            Sw >= np.pi*R**2,
            #Vol >= Vol_center + Vol_nose + Vol_tail,
            Vol_nose >= Vol_payload,
            W/mfac >= self.skin["W"],
            fnose <= fnose_max,
            2.*R >= 2.*Rm+t_ins+t_TEG+t_margin,
            #TCS([l >= lnose + ltail + lbody]),
            ]
        with SignomialsEnabled():
            constraints.extend([SignomialEquality(l,lnose + ltail + lbody),
                                ])

        return self.components, constraints



