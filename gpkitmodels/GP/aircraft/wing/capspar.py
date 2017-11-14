" cap spar "
from gpkit import Model, parse_variables
from .sparloading import SparLoading
from .gustloading import GustL

#pylint: disable=exec-used, undefined-variable, unused-argument, invalid-name

class CapSpar(Model):
    """ Cap Spar Model

    Scalar Variables
    ----------------
    rhocfrp         1.6     [g/cm^3]    density of CFRP
    E               2e7     [psi]       Young modulus of CFRP
    W                       [lbf]       spar weight
    wlim            0.15    [-]         spar width to chord ratio
    tshearmin       0.012   [in]        min gauge shear web thickness
    g               9.81    [m/s^2]     gravitational acceleration
    mfac            0.97    [-]         curvature knockdown factor
    rhocore         0.036   [g/cm^3]    foam core density

    Variables of length N-1
    -----------------------
    hin                     [in]        height between caps
    I                       [m^4]       spar x moment of inertia
    Sy                      [m^3]       section modulus
    dm                      [kg]        segment spar mass
    w                       [in]        spar width
    t                       [in]        spar cap thickness
    tshear                  [in]        shear web thickness

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    Sy

    LaTex Strings
    -------------
    rhocfrp                 \\rho_{\\mathrm{CFRP}}
    wlim                    w_{\\mathrm{lim}}
    tshearmin               t_{\\mathrm{shear-min}}
    mfac                    m_{\\mathrm{fac}}
    rhocore                 \\rho_{\\mathrm{core}}
    hin                     h_{\\mathrm{in}_i}
    I                       I_i
    Sy                      S_{y_i}
    dm                      \\Delta{m}
    w                       w_i
    t                       t_i
    tshear                  t_{\\mathrm{shear}_i}

    """
    loading = SparLoading
    gustloading = GustL

    def setup(self, N, surface):
        exec parse_variables(CapSpar.__doc__)

        return [I/mfac <= 2*w*t*(hin/2)**2,
                dm >= (rhocfrp*(2*w*t + 2*tshear*(hin + 2*t))
                       + rhocore*w*hin)*surface.b/2*surface.deta,
                W >= 2*dm.sum()*g,
                w <= wlim*surface.cave,
                surface.cave*surface.tau >= hin + 2*t,
                Sy*(hin/2 + t) <= I,
                tshear >= tshearmin
               ]

