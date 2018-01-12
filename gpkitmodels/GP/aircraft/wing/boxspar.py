" box spar "
from gpkit import Model, parse_variables
from .sparloading import SparLoading
from .gustloading import GustL
from gpkitmodels.GP.materials import cfrpud, cfrpfabric, foamhd
from gpkitmodels import g

#pylint: disable=exec-used, undefined-variable, unused-argument, invalid-name

class BoxSpar(Model):
    """ Box Spar Model

    Scalar Variables
    ----------------
    W                       [lbf]       spar weight
    wlim            0.15    [-]         spar width to chord ratio
    mfac            0.97    [-]         curvature knockdown factor
    tcoret          0.02    [-]         core to thickness ratio

    Variables of length N-1
    -----------------------
    hin                     [in]        height between caps
    I                       [m^4]       spar x moment of inertia
    Sy                      [m^3]       section modulus
    dm                      [kg]        segment spar mass
    w                       [in]        spar width
    d                       [in]        cross sectional diameter
    t                       [in]        spar cap thickness
    tshear                  [in]        shear web thickness
    tcore                   [in]        core thickness

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    Sy, b

    LaTex Strings
    -------------
    wlim                    w_{\\mathrm{lim}}
    mfac                    m_{\\mathrm{fac}}
    hin                     h_{\\mathrm{in}_i}
    I                       I_i
    Sy                      S_{y_i}
    dm                      \\Delta{m}
    w                       w_i
    t                       t_i
    tshear                  t_{\\mathrm{shear}_i}
    tcoret                  (t_{\\mathrm{core}}/t)

    """
    loading = SparLoading
    gustloading = GustL
    material = cfrpud
    shearMaterial = cfrpfabric
    coreMaterial = foamhd

    def setup(self, N, surface):
        exec parse_variables(BoxSpar.__doc__)

        b = self.b = surface.b
        cave = surface.cave
        tau = surface.tau
        deta = surface.deta
        rho = self.material.rho
        rhoshear = self.shearMaterial.rho
        rhocore = self.coreMaterial.rho
        tshearmin = self.shearMaterial.tmin

        self.weight = W >= 2*dm.sum()*g

        return [I/mfac <= w*t*hin**2,
                dm >= (rho*(4*w*t) + 2*tshear*rhoshear*(hin + 2*tcore + 4*t)
                       + rhocore*w*tcore*2)*b/2*deta,
                w <= wlim*cave,
                cave*tau >= hin + 4*t + 2*tcore,
                self.weight,
                Sy*(hin/2 + 2*t + tcore) <= I,
                tshear >= tshearmin,
                tcore >= tcoret*cave*tau,
                d == w,
               ]
