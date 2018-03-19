" wing skin "
from gpkit import Model, parse_variables
from gpkitmodels.GP.materials import cfrpfabric
from gpkitmodels import g

class WingSkin(Model):
    """ Wing Skin model

    Variables
    ---------
    W                           [lbf]           wing skin weight
    t                           [in]            wing skin thickness
    Jtbar           0.01114     [1/mm]          torsional moment of inertia
    Cmw             0.121       [-]             negative wing moment coeff
    rhosl           1.225       [kg/m^3]        sea level air density
    Vne             45          [m/s]           never exceed vehicle speed

    Upper Unbounded
    ---------------
    W, surface.croot

    Lower Unbounded
    ---------------
    surface.S

    LaTex Strings
    -------------
    W       W_{\\mathrm{skin}}
    t       t_{\\mathrm{skin}}
    Jtbar   \\bar{J/t}
    Cmw     C_{m_w}
    rhosl   \\rho_{\\mathrm{SL}}
    Vne     V_{\\mathrm{NE}}

    """
    material = cfrpfabric

    def setup(self, surface):
        self.surface = surface
        exec parse_variables(WingSkin.__doc__)

        croot = surface.croot
        S = surface.S
        rho = self.material.rho
        tau = self.material.tau
        tmin = self.material.tmin

        return [W >= rho*surface.S*2*t*g,
                t >= tmin,
                tau >= 1/Jtbar/croot**2/t*Cmw*S*rhosl*Vne**2]

class WingSecondStruct(Model):
    """ Wing Skin model

    Variables
    ---------
    W                           [lbf]           wing skin weight
    rhoA            0.35        [kg/m^2]        total aerial density

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    surface.S

    LaTex Strings
    -------------
    W       W_{\\mathrm{skin}}
    rhoA    \\rho_{A}

    """
    def setup(self, surface):
        self.surface = surface
        exec parse_variables(WingSecondStruct.__doc__)

        S = surface.S

        return [W >= rhoA*S*g]
