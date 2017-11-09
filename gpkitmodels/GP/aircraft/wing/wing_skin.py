" wing skin "
from gpkit import Model, Variable, parse_variables

class WingSkin(Model):
    """ Wing Skin model

    Variables
    ---------
    rhocfrp         1.6         [g/cm^3]        density of CFRP
    W                           [lbf]           wing skin weight
    g               9.81        [m/s^2]         gravitational acceleration
    t                           [in]            wing skin thickness
    tmin            0.012       [in]            wing skin min gauge
    Jtbar           0.01114     [1/mm]          torsional moment of inertia
    taucfrp         570         [MPa]           torsional stress limit
    Cmw             0.121       [-]             negative wing moment coeff
    rhosl           1.225       [kg/m^3]        sea level air density
    Vne             45          [m/s]           never exceed vehicle speed

    Upper Unbounded
    ---------------
    W

    LaTex Strings
    -------------
    W       W_{\\mathrm{skin}}
    g       g
    t       t_{\\mathrm{skin}}
    t       t_{\\mathrm{min}}
    Jtbar   \\bar{J/t}
    taucfrp \\tau_{\\mathrm{CFRP}}
    Cmw     C_{m_w}
    rhosl   \\rho_{\\mathrm{SL}}
    Vne     V_{\\mathrm{NE}}

    """
    def setup(self, surface):
        exec parse_variables(WingSkin.__doc__)

        self.loading = WingSkinL
        return [W >= rhocfrp*surface["S"]*2*t*g,
                t >= tmin,
                taucfrp >= (1/Jtbar/(surface["c_{root}"])**2/t*Cmw
                            * surface["S"]*rhosl*Vne**2)
                ]

class WingSkinL(Model):
    "wing skin loading model for torsional loads in skin"
    def setup(self, static):

        taucfrp = Variable("\\tau_{CFRP}", 570, "MPa", "torsional stress limit")
        Cmw = Variable("C_{m_w}", 0.121, "-", "negative wing moment coefficent")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 45, "m/s", "never exceed vehicle speed")

        return []
