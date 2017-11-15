" tail boom model "
import numpy as np
from gpkit import Variable, Model, parse_variables

class TailBoomAero(Model):
    """ Tail Boom Aero Model

    Variables
    ---------
    Cf          [-]     tail boom skin friction coefficient
    Re          [-]     tail boom reynolds number

    Upper Unbounded
    ---------------
    Re, l, Cf

    LaTex Strings
    -------------
    Cf      C_f

    """
    def setup(self, static, state):
        exec parse_variables(TailBoomAero.__doc__)

        l = self.l = static.l

        return [Re == (state["V"]*state["\\rho"]*l/state["\\mu"]),
                Cf >= 0.455/Re**0.3,
                ]

class TailBoomState(Model):
    """ Tail Boom Loading State

    Variables
    ---------
    rhosl           1.225           [kg/m^3]    air density at sea level
    Vne             40              [m/s]       never exceed vehicle speed

    LaTex Strings
    -------------
    rhosl           \\rho_{\\mathrm{sl}}
    Vne             V_{\\mathrm{NE}}

    """
    def setup(self):
        exec parse_variables(TailBoomState.__doc__)


class VerticalBoomTorsion(Model):
    """ Tail Boom Torsion from Vertical Tail

    Variables
    ---------
    T                           [N*m]       vertical tail moment
    taucfrp         210         [MPa]       torsional stress limit of carbon

    Upper Unbounded
    ---------------
    J

    Lower Unbounded
    ---------------
    d0, b, S

    LaTex Strings
    -------------
    taucfrp     \\tau_{\\mathrm{CFRP}}

    """
    def setup(self, tailboom, vtail, state):
        exec parse_variables(VerticalBoomTorsion.__doc__)

        J = self.J = tailboom.J
        d0 = self.d0 = tailboom.d0
        b = self.b = vtail.planform.b
        S = self.S = vtail.planform.S
        rhosl = self.rhosl = state.rhosl
        Vne = self.Vne = state.Vne
        CLmax = vtail.planform.CLmax

        return [T >= 0.5*rhosl*Vne**2*S*CLmax*b,
                taucfrp >= T*d0/2/J
                ]

class VerticalBoomBending(Model):
    """ Tail Boom Bending from Vtail Deflection

    Variables
    ---------
    F                       [N]     vertical tail force
    th                      [-]     tail boom deflection angle
    thmax           0.1     [-]     max tail boom deflection angle

    Upper Unbounded
    ---------------
    I0

    Lower Unbounded
    ---------------
    S, l

    LaTex Strings
    -------------
    th      \\theta
    thmax   \\theta_{\\mathrm{max}}

    """
    def setup(self, tailboom, vtail, state):
        exec parse_variables(VerticalBoomBending.__doc__)

        I0 = self.I0 = tailboom.I0
        l = self.l = tailboom.l
        S = self.S = vtail.planform.S
        E = self.E = tailboom.E
        k = self.k = tailboom.k
        rhosl = self.rhosl = state.rhosl
        Vne = self.Vne = state.Vne
        CLmax = vtail.planform.CLmax

        return [F >= 0.5*rhosl*Vne**2*S*CLmax,
                th >= F*l**2/E/I0*(1+k)/2,
                th <= thmax,
               ]

class HorizontalBoomBending(Model):
    """ Tail Boom Bending from Htail Deflection

    Variables
    ---------
    F                       [N]     horizontal tail force
    th                      [-]     tail boom deflection angle
    thmax           0.1     [-]     max tail boom deflection angle

    Upper Unbounded
    ---------------
    I0

    Lower Unbounded
    ---------------
    S, l

    LaTex Strings
    -------------
    th      \\theta
    thmax   \\theta_{\\mathrm{max}}

    """
    def setup(self, tailboom, htail, state):
        exec parse_variables(HorizontalBoomBending.__doc__)

        I0 = self.I0 = tailboom.I0
        l = self.l = tailboom.l
        S = self.S = htail.planform.S
        E = self.E = tailboom.E
        k = self.k = tailboom.k
        rhosl = self.rhosl = state.rhosl
        Vne = self.Vne = state.Vne
        CLmax = htail.planform.CLmax

        return [F >= 0.5*rhosl*Vne**2*S*CLmax,
                th >= F*l**2/E/I0*(1+k)/2,
                th <= thmax,
               ]

class TailBoom(Model):
    """ Tail Boom Model

    Variables
    ---------
    l                           [ft]        tail boom length
    E           150e9           [N/m^2]     Youngs modulus for carbon fiber
    k           0.8             [-]         tail boom taper index
    kfac        self.minusk2    [-]         (1-k/2)
    I0                          [m^4]       tail boom moment of inertia
    d0                          [in]        tail boom diameter
    t0                          [in]        tail boom thickness
    tmin        0.25            [mm]        minimum tail boom thickness
    rhocfrp     1.6             [g/cm^3]    density of CFRP
    g           9.81            [m/s^2]     gravitational acceleration
    W                           [lbf]       tail boom weight
    J                           [m^4]       tail boom polar moment of inertia
    S                           [ft^2]      tail boom surface area
    mfac        1.0             [-]         tail boom margin factor

    Upper Unbounded
    ---------------
    W

    Lower Unbounded
    ---------------
    S, J, l, I0

    LaTex Strings
    -------------
    kfac        (1-k/2)
    I0          I_0
    d0          d_0
    t0          t_0
    tmin        t_{\\mathrm{min}}
    rhocfrp     \\rho_{\\mathrm{CFRP}}
    mfac        m_{\\mathrm{fac}}

    """

    minusk2 = lambda self, c: 1-c[self.k]/2.
    flight_model = TailBoomAero
    hbending = HorizontalBoomBending
    vbending = VerticalBoomBending
    vtorsion = VerticalBoomTorsion

    def setup(self):
        exec parse_variables(TailBoom.__doc__)

        return [I0 <= np.pi*t0*d0**3/8.0,
                W/mfac >= np.pi*g*rhocfrp*d0*l*t0*kfac,
                t0 >= tmin,
                J <= np.pi/8.0*d0**3*t0,
                S == l*np.pi*d0,
               ]

