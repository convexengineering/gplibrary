" tail boom model "
import numpy as np
from gpkit import Variable, Model, parse_variables

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

    def setup(self):
        exec parse_variables(TailBoom.__doc__)

        return [I0 <= np.pi*t0*d0**3/8.0,
                W/mfac >= np.pi*g*rhocfrp*d0*l*t0*kfac,
                t0 >= tmin,
                J <= np.pi/8.0*d0**3*t0,
                S == l*np.pi*d0,
               ]

    def flight_model(self, state):
        return TailBoomAero(self, state)

    def horizontalbending(self, htail, state):
        return HorizontalBoomBending(self, htail, state)

    def verticalbending(self, vtail, state):
        return VerticalBoomBending(self, vtail, state)

    def verticaltorsion(self, vtail, state):
        return VerticalBoomTorsion(self, vtail, state)

class TailBoomAero(Model):
    """horizontal tail aero model

    Upper Unbounded
    ---------------
    Re, l, Cf

    """
    def setup(self, static, state):

        Cf = self.Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = self.Re = Variable("Re", "-", "fuselage reynolds number")

        self.l = static["l"]

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["l"]/state["\\mu"]),
            Cf >= 0.455/Re**0.3,
            ]

        return constraints

class TailBoomState(Model):
    "tail boom design state"
    def setup(self):

        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 40, "m/s", "never exceed vehicle speed")

        constraints = [rhosl == rhosl,
                       Vne == Vne]

        return constraints

class VerticalBoomTorsion(Model):
    """tail boom torsion case

    Upper Unbounded
    ---------------
    J

    Lower Unbounded
    ---------------
    d0, b, S

    """
    def setup(self, tailboom, vtail, state):

        T = Variable("T", "N*m", "vertical tail moment")
        taucfrp = Variable("\\tau_{CFRP}", 210, "MPa", "torsional stress limit")

        J = self.J = tailboom.J
        d0 = self.d0 = tailboom.d0
        b = self.b = vtail.planform.b
        S = self.S = vtail.planform.S
        CLmax = vtail.planform.CLmax

        constraints = [
            T >= 0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*S*CLmax*b,
            taucfrp >= T*d0/2/J
            ]

        return constraints

class VerticalBoomBending(Model):
    """tail boom bending loading case

    Upper Unbounded
    ---------------
    I0

    Lower Unbounded
    ---------------
    S, l

    """
    def setup(self, tailboom, vtail, state):

        F = Variable("F", "N", "vertical tail force")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.1, "-",
                         "max tail boom deflection angle")

        I0 = self.I0 = tailboom.I0
        l = self.l = tailboom.l
        S = self.S = vtail.planform.S
        CLmax = vtail.planform.CLmax

        constraints = [
            F >= 0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*S*CLmax,
            th >= F*l**2/tailboom["E"]/I0 * (1+tailboom["k"])/2,
            th <= thmax,
            ]

        return constraints

class HorizontalBoomBending(Model):
    """tail boom bending loading case

    Upper Unbounded
    ---------------
    I0

    Lower Unbounded
    ---------------
    S, l

    """
    def setup(self, tailboom, htail, state):

        F = Variable("F", "N", "horizontal tail force")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.1, "-",
                         "max tail boom deflection angle")

        I0 = self.I0 = tailboom.I0
        l = self.l = tailboom.l
        S = self.S = htail.planform.S
        CLmax = htail.planform.CLmax

        constraints = [
            F >= 0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*S*CLmax,
            th >= F*l**2/tailboom["E"]/I0 * (1+tailboom["k"])/2,
            th <= thmax,
            ]

        return constraints
