" tail boom model "
from numpy import pi
from gpkit import Model, parse_variables
from .tube_spar import TubeSpar

#pylint: disable=exec-used, undefined-variable, invalid-name
#pylint: disable=attribute-defined-outside-init

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
        rho = self.rho = state.rho
        V = self.V = state.V
        mu = self.mu = state.mu

        return [Re == V*rho*l/mu,
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
        d0 = self.d0 = tailboom.d
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

        I0 = self.I0 = tailboom.I
        l = self.l = tailboom.l
        S = self.S = vtail.planform.S
        E = self.E = tailboom.material.E
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

        I0 = self.I0 = tailboom.I
        l = self.l = tailboom.l
        S = self.S = htail.planform.S
        E = self.E = tailboom.material.E
        k = self.k = tailboom.k
        rhosl = self.rhosl = state.rhosl
        Vne = self.Vne = state.Vne
        CLmax = htail.planform.CLmax

        return [F >= 0.5*rhosl*Vne**2*S,
                th >= F*l**2/E/I0*(1+k)/2,
                th*CLmax <= thmax,
               ]

class TailBoom(TubeSpar):
    """ Tail Boom Model

    Variables
    ---------
    l                           [ft]        tail boom length
    S                           [ft^2]      tail boom surface area
    deta          1./(N-1)      [-]         normalized segment length

    """

    flight_model = TailBoomAero
    hbending = HorizontalBoomBending
    vbending = VerticalBoomBending
    vtorsion = VerticalBoomTorsion

    def setup(self, N=2):
        exec parse_variables(TailBoom.__doc__)

        self.spar = TubeSpar.setup(self, N, self)

        d0 = self.d0 = self.d[0]

        return self.spar, S == l*pi*d0


