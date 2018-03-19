" tail boom model "
from numpy import pi
from gpkit import Model, parse_variables, Variable, VectorVariable, units
from .tube_spar import TubeSpar
from gpkitmodels.GP.beam.beam import Beam
from gpkitmodels import g

#pylint: disable=exec-used, undefined-variable, invalid-name
#pylint: disable=attribute-defined-outside-init

class TailBoomAero(Model):
    """ Tail Boom Aero Model

    Variables
    ---------
    Cf          [-]     tail boom skin friction coefficient
    Re          [-]     tail boom reynolds number

    Upper Bounded by state
    ----------------------
    \\rho, \\mu

    Lower Bounded by state
    ----------------------
    \\rho, \\mu, V

    Upper Unbounded
    ---------------
    Re, Cf

    Lower Unbounded
    ---------------
    l

    LaTex Strings
    -------------
    Cf      C_f

    """
    def setup(self, static, state):
        exec parse_variables(TailBoomAero.__doc__)

        l = self.l = static.l
        self.state = state
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

class TailBoomBending(Model):
    """ Tail Boom Bending

    Variables
    ---------
    F                       [N]     tail force
    th                      [-]     tail boom deflection angle
    kappa           0.1     [-]     max tail boom deflection
    Nsafety         1.0     [-]     safety load factor

    Variables of length N-1
    -----------------------
    Mr                      [N*m]   section root moment

    Lower Bounded by htail
    ----------------------
    CLmax

    Lower Bounded by tailboom
    ----------------------
    deta

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
        N = tailboom.N
        exec parse_variables(TailBoomBending.__doc__)

        Beam.qbarFun = [1e-10]*N
        Beam.SbarFun = [1.]*N
        beam = Beam(N)

        I = self.I = tailboom.I
        self.I0 = I[0]
        l = self.l = tailboom.l
        S = self.S = htail.planform.S
        E = self.E = tailboom.material.E
        Sy = self.Sy = tailboom.Sy
        qne = self.qne = state.qne
        self.htail = htail
        CLmax = htail.planform.CLmax
        self.tailboom = tailboom
        deta = tailboom.deta
        sigma = tailboom.material.sigma

        constraints = [beam["dx"] == deta,
                       F >= qne*S,
                       beam["\\bar{EI}"] <= E*I/F/l**2/2,
                       Mr >= beam["\\bar{M}"][:-1]*F*l,
                       sigma >= Mr/Sy,
                       th == beam["\\theta"][-1],
                       beam["\\bar{\\delta}"][-1]*CLmax*Nsafety <= kappa]

        if hasattr(tailboom, "J"):
            constraints.append(tailboom.J >= 1e-50*units("m^4"))

        return constraints, beam

class TailBoom(TubeSpar):
    """ Tail Boom Model

    Variables
    ---------
    l                           [ft]        tail boom length
    S                           [ft^2]      tail boom surface area
    b                           [ft]        twice tail boom length
    deta          1./(N-1)      [-]         normalized segment length
    tau           1.0           [-]         thickness to width ratio
    rhoA          0.15          [kg/m^2]    total aerial density

    Variables of length N-1
    -----------------------
    cave                        [in]        average segment width

    """

    flight_model = TailBoomAero
    tailLoad = TailBoomBending
    secondaryWeight = None

    def setup(self, N=5):
        self.N = N
        exec parse_variables(TailBoom.__doc__)
        self.spar = super(TailBoom, self).setup(N, self)

        if self.secondaryWeight:
            self.weight.right += rhoA*g*S

        d0 = self.d0 = self.d[0]

        return self.spar, [S == l*pi*d0, b == 2*l]
