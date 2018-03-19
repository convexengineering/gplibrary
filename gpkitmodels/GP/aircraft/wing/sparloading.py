" spar loading "
from gpkit import Model, parse_variables
from gpkitmodels.GP.beam.beam import Beam
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoomState
from numpy import pi

#pylint: disable=no-member, unused-argument, exec-used, invalid-name
#pylint: disable=undefined-variable, attribute-defined-outside-init

class SparLoading(Model):
    """ Spar Loading Model

    Variables
    ---------
    Nmax            5              [-]     max loading
    Nsafety         1.0            [-]     safety load factor
    kappa           0.2            [-]     max tip deflection ratio
    W                              [lbf]   loading weight
    twmax           30.*pi/180     [-]     max tip twist

    Variables of length wing.N-1
    ----------------------------
    Mr                      [N*m]   wing section root moment
    M                       [N*m]   local moment
    theta                   [-]     twist deflection

    Upper Unbounded
    ---------------
    I, Sy, J (if wingSparJ)
    theta (if not wingSparJ), M (if not wingSparJ)

    Lower Unbounded
    ---------------
    b, W, cave, qne (if wingSparJ)
    theta (if not wingSparJ), M (if not wingSparJ)

    LaTex Strings
    -------------
    Nmax                N_{\\mathrm{max}}
    kappa               \\kappa
    Mr                  M_r

    """
    def new_qbarFun(self, c):
        " define qbar model for chord loading "
        barc = self.wing.planform.cbar
        return [f(c) for f in self.wing.substitutions[barc]]

    new_SbarFun = None

    def setup(self, wing, state):
        self.wing = wing
        exec parse_variables(SparLoading.__doc__)

        Beam.qbarFun = self.new_qbarFun
        Beam.SbarFun = self.new_SbarFun
        self.beam = Beam(self.wing.N)

        b = self.b = self.wing.planform.b
        I = self.I = self.wing.spar.I
        Sy = self.Sy = self.wing.spar.Sy
        cave = self.cave = self.wing.planform.cave
        E = self.wing.spar.material.E
        sigma = self.wing.spar.material.sigma
        deta = self.wing.planform.deta

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            self.beam["dx"] == deta,
            self.beam["\\bar{EI}"] <= 8*E*I/Nmax/Nsafety/W/b**2,
            Mr >= self.beam["\\bar{M}"][:-1]*W*Nmax*Nsafety*b/4,
            sigma >= Mr/Sy,
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            ]

        self.wingSparJ = hasattr(self.wing.spar, "J")

        if self.wingSparJ:
            qne = self.qne = state.qne
            J = self.J = self.wing.spar.J
            G = self.wing.spar.shearMaterial.G
            cm = self.wing.planform.CM
            constraints.extend([
                M >= cm*cave**2*qne*deta*b/2*Nsafety,
                theta[0] >= M[0]/G/J[0]*deta[0]*b/2,
                theta[1:] >= theta[:-1] + M[1:]/G/J[1:]*deta[1:]*b/2,
                twmax >= theta[-1]
                ])
        return self.beam, constraints
