" spar loading "
from gpkit import Model, parse_variables
from gpkitmodels.GP.beam.beam import Beam
from gpkitmodels.GP.aircraft.tail.tail_boom import TailBoomState

#pylint: disable=no-member, unused-argument, exec-used, invalid-name
#pylint: disable=undefined-variable, attribute-defined-outside-init

class SparLoading(Model):
    """ Spar Loading Model

    Variables
    ---------
    Nmax            5       [-]     max loading
    Nsafety         1.0     [-]     safety load factor
    kappa           0.2     [-]     max tip deflection ratio
    W                       [lbf]   loading weight

    Variables of length wing.N-1
    ----------------------------
    Mr                      [N*m]   wing section root moment
    M                       [N*m]   local moment
    theta                   [-]     twist deflection

    Upper Unbounded
    ---------------
    I, tshear, Sy, cave

    Lower Unbounded
    ---------------
    b, W

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

    def setup(self, wing):
        self.wing = wing
        exec parse_variables(SparLoading.__doc__)

        Beam.qbarFun = self.new_qbarFun
        Beam.SbarFun = self.new_SbarFun
        self.beam = Beam(self.wing.N)

        b = self.b = self.wing.planform.b
        I = self.I = self.wing.spar.I
        Sy = self.Sy = self.wing.spar.Sy
        cave = self.cave = self.wing.planform.cave
        tshear = self.tshear = self.wing.spar.tshear
        E = self.wing.spar.material.E
        sigma = self.wing.spar.material.sigma
        tau = self.wing.planform.tau
        taumat = self.wing.spar.shearMaterial.tau
        deta = self.wing.planform.deta
        cm = self.wing.planform.CM

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            self.beam["dx"] == deta,
            self.beam["\\bar{EI}"] <= 8*E*I/Nmax/Nsafety/W/b**2,
            Mr >= self.beam["\\bar{M}"][:-1]*W*Nmax*Nsafety*b/4,
            sigma >= Mr/Sy,
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            taumat >= self.beam["\\bar{S}"][-1]*W*Nmax*Nsafety/4/tshear/cave/tau
            ]

        if hasattr(self.wing.spar, "J"):
            state = TailBoomState()
            rho = state.rhosl
            V = state.Vne
            constraints.extend([
                M >= 0.5*cm*cave**2*rho*V**2*deta*b/2,
                theta[0] >= M[0]/G/J*deta[0]*b/2,
                theta[1:] >= theta[:-1] + M[1:]/G/J*deta[1:]*b/2,
                ])
        return self.beam, constraints
