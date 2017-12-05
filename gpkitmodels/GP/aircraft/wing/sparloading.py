" spar loading "
from gpkit import Model, parse_variables
from gpkitmodels.GP.beam.beam import Beam

#pylint: disable=no-member, unused-argument, exec-used, invalid-name
#pylint: disable=undefined-variable, attribute-defined-outside-init

class SparLoading(Model):
    """ Spar Loading Model

    Variables
    ---------
    Nmax            5       [-]     max loading
    kappa           0.2     [-]     max tip deflection ratio
    W                       [lbf]   loading weight

    Variables of length wing.N-1
    ----------------------------
    Mr                      [N*m]   wing section root moment

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

    def setup(self, wing):
        self.wing = wing
        exec parse_variables(SparLoading.__doc__)

        Beam.qbarFun = self.new_qbarFun
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

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            self.beam["dx"] == self.wing.planform.deta,
            self.beam["\\bar{EI}"] <= 8*E*I/Nmax/W/b**2,
            Mr >= self.beam["\\bar{M}"][:-1]*W*Nmax*b/4,
            sigma >= Mr/Sy,
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            taumat >= self.beam["\\bar{S}"][-1]*W*Nmax/4/tshear/cave/tau
            ]

        return self.beam, constraints
