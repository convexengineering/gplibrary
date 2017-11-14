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
    sigmacfrp       1700e6  [Pa]    CFRP max stress
    taucfrp         450e6   [Pa]    CFRP fabric stress
    kappa           0.2     [-]     max tip deflection ratio
    W                       [lbf]   loading weight

    Variables of length wing.N-1
    ----------------------------
    Mr                      [N*m]   wing section root moment

    Lower Unbounded
    ---------------
    W

    LaTex Strings
    -------------
    Nmax                N_{\\mathrm{max}}
    sigmacfrp           \\sigma_{\\mathrm{CFRP}}
    taucfrp             \\tau_{\\mathrm{CFRP}}
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

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            self.beam["dx"] == self.wing.planform.deta,
            self.beam["\\bar{EI}"] <= (8*self.wing.spar.E*self.wing.spar.I/Nmax
                                       / W/self.wing.planform.b**2),
            Mr >= (self.beam["\\bar{M}"][:-1]*W*Nmax
                   * self.wing.planform.b/4),
            sigmacfrp >= Mr/self.wing.spar.Sy,
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            taucfrp >= (self.beam["\\bar{S}"][-1]*W*Nmax/4
                        / self.wing.spar.tshear/self.wing.planform.cave
                        / self.wing.planform.tau)
            ]

        return self.beam, constraints
