" cap spar "
from numpy import flip
from gpkit import Model, Variable, Vectorize
from gpkitmodels.GP.beam.beam import Beam

#pylint: disable=invalid-name

class SparLoading(Model):
    "spar loading model"

    def new_qbarFun(self, c):
        " define qbar model for chord loading "
        barc = self.wing.planform.cbar
        return [f(c) for f in self.wing.substitutions[barc]]

    def setup(self, wing):

        self.wing = wing
        Nmax = Variable("N_{max}", 5, "-", "max loading")

        sigmacfrp = Variable("\\sigma_{CFRP}", 1700e6, "Pa", "CFRP max stress")
        taucfrp = Variable("\\tau_{CFRP}", 450e6, "Pa", "CFRP fabric stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")
        self.W = Variable("W", "lbf", "loading weight")

        with Vectorize(self.wing.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        Beam.qbarFun = self.new_qbarFun
        self.beam = Beam(self.wing.N)

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            self.beam["dx"] == self.wing.planform.deta,
            self.beam["\\bar{EI}"] <= (8*self.wing.spar.E*self.wing.spar.I/Nmax
                                       / self.W/self.wing.planform.b**2),
            Mr >= (self.beam["\\bar{M}"][:-1]*self.W*Nmax
                   * self.wing.planform.b/4),
            sigmacfrp >= Mr/self.wing.spar.Sy,
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            taucfrp >= (self.beam["\\bar{S}"][-1]*self.W*Nmax/4
                        / self.wing.spar.tshear/self.wing.planform.cave
                        / self.wing.planform.tau)
            ]

        return self.beam, constraints
