" cap spar "
from numpy import flip
from gpkit import Model, Variable, Vectorize
from gpkitmodels.GP.beam.beam import Beam

#pylint: disable=invalid-name

class SparLoading(Model):
    "spar loading model"

    def new_qbarFun(self, c):
        " define qbar model for chord loading "
        barc = self.static.planform.cbar
        return [f(c) for f in self.static.substitutions[barc]]

    def setup(self, static):

        self.static = static
        Nmax = Variable("N_{max}", 5, "-", "max loading")

        sigmacfrp = Variable("\\sigma_{CFRP}", 1700e6, "Pa", "CFRP max stress")
        taucfrp = Variable("\\tau_{CFRP}", 450e6, "Pa", "CFRP fabric stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")
        self.W = Variable("W", "lbf", "loading weight")

        with Vectorize(self.static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")


        Beam.qbarFun = self.new_qbarFun
        self.beam = Beam(self.static.N)

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            self.beam["dx"] == self.static.planform.deta,
            self.beam["\\bar{EI}"] <= (8*self.static["E"]*self.static["I"]/Nmax
                                       / self.W/self.static.planform.b**2),
            Mr == (self.beam["\\bar{M}"][:-1]*self.W*Nmax
                   * self.static.planform.b/4),
            sigmacfrp >= Mr/self.static["S_y"],
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            taucfrp >= (self.beam["\\bar{S}"][-1]*self.W*Nmax/4
                        / self.static["t_{shear}"]/self.static.planform.cave
                        / self.static.planform.tau)
            ]

        return self.beam, constraints
