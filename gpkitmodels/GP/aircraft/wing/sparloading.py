" cap spar "
from numpy import flip
from gpkit import Model, Variable, Vectorize
from gpkitmodels.GP.beam.beam import Beam

#pylint: disable=invalid-name

class ChordSparL(Model):
    " place holder so dependent repos don't break "
    def setup(self, static, Wcent):

        sploading = SparLoading(static, Wcent)

        return sploading


class SparLoading(Model):
    "spar loading model"

    def new_qbarFun(self, c):
        " define qbar model for chord loading "
        barc = self.static["\\bar{c}"]
        return [f(c) for f in self.static.substitutions[barc]]

    def setup(self, static, Wcent=None):

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
            self.beam["dx"] == self.static["d\\eta"],
            self.beam["\\bar{EI}"] <= (8*self.static["E"]*self.static["I"]/Nmax
                                       / self.W/self.static["b"]**2),
            Mr == (self.beam["\\bar{M}"][:-1]*self.W*Nmax*self.static["b"]/4),
            sigmacfrp >= Mr/self.static["S_y"],
            self.beam["\\bar{\\delta}"][-1] <= kappa,
            taucfrp >= (self.beam["\\bar{S}"][-1]*self.W*Nmax/4
                        / self.static["t_{shear}"]/self.static["c_{ave}"]
                        / self.static["\\tau"])
            ]

        if Wcent:
            constraints.extend([self.W == Wcent])

        return self.beam, constraints
