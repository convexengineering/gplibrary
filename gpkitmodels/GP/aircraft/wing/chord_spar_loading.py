" cap spar "
from numpy import flip
from gpkit import Model, Variable, Vectorize
from gpkitmodels.GP.beam.beam import Beam

#pylint: disable=invalid-name

class ChordSparL(Model):
    "spar loading model"

    def new_qbarFun(self, c):
        " define qbar model for chord loading "
        return [f(c) for f in self.static.substitutions["\\bar{c}"]]

    def new_SbarFun(self, c):
        " define Sbar model for chord loading "
        Sb = [1e-10]*self.static.N
        for i in flip(range(self.static.N-1), 0):
            Sb[i] = (Sb[i+1] + self.static.substitutions["d\\eta"][i](c)
                     * (self.new_qbarFun(c)[i]
                        + self.new_qbarFun(c)[i+1])/2)
        return Sb

    def setup(self, static, Wcent):

        self.static = static
        Nmax = Variable("N_{max}", 5, "-", "max loading")

        sigmacfrp = Variable("\\sigma_{CFRP}", 1700e6, "Pa", "CFRP max stress")
        taucfrp = Variable("\\tau_{CFRP}", 450e6, "Pa", "CFRP fabric stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")

        with Vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")


        Beam.qbarFun = self.new_qbarFun
        Beam.SbarFun = self.new_SbarFun
        beam = Beam(static.N)

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            beam["dx"] == static["d\\eta"],
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wcent/static["b"]**2),
            Mr == (beam["\\bar{M}"][:-1]*Wcent*Nmax*static["b"]/4),
            sigmacfrp >= Mr/static["S_y"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            taucfrp >= (beam["\\bar{S}"][-1]*Wcent*Nmax/4/static["t_{shear}"]
                        / static["c_{ave}"]/static["\\tau"])
            ]

        return beam, constraints
