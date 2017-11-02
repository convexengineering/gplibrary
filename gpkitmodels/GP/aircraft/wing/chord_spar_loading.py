" cap spar "
from numpy import flip
from gpkit import Model, Variable, Vectorize
from gpkitmodels.GP.beam.beam import Beam

#pylint: disable=invalid-name

class ChordSparL(Model):
    "spar loading model"
    def setup(self, static, Wcent):

        Nmax = Variable("N_{max}", 5, "-", "max loading")

        sigmacfrp = Variable("\\sigma_{CFRP}", 1700e6, "Pa", "CFRP max stress")
        taucfrp = Variable("\\tau_{CFRP}", 450e6, "Pa", "CFRP fabric stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")

        with Vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        def new_qbarFun(_, c):
            " define qbar model for chord loading "
            return [f(c) for f in static.substitutions["\\bar{c}"]]
        def new_SbarFun(bmodel, c):
            " define Sbar model for chord loading "
            Sb = [1e-10]*static.N
            for i in flip(range(static.N-1), 0):
                Sb[i] = (Sb[i+1] + static.substitutions["d\\eta"][i](c)
                         * (new_qbarFun(bmodel, c)[i]
                            + new_qbarFun(bmodel, c)[i+1])/2)
            return Sb

        Beam.qbarFun = new_qbarFun
        Beam.SbarFun = new_SbarFun
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
