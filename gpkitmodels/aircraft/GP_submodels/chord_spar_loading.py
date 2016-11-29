" cap spar "
from gpkit import Model, Variable, Vectorize
from constant_taper_chord import c_bar
from beam import Beam

class ChordSparL(Model):
    "spar loading model"
    def __init__(self, static, Wcent, **kwargs):

        Nmax = Variable("N_{max}", 5, "-", "max loading")
        cbar = c_bar(0.5, static.N)
        sigmacfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")
        with Vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        beam = Beam(static.N, cbar)

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wcent/static["b"]**2),
            Mr == (beam["\\bar{M}"][:-1]*Wcent*Nmax*static["b"]/4),
            sigmacfrp >= Mr/static["S_y"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            ]

        Model.__init__(self, None, [beam, constraints], **kwargs)
