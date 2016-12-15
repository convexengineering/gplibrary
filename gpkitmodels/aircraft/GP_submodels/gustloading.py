" spar loading for gust case "
from gpkit import Model, Variable, Vectorize
from constant_taper_chord import c_bar
from beam import Beam
from numpy import pi
import numpy as np

class GustL(Model):
    "spar loading model"
    def setup(self, static, Wcent, Wwing, V, CL):

        Nmax = Variable("N_{max}", 5, "-", "max loading")
        cbar, eta = c_bar(0.5, static.N)
        sigmacfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")

        with Vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        vgust = Variable("V_{gust}", 10, "m/s", "gust velocity")
        agust = Variable("\\alpha_{gust}", "-", "gust angle of attack")

        with Vectorize(static.N):
            qbar = Variable("\\bar{q}", "-", "normalized loading")
            cosminus1 = Variable("1-cos(\\eta)", (1-np.cos(eta)/2)**2, "-",
                                 "1 minus cosine factor")

        beam = Beam(static.N, qbar)

        constraints = [
            # fit for arctan from 0 to 1, RMS = 0.055
            agust == 0.874071*(vgust/V)**0.958316,
            qbar >= 2*pi/CL*(1+Wcent/Wwing)*cosminus1*agust,
            # dimensionalize moment of inertia and young's modulus
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wwing/static["b"]**2),
            Mr == (beam["\\bar{M}"][:-1]*Wwing*Nmax*static["b"]/4),
            sigmacfrp >= Mr/static["S_y"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            ]

        return beam, constraints
