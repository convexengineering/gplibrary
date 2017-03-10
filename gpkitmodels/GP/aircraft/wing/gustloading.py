" spar loading for gust case "
from gpkit import Model, Variable, Vectorize
from constant_taper_chord import c_bar
from gpkitmodels.GP.beam.beam import Beam
from numpy import pi
import numpy as np

class GustL(Model):
    "spar loading model"
    def setup(self, static, Wcent, Wwing, V, CL):
    # def setup(self, static, Wcent, rho, V, S):

        Nmax = Variable("N_{max}", 2, "-", "load safety factor")
        cbar, eta = c_bar(0.5, static.N)
        sigmacfrp = Variable("\\sigma_{CFRP}", 570e6, "Pa", "CFRP max stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")

        with Vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        vgust = Variable("V_{gust}", 10, "m/s", "gust velocity")

        with Vectorize(static.N):
            agust = Variable("\\alpha_{gust}", "-", "gust angle of attack")
            qbar = Variable("\\bar{q}", "-", "normalized loading")
            cosminus1 = Variable("1-cos(\\eta)",
                                 np.hstack([1e-10, 1-np.cos(eta[1:]*pi/2)]),
                                 "-", "1 minus cosine factor")

        beam = Beam(static.N, qbar)

        constraints = [
            # fit for arctan from 0 to 1, RMS = 0.044
            agust == 0.905329*(cosminus1*vgust/V)**0.961857,
            qbar >= cbar*(1 + 2*pi*agust/CL*(1+Wwing/Wcent)),
            # dimensionalize moment of inertia and young's modulus
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wcent/static["b"]**2),
            Mr == (beam["\\bar{M}"][:-1]*Wcent*Nmax*static["b"]/4),
            sigmacfrp >= Mr/static["S_y"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            ]

        return beam, constraints
