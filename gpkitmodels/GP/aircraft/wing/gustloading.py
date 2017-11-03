" spar loading for gust case "
import os
from numpy import pi, hstack, cos
import pandas as pd
from gpkit import Model, Variable, Vectorize
from gpkitmodels.GP.beam.beam import Beam
from gpfit.fit_constraintset import FitCS

#pylint: disable=invalid-name, no-member, too-many-locals

class GustL(Model):
    "spar loading model"
    def setup(self, static, Wcent, Wwing, V, CL):

        Nmax = Variable("N_{max}", 2, "-", "load safety factor")

        sigmacfrp = Variable("\\sigma_{CFRP}", 1700e6, "Pa", "CFRP max stress")
        taucfrp = Variable("\\tau_{CFRP}", 450e6, "Pa", "CFRP fabric stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")

        with Vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        vgust = Variable("V_{gust}", 10, "m/s", "gust velocity")

        with Vectorize(static.N):
            agust = Variable("\\alpha_{gust}", "-", "gust angle of attack")
            return_cosm1 = lambda c: hstack(
                [1e-10, 1-cos(c[static["\\eta"]][1:]*pi/2)])
            cosminus1 = Variable("1-cos(\\eta)", return_cosm1,
                                 "-", "1 minus cosine factor")

        Beam.qbarFun = None
        Beam.SbarFun = None
        beam = Beam(static.N)
        path = os.path.abspath(__file__).replace(os.path.basename(__file__), "")
        df = pd.read_csv(path + os.sep + "arctan_fit.csv").to_dict(
            orient="records")[0]

        constraints = [
            # fit for arctan from 0 to 1, RMS = 0.044
            FitCS(df, agust, [cosminus1*vgust/V]),
            beam["\\bar{q}"] >= static["\\bar{c}"]*(
                1 + 2*pi*agust/CL*(1+Wwing/Wcent)),
            beam["dx"] == static["d\\eta"],
            # dimensionalize moment of inertia and young's modulus
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wcent/static["b"]**2),
            Mr == (beam["\\bar{M}"][:-1]*Wcent*Nmax*static["b"]/4),
            sigmacfrp >= Mr/static["S_y"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            taucfrp >= (beam["\\bar{S}"][-1]*Wcent*Nmax/4/static["t_{shear}"]
                        / static["c_{ave}"]/static["\\tau"])
            ]

        return beam, constraints
