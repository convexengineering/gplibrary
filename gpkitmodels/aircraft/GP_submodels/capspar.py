" cap spar "
from gpkit import Model, Variable, Vectorize
from chord_spar_loading import ChordSparL

class CapSpar(Model):
    "cap spar model"
    def setup(self, b, cave, tau, N=5):
        self.N = N

        # phyiscal properties
        rhocfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        E = Variable("E", 2e7, "psi", "Youngs modulus of CFRP")

        with Vectorize(self.N-1):
            t = Variable("t", "in", "spar cap thickness")
            hin = Variable("h_{in}", "in", "inner spar height")
            w = Variable("w", "in", "spar width")
            I = Variable("I", "m^4", "spar x moment of inertia")
            Sy = Variable("S_y", "m**3", "section modulus")
            dm = Variable("dm", "kg", "segment spar mass")

        W = Variable("W", "lbf", "spar weight")
        w_lim = Variable("w_{lim}", 0.15, "-", "spar width to chord ratio")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [I <= 2*w*t*(hin/2)**2,
                       dm >= rhocfrp*w*t*b/(self.N-1),
                       W >= 2*dm.sum()*g,
                       w <= w_lim*cave,
                       cave*tau >= hin + 2*t,
                       Sy*(hin + t) <= I,
                      ]

        return constraints

    def loading(self, Wcent):
        return ChordSparL(self, Wcent)

