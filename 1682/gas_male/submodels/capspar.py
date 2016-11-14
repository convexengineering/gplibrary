" cap spar "
from gpkit import Model, Variable, vectorize
from constant_taper_chord import c_bar
from beam import Beam

class CapSpar(Model):
    "cap spar model"
    def __init__(self, b, cave, tau, N=5, **kwargs):
        self.N = N

        # phyiscal properties
        rhocfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        E = Variable("E", 2e7, "psi", "Youngs modulus of CFRP")

        with vectorize(self.N-1):
            t = Variable("t", "in", "spar cap thickness")
            hin = Variable("h_{in}", "in", "inner spar height")
            w = Variable("w", "in", "spar width")
            I = Variable("I", "m^4", "spar x moment of inertia")
            dm = Variable("dm", "kg", "segment spar mass")

        W = Variable("W", "lbf", "spar weight")
        w_lim = Variable("w_{lim}", 0.15, "-", "spar width to chord ratio")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        self.loading = CapSparL

        constraints = [I <= 2*w*t*(hin/2)**2,
                       dm >= rhocfrp*w*t*b/(self.N-1),
                       W >= 2*dm.sum()*g,
                       w <= w_lim*cave,
                       cave*tau >= hin + 2*t,
                       E == E,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class CapSparL(Model):
    "spar loading model"
    def __init__(self, static, Wcent, **kwargs):

        Nmax = Variable("N_{max}", 5, "-", "max loading")
        cbar = c_bar(0.5, static.N)
        sigmacfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")
        with vectorize(static.N-1):
            Mr = Variable("M_r", "N*m", "wing section root moment")

        beam = Beam(static.N, cbar)

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wcent/static["b"]**2),
            Mr == (beam["\\bar{M}"][:-1]*Wcent*Nmax*static["b"]/4),
            sigmacfrp >= Mr*(static["h_{in}"]+static["t"])/static["I"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            ]

        Model.__init__(self, None, [beam, constraints], **kwargs)
