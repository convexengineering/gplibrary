" wing.py "
import numpy as np
from gpkit import Variable, Model, vectorize

class Wing(Model):
    "The thing that creates the lift"
    def __init__(self, N=5, lam=0.5, **kwargs):

        W = Variable("W", "lbf", "weight")
        mfac = Variable("m_{fac}", 1.2, "-", "wing weight margin factor")
        S = Variable("S", "ft^2", "surface area")
        A = Variable("A", "-", "aspect ratio")
        b = Variable("b", "ft", "wing span")
        tau = Variable("\\tau", 0.115, "-", "airfoil thickness ratio")
        CLmax = Variable("C_{L_{max}}", 1.39, "-", "maximum CL of JHO1")
        CM = Variable("C_M", 0.14, "-", "wing moment coefficient")
        mw = Variable("m_w", 2.0*np.pi/(1+2.0/23), "-",
                      "assumed span wise effectiveness")
        croot = Variable("c_{root}", "ft", "root chord")
        cmac = Variable("c_{MAC}", "ft", "mean aerodynamic chord")
        cb = c_bar(lam, N)
        with vectorize(N):
            cbar = Variable("\\bar{c}", cb, "-",
                            "normalized chord at mid element")
        with vectorize(N-1):
            cave = Variable("c_{ave}", "ft", "mid section chord")

        self.flight_model = WingAero

        constraints = [b**2 == S*A,
                       tau == tau,
                       CLmax == CLmax,
                       CM == CM,
                       mw == mw,
                       cbar == cbar,
                       cave == (cb[1:] + cb[:-1])/2*S/b,
                       croot == S/b*cb[0],
                       cmac == S/b]

        self.capspar = CapSpar(b, cave, tau, N)
        self.wingskin = WingSkin(S, croot, b)
        self.winginterior = WingInterior(cave, b, N)
        self.components = [self.capspar, self.wingskin, self.winginterior]
        self.loading = WingLoading

        constraints.extend([W/mfac >= sum(c["W"] for c in self.components)])

        Model.__init__(self, None, [self.components, constraints],
                       **kwargs)

def c_bar(lam, N):
    "returns wing chord lengths for constant taper wing"
    eta = np.linspace(0, 1, N)
    c = 2/(1+lam)*(1+(lam-1)*eta)
    return c

class WingLoading(Model):
    "wing loading cases"
    def __init__(self, wing, Wcent, **kwargs):

        skinloading = wing.wingskin.loading(wing.wingskin)
        caploading = wing.capspar.loading(wing.capspar, Wcent)

        Model.__init__(self, None, [skinloading, caploading], **kwargs)

class WingInterior(Model):
    "wing interior model"
    def __init__(self, cave, b, N, **kwargs):

        W = Variable("W", "lbf", "interior mass of wing")
        rhofoam = Variable("\\rho_{foam}", 0.036, "g/cm^3", "foam density")
        Abar = Variable("\\bar{A}_{jh01}", 0.0753449, "-",
                        "jh01 non dimensional area")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [W >= 2*(g*rhofoam*Abar*cave**2*(b/2)/(N-1)).sum()]

        Model.__init__(self, None, constraints, **kwargs)

class WingSkin(Model):
    "wing skin model"
    def __init__(self, S, croot, b, **kwargs):

        rhocfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        W = Variable("W", "lbf", "wing skin weight")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")
        t = Variable("t", "in", "wing skin thickness")
        tmin = Variable("t_{min}", 0.012, "in",
                        "minimum gague wing skin thickness")
        Jtbar = Variable("\\bar{J/t}", 0.01114, "1/mm",
                         "torsional moment of inertia")

        self.loading = WingSkinL

        constraints = [W >= rhocfrp*S*2*t*g,
                       t >= tmin,
                       Jtbar == Jtbar,
                       b == b,
                       croot == croot]

        Model.__init__(self, None, constraints, **kwargs)

class WingSkinL(Model):
    "wing skin loading model for torsional loads in skin"
    def __init__(self, static, **kwargs):

        taucfrp = Variable("\\tau_{CFRP}", 570, "MPa", "torsional stress limit")
        Cmw = Variable("C_{m_w}", 0.121, "-", "negative wing moment coefficent")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 45, "m/s", "never exceed vehicle speed")

        constraints = [
            taucfrp >= (1/static["\\bar{J/t}"]/(static["c_{root}"])**2
                        / static["t"]*Cmw*static["S"]*rhosl*Vne**2)]

        Model.__init__(self, None, constraints, **kwargs)

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

class Beam(Model):
    "discretized beam bending model"
    def __init__(self, N, q, **kwargs):

        with vectorize(N-1):
            EIbar = Variable("\\bar{EI}", "-",
                             "normalized YM and moment of inertia")

        with vectorize(N):
            qbar = Variable("\\bar{q}", q, "-", "normalized loading")
            Sbar = Variable("\\bar{S}", "-", "normalized shear")
            Mbar = Variable("\\bar{M}", "-", "normalized moment")
            th = Variable("\\theta", "-", "deflection slope")
            dbar = Variable("\\bar{\\delta}", "-", "normalized displacement")


        Sbartip = Variable("\\bar{S}_{tip}", 1e-10, "-", "Tip loading")
        Mbartip = Variable("\\bar{M}_{tip}", 1e-10, "-", "Tip moment")
        throot = Variable("\\theta_{root}", 1e-10, "-", "Base angle")
        dbarroot = Variable("\\bar{\\delta}_{root}", 1e-10, "-",
                            "Base deflection")
        dx = Variable("dx", "-", "normalized length of element")

        constraints = [
            Sbar[:-1] >= Sbar[1:] + 0.5*dx*(qbar[:-1] + qbar[1:]),
            Sbar[-1] >= Sbartip,
            Mbar[:-1] >= Mbar[1:] + 0.5*dx*(Sbar[:-1] + Sbar[1:]),
            Mbar[-1] >= Mbartip,
            th[0] >= throot,
            th[1:] >= th[:-1] + 0.5*dx*(Mbar[1:] + Mbar[:-1])/EIbar,
            dbar[0] >= dbarroot,
            dbar[1:] >= dbar[:-1] + 0.5*dx*(th[1:] + th[:-1]),
            1 == (N-1)*dx,
            ]

        Model.__init__(self, None, constraints, **kwargs)


class WingAero(Model):
    "wing aerodynamic model with profile and induced drag"
    def __init__(self, static, state, **kwargs):
        "wing drag model"
        Cd = Variable("C_d", "-", "wing drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        e = Variable("e", 0.9, "-", "Oswald efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        cdp = Variable("c_{dp}", "-", "wing profile drag coeff")

        constraints = [
            Cd >= cdp + CL**2/np.pi/static["A"]/e,
            cdp**3.72 >= (0.0247*CL**2.49*Re**-1.11
                          + 2.03e-7*CL**12.7*Re**-0.338
                          + 6.35e10*CL**-0.243*Re**-3.43
                          + 6.49e-6*CL**-1.9*Re**-0.681),
            Re == state["\\rho"]*state["V"]*static["c_{MAC}"]/state["\\mu"],
            ]

        Model.__init__(self, None, constraints, **kwargs)

