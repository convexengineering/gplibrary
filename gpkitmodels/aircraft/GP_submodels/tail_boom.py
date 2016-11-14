" tail boom model "
import numpy as np
from gpkit import Variable, Model

class TailBoom(Model):
    "tail boom model"
    def __init__(self, **kwargs):

        l = Variable("l", "ft", "tail boom length")
        E = Variable("E", 150e9, "N/m^2", "young's modulus carbon fiber")
        k = Variable("k", 0.8, "-", "tail boom inertia value")
        kfac = Variable("(1-k/2)", 1-k.value/2, "-", "(1-k/2)")
        I0 = Variable("I_0", "m^4", "tail boom moment of inertia")
        d0 = Variable("d_0", "ft", "tail boom diameter")
        t0 = Variable("t_0", "mm", "tail boom thickness")
        tmin = Variable("t_{min}", 0.25, "mm", "minimum tail boom thickness")
        rhocfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        W = Variable("W", "lbf", "tail boom weight")
        J = Variable("J", "m^4", "tail boom polar moment of inertia")
        S = Variable("S", "ft**2", "tail boom surface area")

        self.flight_model = TailBoomAero
        self.horizontalbending = HorizontalBoomBending
        self.verticalbending = VerticalBoomBending
        self.verticaltorsion = VerticalBoomTorsion

        constraints = [
            I0 <= np.pi*t0*d0**3/8.0,
            W >= np.pi*g*rhocfrp*d0*l*t0*kfac,
            t0 >= tmin,
            J <= np.pi/8.0*d0**3*t0,
            S == l*np.pi*d0,
            k == k,
            E == E
            ]

        Model.__init__(self, None, constraints, **kwargs)

class TailBoomAero(Model):
    "horizontal tail aero model"
    def __init__(self, static, state, **kwargs):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["l"]/state["\\mu"]),
            Cf >= 0.455/Re**0.3,
            ]

        Model.__init__(self, None, constraints, **kwargs)

class TailBoomState(Model):
    "tail boom design state"
    def __init__(self, **kwargs):

        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 40, "m/s", "never exceed vehicle speed")

        constraints = [rhosl == rhosl,
                       Vne == Vne]

        Model.__init__(self, None, constraints, **kwargs)

class VerticalBoomTorsion(Model):
    "tail boom torison case"
    def __init__(self, tailboom, vtail, state, **kwargs):

        T = Variable("T", "N*m", "vertical tail moment")
        taucfrp = Variable("\\tau_{CFRP}", 210, "MPa", "torsional stress limit")

        constraints = [
            T >= (0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*vtail["S"]
                  * vtail["C_{L_{max}}"]*vtail["b_v"]),
            taucfrp >= T*tailboom["d_0"]/2/tailboom["J"]
            ]

        Model.__init__(self, None, constraints, **kwargs)

class VerticalBoomBending(Model):
    "tail boom bending loading case"
    def __init__(self, tailboom, vtail, state, **kwargs):

        F = Variable("F", "N", "vertical tail force")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.3, "-",
                         "max tail boom deflection angle")

        constraints = [
            F >= (0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*vtail["S"]
                  * vtail["C_{L_{max}}"]),
            th >= (F*tailboom["l"]**2/tailboom["E"]/tailboom["I_0"]
                   * (1+tailboom["k"])/2),
            th <= thmax,
            ]

        Model.__init__(self, None, constraints, **kwargs)

class HorizontalBoomBending(Model):
    "tail boom bending loading case"
    def __init__(self, tailboom, htail, state, **kwargs):

        F = Variable("F", "N", "horizontal tail force")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.3, "-",
                         "max tail boom deflection angle")

        constraints = [
            F >= (0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*htail["S"]
                  * htail["C_{L_{max}}"]),
            th >= (F*tailboom["l"]**2/tailboom["E"]/tailboom["I_0"]
                   * (1+tailboom["k"])/2),
            th <= thmax,
            ]

        Model.__init__(self, None, constraints, **kwargs)
