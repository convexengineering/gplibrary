" tail boom model "
import numpy as np
from gpkit import Variable, Model

class TailBoom(Model):
    "tail boom model"
    def setup(self):

        l = Variable("l", "ft", "tail boom length")
        E = Variable("E", 150e9, "N/m^2", "young's modulus carbon fiber")
        k = Variable("k", 0.8, "-", "tail boom inertia value")
        kfac = Variable("(1-k/2)", 1-k.value/2, "-", "(1-k/2)")
        I0 = Variable("I_0", "m^4", "tail boom moment of inertia")
        d0 = Variable("d_0", "in", "tail boom diameter")
        t0 = Variable("t_0", "mm", "tail boom thickness")
        tmin = Variable("t_{min}", 0.25, "mm", "minimum tail boom thickness")
        rhocfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        W = Variable("W", "lbf", "tail boom weight")
        J = Variable("J", "m^4", "tail boom polar moment of inertia")
        S = Variable("S", "ft**2", "tail boom surface area")
        mfac = Variable("m_{fac}", 1.0, "-", "tail boom margin factor")

        constraints = [
            I0 <= np.pi*t0*d0**3/8.0,
            W/mfac >= np.pi*g*rhocfrp*d0*l*t0*kfac,
            t0 >= tmin,
            J <= np.pi/8.0*d0**3*t0,
            S == l*np.pi*d0,
            k == k,
            E == E
            ]

        return constraints

    def flight_model(self, state):
        return TailBoomAero(self, state)

    def horizontalbending(self, htail, state):
        return HorizontalBoomBending(self, htail, state)

    def verticalbending(self, vtail, state):
        return VerticalBoomBending(self, vtail, state)

    def verticaltorsion(self, vtail, state):
        return VerticalBoomTorsion(self, vtail, state)

class TailBoomAero(Model):
    "horizontal tail aero model"
    def setup(self, static, state):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["l"]/state["\\mu"]),
            Cf >= 0.455/Re**0.3,
            ]

        return constraints

class TailBoomState(Model):
    "tail boom design state"
    def setup(self):

        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vne = Variable("V_{NE}", 40, "m/s", "never exceed vehicle speed")

        constraints = [rhosl == rhosl,
                       Vne == Vne]

        return constraints

class VerticalBoomTorsion(Model):
    "tail boom torison case"
    def setup(self, tailboom, vtail, state):

        T = Variable("T", "N*m", "vertical tail moment")
        taucfrp = Variable("\\tau_{CFRP}", 210, "MPa", "torsional stress limit")

        constraints = [
            T >= (0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*vtail.planform.S
                  * vtail.planform.CLmax*vtail.planform.b),
            taucfrp >= T*tailboom["d_0"]/2/tailboom["J"]
            ]

        return constraints

class VerticalBoomBending(Model):
    "tail boom bending loading case"
    def setup(self, tailboom, vtail, state):

        F = Variable("F", "N", "vertical tail force")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.1, "-",
                         "max tail boom deflection angle")

        constraints = [
            F >= (0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*vtail.planform.S
                  * vtail.planform.CLmax),
            th >= (F*tailboom["l"]**2/tailboom["E"]/tailboom["I_0"]
                   * (1+tailboom["k"])/2),
            th <= thmax,
            ]

        return constraints

class HorizontalBoomBending(Model):
    "tail boom bending loading case"
    def setup(self, tailboom, htail, state):

        F = Variable("F", "N", "horizontal tail force")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.1, "-",
                         "max tail boom deflection angle")

        constraints = [
            F >= (0.5*state["\\rho_{sl}"]*state["V_{NE}"]**2*htail.planform.S
                  * htail.planform.CLmax),
            th >= (F*tailboom["l"]**2/tailboom["E"]/tailboom["I_0"]
                   * (1+tailboom["k"])/2),
            th <= thmax,
            ]

        return constraints
