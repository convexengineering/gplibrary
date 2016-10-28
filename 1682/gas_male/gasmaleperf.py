"""Jungle Hawk Owl Concept"""
import numpy as np
from gpkit import Model, Variable, vectorize
from gpkit.tools import te_exp_minus1

# pylint: disable=invalid-name

class Aircraft(Model):
    "the JHO vehicle"
    def __init__(self, **kwargs):
        self.dynamic_model = AircraftP
        self.fuse = Fuselage()
        self.wing = Wing()

        self.components = [self.fuse, self.wing]

        Wzfw = Variable("W_{zfw}", "lbf", "weight")
        constraints = [Wzfw >= sum(c["W"] for c in self.components)]

        super(Aircraft, self).__init__(None, self.components + constraints, **kwargs)


class AircraftP(Model):
    def __init__(self, static, state, **kwargs):

        Pshaft = Variable("P_{shaft}", "W", "shaft power")
        Wend = Variable("W_{end}", "lbf", "vector-end weight")
        Wstart = Variable("W_{start}", "lbf", "vector-begin weight")

        constraints = [Pshaft == Pshaft,
                       Wend == Wend,
                       Wstart == Wstart]

        self.wing = static.components[1].dynamic_model(static.wing, state)
        Model.__init__(self, None, [self.wing, constraints], **kwargs)

class FlightState(Model):
    "One chunk of a mission"
    def __init__(self, **kwargs):
        V = Variable("V", 40, "knots", "true airspeed")
        mu = Variable("\\mu", 1.628e-5, "N*s/m^2", "dynamic viscosity")
        rho = Variable("\\rho", 0.74, "kg/m^3", "air density")
        constraints = [V == V,
                       mu == mu,
                       rho == rho]
        super(FlightState, self).__init__(None, constraints, **kwargs)


class FlightSegment(Model):
    def __init__(self, aircraft, **kwargs):
        fs = FlightState()
        aircraftP = aircraft.dynamic_model(aircraft, fs)
        slf = SteadyLevelFlight(fs, aircraft, aircraftP)
        be = BreguetEndurance(aircraftP)

        Model.__init__(self, None, [fs, aircraft, aircraftP, slf, be], **kwargs)

class BreguetEndurance(Model):
    def __init__(self, aircraftP, **kwargs):
        z_bre = Variable("z_{bre}", "-", "Breguet coefficient")
        t = Variable("t", 1, "days", "Time per flight segment")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        bsfc = Variable("BSFC", "lb/hr/hp", 0.6,
                        "Brake specific fuel consumption")
        Wfuel = Variable("W_{fuel}", "lbf", "Segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")

        constraints = [
            z_bre >= (aircraftP["P_{shaft}"]*t*bsfc*g/
                      (aircraftP["W_{end}"]*aircraftP["W_{start}"])**0.5),
            # TCS([z_bre >= P_shafttot*t*bsfc*g/(Wend*Wstart)**0.5]),
            # TCS([z_bre >= P_shafttot*t*bsfc*g/Wend]),
            f_fueloil*Wfuel/aircraftP["W_{end}"] >= te_exp_minus1(z_bre, 3),
            aircraftP["W_{start}"] >= aircraftP["W_{end}"] + Wfuel
            ]

        Model.__init__(self, None, constraints, **kwargs)

class SteadyLevelFlight(Model):
    def __init__(self, state, aircraft, aircraftP, **kwargs):
        T = Variable("T", "N", "thrust")
        etaprop = Variable("\\eta_{prop}", 0.7, "-", "propulsive efficiency")

        constraints = [
            (aircraftP["W_{end}"]*aircraftP["W_{start}"])**0.5 <= (
                0.5*state["\\rho"]*state["V"]**2*aircraftP.wing["C_L"]
                * aircraft.wing["S"]),
            T == (0.5*state["\\rho"]*state["V"]**2*aircraftP.wing["C_D"]
                  *aircraft.wing["S"]),
            aircraftP["P_{shaft}"] == T*state["V"]/etaprop]

        Model.__init__(self, None, constraints, **kwargs)


class Wing(Model):
    "The thing that creates the lift"
    def __init__(self, **kwargs):
        W = Variable("W", "lbf", "weight")
        S = Variable("S", 190, "ft^2", "surface area")
        rho = Variable("\\rho", 1, "lbf/ft^2", "areal density")
        A = Variable("A", 27, "-", "aspect ratio")
        c = Variable("c", "ft", "mean chord")
        self.dynamic_model = WingP

        constraints = [W >= S*rho,
                       c == (S/A)**0.5]
        super(Wing, self).__init__(None, constraints, **kwargs)


class WingP(Model):
    def __init__(self, static, state, **kwargs):
        CD = Variable("C_D", "-", "drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        e = Variable("e", 0.9, "-", "Oswald efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        constraints = [
            CD >= (0.074/Re**0.2 + CL**2/np.pi/static["A"]/e),
            Re == state["\\rho"]*state["V"]*static["c"]/state["\\mu"],
            ]
        Model.__init__(self, None, constraints, **kwargs)


class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def __init__(self, **kwargs):
        V = Variable("V", 16, "gal", "volume")
        d = Variable("d", 12, "in", "diameter")
        # S = Variable("S", "ft^2", "wetted area")
        cd = Variable("c_d", .0047, "-", "drag coefficient")
        CDA = Variable("CDA", "ft^2", "drag area")
        W = Variable("W", 100, "lbf", "weight")

        constraints = [  # CDA >= cd*4*V/d,
            W == W,  # todo replace with model
            ]

        super(Fuselage, self).__init__(None, constraints, **kwargs)

class Mission(Model):
    def __init__(self, **kwargs):
        JHO = Aircraft()
        N = 4
        with vectorize(N):
            loiter = FlightSegment(JHO)

        loiter.substitutions["V"] = np.linspace(20, 40, N)
        mtow = Variable("MTOW", "lbf", "max-take off weight")
        Wfueltot = Variable("W_{fuel-tot}", "lbf", "total fuel weight")

        constraints = [mtow >= loiter["W_{start}"][0],
                       mtow >= JHO["W_{zfw}"] + Wfueltot,
                       Wfueltot >= loiter["W_{fuel}"].sum(),
                       loiter["W_{end}"][-1] >= JHO["W_{zfw}"],
                       loiter["W_{end}"][:-1] >= loiter["W_{start}"][1:]
                      ]

        Model.__init__(self, mtow, [JHO, loiter, constraints], **kwargs)


if __name__ == "__main__":
    M = Mission()
    # JHO.debug(solver="mosek")
    sol = M.solve("mosek")
    print sol.table()
