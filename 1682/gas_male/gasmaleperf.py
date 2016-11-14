"""Jungle Hawk Owl"""
import numpy as np
from submodels.breguet_endurance import BreguetEndurance
from submodels.flight_state import FlightState
from submodels.gas_engine import Engine
from submodels.wing import Wing
from submodels.fuselage import Fuselage
from submodels.empennage import Empennage, TailBoomState, TailBoomFlexibility
from helpers import summing_vars
from gpkit import Model, Variable, vectorize, units

# pylint: disable=invalid-name

class Aircraft(Model):
    "the JHO vehicle"
    def __init__(self, Wfueltot, DF70=False, **kwargs):
        self.flight_model = AircraftPerf
        self.fuselage = Fuselage(Wfueltot)
        self.wing = Wing()
        self.engine = Engine(DF70)
        self.empennage = Empennage()

        components = [self.fuselage, self.wing, self.engine, self.empennage]
        self.smeared_loads = [self.fuselage, self.engine]

        self.loading = AircraftLoading

        Wzfw = Variable("W_{zfw}", "lbf", "zero fuel weight")
        Wpay = Variable("W_{pay}", 10, "lbf", "payload weight")
        Wavn = Variable("W_{avn}", 8, "lbf", "avionics weight")

        constraints = [
            Wzfw >= sum(summing_vars(components, "W")) + Wpay + Wavn,
            self.empennage.horizontaltail["V_h"] <= (
                self.empennage.horizontaltail["S"]
                * self.empennage.horizontaltail["l_h"]/self.wing["S"]**2
                * self.wing["b"]),
            self.empennage.verticaltail["V_v"] <= (
                self.empennage.verticaltail["S"]
                * self.empennage.verticaltail["l_v"]/self.wing["S"]
                / self.wing["b"]),
            self.wing["C_{L_{max}}"]/self.wing["m_w"] <= (
                self.empennage.horizontaltail["C_{L_{max}}"]
                / self.empennage.horizontaltail["m_h"])
            ]

        Model.__init__(self, None, [components, constraints],
                       **kwargs)

class AircraftLoading(Model):
    "aircraft loading model"
    def __init__(self, aircraft, Wcent, **kwargs):

        loading = [aircraft.wing.loading(aircraft.wing, Wcent)]
        loading.append(aircraft.empennage.loading(aircraft.empennage))
        loading.append(aircraft.fuselage.loading(aircraft.fuselage, Wcent))

        tbstate = TailBoomState()
        loading.append(TailBoomFlexibility(aircraft.empennage.horizontaltail,
                                           aircraft.empennage.tailboom,
                                           aircraft.wing, tbstate, **kwargs))

        Model.__init__(self, None, loading, **kwargs)

class AircraftPerf(Model):
    "performance model for aircraft"
    def __init__(self, static, state, **kwargs):

        self.wing = static.wing.flight_model(static.wing, state)
        self.fuselage = static.fuselage.flight_model(static.fuselage, state)
        self.engine = static.engine.flight_model(static.engine, state)
        self.htail = static.empennage.horizontaltail.flight_model(
            static.empennage.horizontaltail, state)
        self.vtail = static.empennage.verticaltail.flight_model(
            static.empennage.verticaltail, state)
        self.tailboom = static.empennage.tailboom.flight_model(
            static.empennage.tailboom, state)

        self.dynamicmodels = [self.wing, self.fuselage, self.engine,
                              self.htail, self.vtail, self.tailboom]
        areadragmodel = [self.fuselage, self.htail, self.vtail, self.tailboom]
        areadragcomps = [static.fuselage, static.empennage.horizontaltail,
                         static.empennage.verticaltail,
                         static.empennage.tailboom]

        Wend = Variable("W_{end}", "lbf", "vector-end weight")
        Wstart = Variable("W_{start}", "lbf", "vector-begin weight")
        CD = Variable("C_D", "-", "drag coefficient")
        CDA = Variable("CDA", "-", "area drag coefficient")
        mfac = Variable("m_{fac}", 1.7, "-", "drag margin factor")

        dvars = []
        for dc, dm in zip(areadragcomps, areadragmodel):
            if "C_f" in dm.varkeys:
                dvars.append(dm["C_f"]*dc["S"]/static.wing["S"])

        constraints = [Wend == Wend,
                       Wstart == Wstart,
                       CDA/mfac >= sum(dvars),
                       CD >= CDA + self.wing["C_d"]]

        Model.__init__(self, None, [self.dynamicmodels, constraints], **kwargs)

class FlightSegment(Model):
    "creates flight segment for aircraft"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 etap=0.7, **kwargs):

        self.aircraft = aircraft

        with vectorize(N):
            self.fs = FlightState(alt, onStation, wind)
            self.aircraftPerf = self.aircraft.flight_model(self.aircraft,
                                                           self.fs)
            self.slf = SteadyLevelFlight(self.fs, self.aircraft,
                                         self.aircraftPerf, etap)
            self.be = BreguetEndurance(self.aircraftPerf)

        self.submodels = [self.fs, self.aircraftPerf, self.slf, self.be]

        Wfuelfs = Variable("W_{fuel-fs}", "lbf", "flight segment fuel weight")

        self.constraints = [Wfuelfs >= self.be["W_{fuel}"].sum()]

        if N > 1:
            self.constraints.extend([self.aircraftPerf["W_{end}"][:-1] >=
                                     self.aircraftPerf["W_{start}"][1:]])

        Model.__init__(self, None, [self.aircraft, self.submodels,
                                    self.constraints], **kwargs)

class Loiter(Model):
    "make a loiter flight segment"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 etap=0.7, **kwargs):
        fs = FlightSegment(N, aircraft, alt, onStation, wind, etap)

        t = Variable("t", 6, "days", "time loitering")
        constraints = [fs.be["t"] >= t/N]

        Model.__init__(self, None, [constraints, fs], **kwargs)

class Cruise(Model):
    "make a cruise flight segment"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 etap=0.7, R=200, **kwargs):
        fs = FlightSegment(N, aircraft, alt, onStation, wind, etap)

        R = Variable("R", R, "nautical_miles", "Range to station")
        constraints = [R/N <= fs["V"]*fs.be["t"]]

        Model.__init__(self, None, [fs, constraints], **kwargs)

class Climb(Model):
    "make a climb flight segment"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 etap=0.7, dh=15000, **kwargs):
        fs = FlightSegment(N, aircraft, alt, onStation, wind, etap)

        with vectorize(N):
            hdot = Variable("\\dot{h}", "ft/min", "Climb rate")

        deltah = Variable("\\Delta_h", dh, "ft", "altitude difference")
        hdotmin = Variable("\\dot{h}_{min}", 100, "ft/min",
                           "minimum climb rate")

        constraints = [
            hdot*fs.be["t"] >= deltah/N,
            hdot >= hdotmin,
            fs.slf["T"] >= (0.5*fs["\\rho"]*fs["V"]**2*fs["C_D"]
                            * fs.aircraft.wing["S"] + fs["W_{start}"]*hdot
                            / fs["V"]),
            ]

        Model.__init__(self, None, [fs, constraints], **kwargs)

class SteadyLevelFlight(Model):
    "steady level flight model"
    def __init__(self, state, aircraft, perf, etap, **kwargs):

        T = Variable("T", "N", "thrust")
        etaprop = Variable("\\eta_{prop}", etap, "-", "propulsive efficiency")

        constraints = [
            (perf["W_{end}"]*perf["W_{start}"])**0.5 <= (
                0.5*state["\\rho"]*state["V"]**2*perf["C_L"]
                * aircraft.wing["S"]),
            T >= (0.5*state["\\rho"]*state["V"]**2*perf["C_D"]
                  *aircraft.wing["S"]),
            perf["P_{shaft}"] >= T*state["V"]/etaprop]

        Model.__init__(self, None, constraints, **kwargs)

class Mission(Model):
    "creates flight profile"
    def __init__(self, DF70=False, **kwargs):

        mtow = Variable("MTOW", "lbf", "max-take off weight")
        Wcent = Variable("W_{cent}", "lbf", "center aircraft weight")
        Wfueltot = Variable("W_{fuel-tot}", "lbf", "total aircraft fuel weight")

        JHO = Aircraft(Wfueltot, DF70)
        loading = JHO.loading(JHO, Wcent)

        climb1 = Climb(10, JHO, alt=np.linspace(0, 15000, 11)[1:], etap=0.508)
        cruise1 = Cruise(1, JHO, etap=0.684, R=180)
        loiter1 = Loiter(5, JHO, etap=0.647, onStation=True)
        cruise2 = Cruise(1, JHO, etap=0.684)
        mission = [climb1, cruise1, loiter1, cruise2]

        constraints = [
            mtow >= JHO["W_{zfw}"] + Wfueltot,
            Wfueltot >= sum(fs["W_{fuel-fs}"] for fs in mission),
            mission[-1]["W_{end}"][-1] >= JHO["W_{zfw}"],
            Wcent >= Wfueltot + sum(summing_vars(JHO.smeared_loads, "W"))
            ]

        for i, fs in enumerate(mission[1:]):
            constraints.extend([
                mission[i]["W_{end}"][-1] == fs["W_{start}"][0]
                ])

        Model.__init__(self, mtow, [JHO, mission, loading, constraints],
                       **kwargs)


if __name__ == "__main__":
    M = Mission(DF70=True)
    # JHO.debug(solver="mosek")
    sol = M.solve("mosek")
    print sol.table()
