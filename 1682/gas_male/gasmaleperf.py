"""Jungle Hawk Owl"""
import numpy as np
from submodels.breguet_endurance import BreguetEndurance
from submodels.flight_state import FlightState
from submodels.gas_engine import Engine
from submodels.wing import Wing
from submodels.fuselage import Fuselage
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

class Empennage(Model):
    "empennage model, consisting of vertical, horizontal and tailboom"
    def __init__(self, **kwargs):
        mfac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "empennage weight")

        self.horizontaltail = HorizontalTail()
        self.verticaltail = VerticalTail()
        self.tailboom = TailBoom()
        self.components = [self.horizontaltail, self.verticaltail,
                           self.tailboom]

        self.loading = EmpennageLoading

        constraints = [
            W/mfac >= (self.horizontaltail["W"] + self.verticaltail["W"]
                       + self.tailboom["W"]),
            self.tailboom["l"] >= self.horizontaltail["l_h"],
            self.tailboom["l"] >= self.verticaltail["l_v"],
            ]

        Model.__init__(self, None, [self.components, constraints],
                       **kwargs)

class HorizontalTail(Model):
    "horizontal tail model"
    def __init__(self, **kwargs):
        Sh = Variable("S", "ft**2", "horizontal tail area")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        ARh = Variable("AR_h", "-", "horizontal tail aspect ratio")
        Abar = Variable("\\bar{A}_{NACA0008}", 0.0548, "-",
                        "cross sectional area of NACA 0008")
        rhofoam = Variable("\\rho_{foam}", 1.5, "lbf/ft^3",
                           "Density of formular 250")
        rhoskin = Variable("\\rho_{skin}", 0.1, "g/cm**2",
                           "horizontal tail skin density")
        bh = Variable("b_h", "ft", "horizontal tail span")
        W = Variable("W", "lbf", "horizontal tail weight")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        lh = Variable("l_h", "ft", "horizontal tail moment arm")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")
        cth = Variable("c_{t_h}", "ft", "horizontal tail tip chord")
        lamhfac = Variable("\\lambda_h/(\\lambda_h+1)", 1.0/(1.0+1), "-",
                           "horizontal tail taper ratio factor")
        CLhtmax = Variable("C_{L_{max}}", "-", "maximum CL of horizontal tail")

        self.flight_model = HorizontalTailAero

        constraints = [
            bh**2 == ARh*Sh,
            mh*(1+2/ARh) <= 2*np.pi,
            W >= g*rhoskin*Sh + rhofoam*Sh**2/bh*Abar,
            cth == 2*Sh/bh*lamhfac,
            lh == lh,
            CLhmin == CLhmin,
            CLhtmax == CLhtmax,
            Vh == Vh,
            ]

        Model.__init__(self, None, constraints, **kwargs)

class HorizontalTailAero(Model):
    "horizontal tail aero model"
    def __init__(self, static, state, **kwargs):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["S"]/static["b_h"]
                   / state["\\mu"]),
            Cf >= 0.455/Re**0.3,
            ]

        Model.__init__(self, None, constraints, **kwargs)

class VerticalTail(Model):
    "vertical tail model"
    def __init__(self, **kwargs):

        W = Variable("W", "lbf", "one vertical tail weight")
        Sv = Variable("S", "ft**2", "total vertical tail surface area")
        Vv = Variable("V_v", 0.025, "-", "vertical tail volume coefficient")
        ARv = Variable("AR_v", "-", "vertical tail aspect ratio")
        bv = Variable("b_v", "ft", "one vertical tail span")
        rhofoam = Variable("\\rho_{foam}", 1.5, "lbf/ft^3",
                           "Density of formular 250")
        rhoskin = Variable("\\rho_{skin}", 0.1, "g/cm**2",
                           "vertical tail skin density")
        Abar = Variable("\\bar{A}_{NACA0008}", 0.0548, "-",
                        "cross sectional area of NACA 0008")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        lv = Variable("l_v", "ft", "horizontal tail moment arm")
        ctv = Variable("c_{t_v}", "ft", "vertical tail tip chord")
        lamvfac = Variable("\\lambda_v/(\\lambda_v+1)", 1.0/(1.0+1), "-",
                           "vertical tail taper ratio factor")
        CLvtmax = Variable("C_{L_{max}}", 1.1, "-",
                           "maximum CL of vertical tail")
        lantenna = Variable("l_{antenna}", 13.4, "in", "antenna length")
        wantenna = Variable("w_{antenna}", 10.2, "in", "antenna width")

        self.flight_model = VerticalTailAero

        constraints = [Vv == Vv,
                       lv == lv,
                       bv**2 == ARv*Sv,
                       W >= rhofoam*Sv**2/bv*Abar + g*rhoskin*Sv,
                       ctv == 2*Sv/bv*lamvfac,
                       ctv >= wantenna*1.3,
                       bv >= lantenna,
                       CLvtmax == CLvtmax,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class VerticalTailAero(Model):
    "horizontal tail aero model"
    def __init__(self, static, state, **kwargs):

        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == (state["V"]*state["\\rho"]*static["S"]/static["b_v"]
                   / state["\\mu"]),
            Cf >= 0.455/Re**0.3,
            ]

        Model.__init__(self, None, constraints, **kwargs)


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

        self.case = TailBoomState()
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

class TailBoomFlexibility(Model):
    "tail boom flexibility model"
    def __init__(self, htail, tailboom, wing, state, **kwargs):

        Fne = Variable("F_{NE}", "-", "tail boom flexibility factor")
        deda = Variable("d\\epsilon/d\\alpha", "-", "wing downwash derivative")
        SMcorr = Variable("SM_{corr}", 0.35, "-", "corrected static margin")

        # signomial helper variables
        sph1 = Variable("sph1", "-", "first term involving $V_h$")
        sph2 = Variable("sph2", "-", "second term involving $V_h$")

        constraints = [
            Fne >= (1 + htail["m_h"]*0.5*state["V_{NE}"]**2*state["\\rho_{sl}"]
                    * htail["S"]*tailboom["l"]**2/tailboom["E"]
                    / tailboom["I_0"]*tailboom["(1-k/2)"]),
            sph1*(wing["m_w"]*Fne/htail["m_h"]/htail["V_h"]) + deda <= 1,
            sph2 <= htail["V_h"]*htail["(C_{L_h})_{min}"]/wing["C_{L_{max}}"],
            (sph1 + sph2).mono_lower_bound({"sph1": .48, "sph2": .52}) >= (
                SMcorr + wing["C_M"]/wing["C_{L_{max}}"]),
            deda >= wing["m_w"]*wing["S"]/wing["b"]/4/np.pi/htail["l_h"]]

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

class EmpennageLoading(Model):
    "tail boom loading case"
    def __init__(self, empennage, **kwargs):
        state = TailBoomState()

        loading = [empennage.tailboom.horizontalbending(
            empennage.tailboom, empennage.horizontaltail, state)]
        loading.append(empennage.tailboom.verticalbending(
            empennage.tailboom, empennage.verticaltail, state))
        loading.append(empennage.tailboom.verticaltorsion(
            empennage.tailboom, empennage.verticaltail, state))

        Model.__init__(self, None, loading, **kwargs)

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
