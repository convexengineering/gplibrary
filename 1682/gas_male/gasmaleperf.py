"""Jungle Hawk Owl"""
import numpy as np
from gpkit import Model, Variable, vectorize, units
from gpkit.tools import te_exp_minus1

# pylint: disable=invalid-name

class Aircraft(Model):
    "the JHO vehicle"
    def __init__(self, DF70=False, **kwargs):
        self.dynamic_model = AircraftP
        self.fuse = Fuselage()
        self.wing = Wing()
        self.engine = Engine(DF70)

        self.components = [self.fuse, self.wing, self.engine]

        Wzfw = Variable("W_{zfw}", "lbf", "weight")
        constraints = [Wzfw >= sum(c["W"] for c in self.components)]

        Model.__init__(self, None, self.components + constraints, **kwargs)

class AircraftP(Model):
    "performance model for aircraft"
    def __init__(self, static, state, **kwargs):

        self.dmodels = []
        self.dragcomps = []
        for c in static.components:
            if "dynamic_model" in c.__dict__.keys():
                dm = c.dynamic_model(c, state)
                self.dmodels.append(dm)
                if "C_d" in dm.varkeys:
                    self.dragcomps.append(dm)

        Wend = Variable("W_{end}", "lbf", "vector-end weight")
        Wstart = Variable("W_{start}", "lbf", "vector-begin weight")
        CD = Variable("C_D", "-", "drag coefficient")

        constraints = [Wend == Wend,
                       Wstart == Wstart,
                       CD >= sum(dm["C_d"] for dm in self.dragcomps)]

        Model.__init__(self, None, [self.dmodels, constraints], **kwargs)

class FlightState(Model):
    "One chunk of a mission"
    def __init__(self, alt, onStation, **kwargs):

        self.onStation = onStation

        V = Variable("V", 25, "m/s", "true airspeed")
        mu = Variable("\\mu", 1.628e-5, "N*s/m^2", "dynamic viscosity")
        rho = Variable("\\rho", "kg/m^3", "air density")
        h = Variable("h", alt, "ft", "altitude")
        href = Variable("h_{ref}", 15000, "ft", "Reference altitude")
        psl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        Latm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        Tsl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        temp = [(t.value - l.value*v.value).magnitude
                for t, v, l in zip(Tsl, h, Latm)]
        Tatm = Variable("t_{atm}", temp, "K", "Air temperature")
        mu = Variable("\\mu", "N*s/m^2", "Dynamic viscosity")
        musl = Variable("\\mu_{sl}", 1.789*10**-5, "N*s/m^2",
                        "Dynamic viscosity at sea level")
        Rspec = Variable("R_{spec}", 287.058, "J/kg/K",
                         "Specific gas constant of air")

        # Atmospheric variation with altitude (valid from 0-7km of altitude)
        constraints = [rho == psl*Tatm**(5.257-1)/Rspec/(Tsl**5.257),
                       (mu/musl)**0.1 == 0.991*(h/href)**(-0.00529),
                       V == V,
                       mu == mu,
                       h == h,
                       href == href]
        Model.__init__(self, None, constraints, **kwargs)


class FlightSegment(Model):
    "creates flight segment for aircraft"
    def __init__(self, N, aircraft, alt=15000, onStation=False, **kwargs):
        with vectorize(N):
            fs = FlightState(alt, onStation)
            aircraftP = aircraft.dynamic_model(aircraft, fs)
            slf = SteadyLevelFlight(fs, aircraft, aircraftP)
            be = BreguetEndurance(aircraftP)

        Wfuelfs = Variable("W_{fuel-fs}", "lbf", "flight segment fuel weight")

        constraints = [Wfuelfs >= be["W_{fuel}"].sum(),
                       aircraftP["W_{end}"][:-1] >= aircraftP["W_{start}"][1:]
                      ]

        Model.__init__(self, None, [fs, aircraft, aircraftP, slf, be,
                                    constraints], **kwargs)

class BreguetEndurance(Model):
    "breguet endurance model"
    def __init__(self, perf, **kwargs):
        z_bre = Variable("z_{bre}", "-", "Breguet coefficient")
        t = Variable("t", 1, "days", "Time per flight segment")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        # # bsfc = Variable("BSFC", "lb/hr/hp", 0.6,
        # #                 "Brake specific fuel consumption")
        Wfuel = Variable("W_{fuel}", "lbf", "Segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")

        constraints = [
            z_bre >= (perf["P_{total}"]*t*perf["BSFC"]*g/
                      (perf["W_{end}"]*perf["W_{start}"])**0.5),
            # TCS([z_bre >= Ptotal*t*bsfc*g/(Wend*Wstart)**0.5]),
            # TCS([z_bre >= Ptotal*t*bsfc*g/Wend]),
            f_fueloil*Wfuel/perf["W_{end}"] >= te_exp_minus1(z_bre, 3),
            perf["W_{start}"] >= perf["W_{end}"] + Wfuel
            ]

        Model.__init__(self, None, constraints, **kwargs)

class SteadyLevelFlight(Model):
    "steady level flight model"
    def __init__(self, state, aircraft, perf, **kwargs):

        T = Variable("T", "N", "thrust")
        etaprop = Variable("\\eta_{prop}", 0.7, "-", "propulsive efficiency")

        constraints = [
            (perf["W_{end}"]*perf["W_{start}"])**0.5 <= (
                0.5*state["\\rho"]*state["V"]**2*perf["C_L"]
                * aircraft.wing["S"]),
            T == (0.5*state["\\rho"]*state["V"]**2*perf["C_D"]
                  *aircraft.wing["S"]),
            perf["P_{shaft}"] == T*state["V"]/etaprop]

        Model.__init__(self, None, constraints, **kwargs)

class Engine(Model):
    "engine model"
    def __init__(self, DF70=False, **kwargs):

        self.DF70 = DF70

        W = Variable("W", "lbf", "Installed/Total engine weight")
        mfac = Variable("m_{fac}", 1.0, "-", "Engine weight margin factor")

        if DF70:
            Wdf70 = Variable("W_{DF70}", 7.1, "lbf",
                             "Installed/Total DF70 engine weight")
            Pslmax = Variable("P_{sl-max}", 5.17, "hp",
                              "Max shaft power at sea level")
            constraints = [W/mfac >= Wdf70]

        else:
            Pref = Variable("P_{ref}", 2.295, "hp", "Reference shaft power")
            Wengref = Variable("W_{eng-ref}", 4.4107, "lbf",
                               "Reference engine weight")
            Weng = Variable("W_{eng}", "lbf", "engine weight")
            Pslmax = Variable("P_{sl-max}", "hp",
                              "Max shaft power at sea level")

            constraints = [
                Weng/Wengref >= 0.5538*(Pslmax/Pref)**1.075,
                W/mfac >= 2.572*Weng**0.922*units("lbf")**0.078]

        self.dynamic_model = EngineP

        Model.__init__(self, None, constraints, **kwargs)

class EngineP(Model):
    "engine performance model"
    def __init__(self, static, state, **kwargs):

        P_shaft = Variable("P_{shaft}", "hp", "Shaft power")
        bsfc = Variable("BSFC", "lb/hr/hp", "Brake specific fuel consumption")
        rpm = Variable("RPM", "rpm", "Engine operating RPM")
        P_avn = Variable("P_{avn}", 40, "watts", "Avionics power")
        P_pay = Variable("P_{pay}", 10, "watts", "Payload power")
        Ptotal = Variable("P_{total}", "hp", "Total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        lfac = [1 - 0.906**(1/0.15)*(v.value/hr.value)**0.92
                for v, hr in zip(state["h"], state["h_{ref}"])]
        h_loss = Variable("h_{loss}", lfac, "-", "shaft power loss factor")
        Pshaftmax = Variable("P_{shaft-max}",
                             "hp", "Max shaft power at altitude")
        mfac = Variable("m_{fac}", 1.0, "-", "BSFC margin factor")

        if static.DF70:
            rpm_max = Variable("RPM_{max}", 7698, "rpm", "Maximum RPM")
            bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr",
                                "Minimum BSFC")

            constraints = [
                (bsfc/mfac/bsfc_min)**35.7 >=
                (2.29*(rpm/rpm_max)**8.02 + 0.00114*(rpm/rpm_max)**-38.3),
                (Ptotal/Pshaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.294,
                ]
        else:
            bsfc_min = Variable("BSFC_{min}", 0.32, "kg/kW/hr",
                                "Minimum BSFC")
            rpm_max = Variable("RPM_{max}", 9000, "rpm", "Maximum RPM")

            constraints = [
                (bsfc/mfac/bsfc_min)**0.129 >=
                (0.972*(rpm/rpm_max)**-0.141 + 0.0268*(rpm/rpm_max)**9.62),
                (Ptotal/Pshaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.292,
                ]

        constraints.extend([Pshaftmax/static["P_{sl-max}"] == h_loss,
                            Pshaftmax >= Ptotal,
                            rpm <= rpm_max,
                           ])

        if state.onStation:
            constraints.extend([
                Ptotal >= P_shaft + (P_avn + P_pay)/eta_alternator])
        else:
            constraints.extend([Ptotal >= P_shaft + P_avn/eta_alternator])


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
        Cd = Variable("C_d", "-", "wing drag coefficient")
        CL = Variable("C_L", "-", "lift coefficient")
        e = Variable("e", 0.9, "-", "Oswald efficiency")
        Re = Variable("Re", "-", "Reynold's number")
        cdp = Variable("c_{dp}", "-", "wing profile drag coeff")
        Reref = Variable("Re_{ref}", 3e5, "-", "Reference Re for cdp")

        constraints = [
            Cd >= (cdp + CL**2/np.pi/static["A"]/e),
            (cdp/(Re/Reref)**-0.4)**0.00544 >= (
                0.33*CL**-0.0809 + 0.645*CL**0.045 + 7.35e-5*CL**12),
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

        loiter = FlightSegment(N, JHO)

        loiter.substitutions["V"] = np.linspace(20, 40, N)
        mtow = Variable("MTOW", "lbf", "max-take off weight")
        Wfueltot = Variable("W_{fuel-tot}", "lbf", "total fuel weight")

        constraints = [mtow >= loiter["W_{start}"][0],
                       mtow >= JHO["W_{zfw}"] + Wfueltot,
                       Wfueltot >= loiter["W_{fuel-fs}"],
                       loiter["W_{end}"][-1] >= JHO["W_{zfw}"],
                      ]

        Model.__init__(self, mtow, [JHO, loiter, constraints], **kwargs)


if __name__ == "__main__":
    M = Mission()
    # JHO.debug(solver="mosek")
    sol = M.solve("mosek")
    print sol.table()
