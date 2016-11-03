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
        self.smeared_loads = [self.fuse, self.engine]

        Wzfw = Variable("W_{zfw}", "lbf", "zero fuel weight")
        Wpay = Variable("W_{pay}", 10, "lbf", "payload weight")
        Wavn = Variable("W_{avn}", 8, "lbf", "avionics weight")

        constraints = [Wzfw >= sum(summing_vars(self.components, "W")) + Wpay]

        Model.__init__(self, None, [self.components, constraints], **kwargs)

def summing_vars(models, varname):
    "take models, a variable name a returns a list of variables with shared vname"
    modelnames = [m.__class__.__name__ for m in models]
    vkeys = np.hstack([list(m.varkeys[varname]) for m in models])
    vkeys = [v for v in vkeys if v.models[-1] in modelnames]
    vrs = [m[v] for m, v in zip(models, vkeys)]
    return vrs

class AircraftP(Model):
    "performance model for aircraft"
    def __init__(self, static, state, **kwargs):

        self.dynamicmodels = []
        self.dragcomps = []
        self.dragmodels = []
        for c in static.components:
            if "dynamic_model" in c.__dict__.keys():
                dm = c.dynamic_model(c, state)
                self.dynamicmodels.append(dm)
                if "C_d" or "C_f" in dm.varkeys:
                    self.dragmodels.append(dm)
                    self.dragcomps.append(c)

        Wend = Variable("W_{end}", "lbf", "vector-end weight")
        Wstart = Variable("W_{start}", "lbf", "vector-begin weight")
        CD = Variable("C_D", "-", "drag coefficient")
        mfac = Variable("m_{fac}", 1.7, "-", "drag margin factor")

        dvars = []
        for dc, dm in zip(self.dragcomps, self.dragmodels):
            if "C_d" in dm.varkeys:
                dvars.append(dm["C_d"])
            elif "C_f" in dm.varkeys:
                dvars.append(dm["C_f"]*dc["S"]/static.wing["S"])

        constraints = [Wend == Wend,
                       Wstart == Wstart,
                       CD/mfac >= sum(dvars)]

        Model.__init__(self, None, [self.dynamicmodels, constraints], **kwargs)

class FlightState(Model):
    "One chunk of a mission"
    def __init__(self, alt, onStation, wind, **kwargs):

        self.onStation = onStation

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
                       h == h,
                       href == href]

        V = Variable("V", "m/s", "true airspeed")

        if wind:

            V_wind = Variable("V_{wind}", 25, "m/s", "Wind speed")
            constraints.extend([V >= V_wind])

        else:

            V_wind = Variable("V_{wind}", "m/s", "Wind speed")
            V_ref = Variable("V_{ref}", 25, "m/s", "Reference wind speed")

            constraints.extend([(V_wind/V_ref) >= 0.6462*(h/href) + 0.3538,
                                V >= V_wind])
        Model.__init__(self, None, constraints, **kwargs)


class FlightSegment(Model):
    "creates flight segment for aircraft"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 **kwargs):

        with vectorize(N):
            self.fs = FlightState(alt, onStation, wind)
            self.aircraftP = aircraft.dynamic_model(aircraft, self.fs)
            self.slf = SteadyLevelFlight(self.fs, aircraft, self.aircraftP)
            self.be = BreguetEndurance(self.aircraftP)

        self.submodels = [self.fs, self.aircraftP, self.slf, self.be]

        Wfuelfs = Variable("W_{fuel-fs}", "lbf", "flight segment fuel weight")

        self.constraints = [Wfuelfs >= self.be["W_{fuel}"].sum()]

        if N > 1:
            self.constraints.extend([self.aircraftP["W_{end}"][:-1] >=
                                     self.aircraftP["W_{start}"][1:]])

        Model.__init__(self, None, [aircraft, self.submodels,
                                    self.constraints], **kwargs)

class Loiter(FlightSegment):
    "make a loiter flight segment"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 **kwargs):
        super(Loiter, self).__init__(N, aircraft, alt, onStation, wind,
                                     **kwargs)

        t = Variable("t", 6, "days", "time loitering")
        self.constraints.extend([self.be["t"] >= t/N])

        Model.__init__(self, None, [aircraft, self.submodels,
                                    self.constraints], **kwargs)

class Cruise(FlightSegment):
    "make a cruise flight segment"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 R=200, **kwargs):
        super(Cruise, self).__init__(N, aircraft, alt, onStation, wind,
                                     **kwargs)

        R = Variable("R", R, "nautical_miles", "Range to station")
        self.constraints.extend([R/N <= self.aircraftP["V"]*self.be["t"]])

        Model.__init__(self, None, [aircraft, self.submodels,
                                    self.constraints], **kwargs)

class Climb(FlightSegment):
    "make a climb flight segment"
    def __init__(self, N, aircraft, alt=15000, onStation=False, wind=False,
                 dh=15000, **kwargs):
        super(Climb, self).__init__(N, aircraft, alt, onStation, wind,
                                    **kwargs)

        with vectorize(N):
            hdot = Variable("\\dot{h}", "ft/min", "Climb rate")

        deltah = Variable("\\Delta_h", dh, "ft", "altitude difference")
        hdotmin = Variable("\\dot{h}_{min}", 100, "ft/min",
                           "minimum climb rate")

        self.constraints.extend([
            hdot*self.be["t"] >= deltah/N,
            hdot >= hdotmin,
            self.slf["T"] >=
            (0.5*self.fs["\\rho"]*self.fs["V"]**2*self.aircraftP["C_D"]
             * aircraft.wing["S"] + self.aircraftP["W_{start}"]*hdot
             / self.aircraftP["V"]),
            ])

        Model.__init__(self, None, [aircraft, self.submodels,
                                    self.constraints], **kwargs)

class BreguetEndurance(Model):
    "breguet endurance model"
    def __init__(self, perf, **kwargs):
        z_bre = Variable("z_{bre}", "-", "Breguet coefficient")
        t = Variable("t", "days", "Time per flight segment")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        Wfuel = Variable("W_{fuel}", "lbf", "Segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

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
            T >= (0.5*state["\\rho"]*state["V"]**2*perf["C_D"]
                  *aircraft.wing["S"]),
            perf["P_{shaft}"] >= T*state["V"]/etaprop]

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

        Pshaft = Variable("P_{shaft}", "hp", "Shaft power")
        bsfc = Variable("BSFC", "lb/hr/hp", "Brake specific fuel consumption")
        rpm = Variable("RPM", "rpm", "Engine operating RPM")
        Pavn = Variable("P_{avn}", 40, "watts", "Avionics power")
        Ppay = Variable("P_{pay}", 10, "watts", "Payload power")
        Ptotal = Variable("P_{total}", "hp", "Total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        lfac = [1 - 0.906**(1/0.15)*(v.value/hr.value)**0.92
                for v, hr in zip(state["h"], state["h_{ref}"])]
        Leng = Variable("L_{eng}", lfac, "-", "shaft power loss factor")
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
                (Ptotal/Pshaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.294
                ]
        else:
            bsfc_min = Variable("BSFC_{min}", 0.32, "kg/kW/hr",
                                "Minimum BSFC")
            rpm_max = Variable("RPM_{max}", 9000, "rpm", "Maximum RPM")

            constraints = [
                (bsfc/mfac/bsfc_min)**0.129 >=
                (0.972*(rpm/rpm_max)**-0.141 + 0.0268*(rpm/rpm_max)**9.62),
                (Ptotal/Pshaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.292
                ]

        constraints.extend([Pshaftmax/static["P_{sl-max}"] == Leng,
                            Pshaftmax >= Ptotal,
                            rpm <= rpm_max])

        if state.onStation:
            constraints.extend([
                Ptotal >= Pshaft + (Pavn + Ppay)/eta_alternator])
        else:
            constraints.extend([Ptotal >= Pshaft + Pavn/eta_alternator])


        Model.__init__(self, None, constraints, **kwargs)

class Wing(Model):
    "The thing that creates the lift"
    def __init__(self, N=5, lam=0.5, **kwargs):
        self.N = N

        W = Variable("W", "lbf", "weight")
        S = Variable("S", "ft^2", "surface area")
        A = Variable("A", "-", "aspect ratio")
        b = Variable("b", "ft", "wing span")
        tau = Variable("\\tau", 0.115, "-", "airfoil thickness ratio")
        croot = Variable("c_{root}", "ft", "root chord")
        cmac = Variable("c_{MAC}", "ft", "mean aerodynamic chord")
        cb = c_bar(lam, N)
        with vectorize(N):
            cbar = Variable("\\bar{c}", cb, "-",
                            "normalized chord at mid element")
        with vectorize(N-1):
            cave = Variable("c_{ave}", "ft", "mid section chord")

        self.dynamic_model = WingP

        constraints = [b**2 == S*A,
                       tau == tau,
                       cbar == cbar,
                       cave == (cb[1:] + cb[:-1])/2*S/b,
                       croot == S/b*cb[0],
                       cmac == S/b]

        Model.__init__(self, None, constraints, **kwargs)

        capspar = CapSpar(self)
        wingskin = WingSkin(self)
        winginterior = WingInterior(self)
        self.components = [capspar, wingskin, winginterior]
        loading = WingLoading(self.components)

        constraints.extend([W >= sum(c["W"] for c in self.components)])

        Model.__init__(self, None, [self.components, constraints, loading],
                       **kwargs)

def c_bar(lam, N):
    "returns wing chord lengths for constant taper wing"
    eta = np.linspace(0, 1, N)
    c = 2/(1+lam)*(1+(lam-1)*eta)
    return c

class WingLoading(Model):
    "wing loading cases"
    def __init__(self, components, **kwargs):
        loadingcases = []
        for c in components:
            if "loading_model" in c.__dict__.keys():
                loadingcases.append(c.loading_model(c))

        Model.__init__(self, None, loadingcases, **kwargs)

class WingInterior(Model):
    "wing interior model"
    def __init__(self, wing, **kwargs):

        W = Variable("W", "lbf", "interior mass of wing")
        rhofoam = Variable("\\rho_{foam}", 0.036, "g/cm^3", "foam density")
        Abar = Variable("\\bar{A}_{jh01}", 0.0753449, "-",
                        "jh01 non dimensional area")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")

        constraints = [
            W >= (g*rhofoam*Abar*(wing["c_{ave}"])**2*(wing["b"]/2)
                  / (wing.N-1)).sum()]

        Model.__init__(self, None, constraints, **kwargs)

class WingSkin(Model):
    "wing skin model"
    def __init__(self, wing, **kwargs):

        rhocfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        W = Variable("W", "lbf", "wing skin weight")
        g = Variable("g", 9.81, "m/s^2", "gravitational acceleration")
        t = Variable("t", "in", "wing skin thickness")
        tmin = Variable("t_{min}", 0.012, "in",
                        "minimum gague wing skin thickness")
        Jtbar = Variable("\\bar{J/t}", 0.01114, "1/mm",
                         "torsional moment of inertia")

        self.loading_model = WingSkinL

        constraints = [W >= rhocfrp*wing["S"]*2*t*g,
                       t >= tmin,
                       Jtbar == Jtbar,
                       wing["c_{root}"] == wing["c_{root}"]]

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
            taucfrp >= (1/static["\\bar{J/t}"]/static["c_{root}"]**2/static["t"]
                        * Cmw*static["S"]*rhosl*Vne**2)]

        Model.__init__(self, None, constraints, **kwargs)

class CapSpar(Model):
    "cap spar model"
    def __init__(self, wing, **kwargs):
        self.N = wing.N

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

        self.loading_model = CapSparL

        constraints = [I <= 2*w*t*(hin/2)**2,
                       dm >= rhocfrp*w*t*wing["b"]/(self.N-1),
                       W >= dm.sum()*g,
                       w <= w_lim*wing["c_{ave}"],
                       wing["c_{ave}"]*wing["\\tau"] >= hin + 2*t,
                       E == E,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class CapSparL(Model):
    def __init__(self, static, **kwargs):

        Nmax = Variable("N_{max}", 5, "-", "max loading")
        Wcent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        cbar = c_bar(0.5, 5)
        sigmacfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")
        kappa = Variable("\\kappa", 0.2, "-", "max tip deflection ratio")
        Mroot = Variable("M_{root}", "N*m", "wing root moment")

        beam = Beam(static.N, cbar)

        constraints = [
            # dimensionalize moment of inertia and young's modulus
            beam["\\bar{EI}"] <= (8*static["E"]*static["I"]/Nmax
                                  / Wcent/static["b"]**2),
            Mroot == (beam["\\bar{M}"][0]*Wcent*Nmax
                      * static["b"]/4),
            sigmacfrp >= Mroot*(static["h_{in}"]+static["t"])/static["I"],
            beam["\\bar{\\delta}"][-1] <= kappa,
            ]

        Model.__init__(self, None, [beam, constraints], **kwargs)

class Beam(Model):
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


        Sbar_tip = Variable("\\bar{S}_{tip}", 1e-10, "-", "Tip loading")
        Mbar_tip = Variable("\\bar{M}_{tip}", 1e-10, "-", "Tip moment")
        th_root = Variable("\\theta_{root}", 1e-10, "-", "Base angle")
        dbar_root = Variable("\\bar{\\delta}_{root}", 1e-10, "-",
                             "Base deflection")
        dx = Variable("dx", "-", "normalized length of element")

        constraints = [
            Sbar[:-1] >= Sbar[1:] + 0.5*dx*(qbar[:-1] + qbar[1:]),
            Sbar[-1] >= Sbar_tip,
            Mbar[:-1] >= Mbar[1:] + 0.5*dx*(Sbar[:-1] + Sbar[1:]),
            Mbar[-1] >= Mbar_tip,
            th[0] >= th_root,
            th[1:] >= th[:-1] + 0.5*dx*(Mbar[1:] + Mbar[:-1])/EIbar,
            dbar[0] >= dbar_root,
            dbar[1:] >= dbar[:-1] + 0.5*dx*(th[1:] + th[:-1]),
            1 == (N-1)*dx,
            ]

        Model.__init__(self, None, constraints, **kwargs)


class WingP(Model):
    def __init__(self, static, state, **kwargs):
        "wing drag model"
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
            Re == state["\\rho"]*state["V"]*static["S"]/static["b"]/state["\\mu"],
            ]
        Model.__init__(self, None, constraints, **kwargs)

class FuelTank(Model):
    """
    Returns the weight of the fuel tank.  Assumes a cylinder shape with some
    fineness ratio
    """
    def __init__(self, **kwargs):

        phi = Variable("\\phi", 6, "-", "fuel tank fineness ratio")
        l = Variable("l", "ft", "fuel tank length")
        Stank = Variable("S_{tank}", "ft^2", "fuel tank surface area")
        W = Variable("W", "lbf", "fuel tank weight")
        Wfueltot = Variable("W_{fuel-tot}", "lbf", "Total fuel weight")
        m_fac = Variable("m_{fac}", 1.1, "-", "fuel volume margin factor")
        rhofuel = Variable("\\rho_{fuel}", 6.01, "lbf/gallon",
                           "density of 100LL")
        rhotank = Variable("\\rho_{fuel-tank}", 0.089, "g/cm^2",
                           "density of plastic fuel tank")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        Voltank = Variable("\\mathcal{V}", "ft^3", "fuel tank volume")

        constraints = [W >= Stank*rhotank*g,
                       Stank/4/phi >= Voltank/l,
                       Voltank/m_fac >= Wfueltot/rhofuel,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Fuselage(Model):
    "The thing that carries the fuel, engine, and payload"
    def __init__(self, **kwargs):
        d = Variable("d", "ft", "fuselage diameter")
        lfuse = Variable("l_{fuse}", "ft", "fuselage length")

        mskin = Variable("m_{skin}", "kg", "fuselage skin mass")
        rhokevlar = Variable("\\rho_{kevlar}", 1.3629, "g/cm**3",
                             "kevlar density")
        Sfuse = Variable("S", "ft^2", "Fuselage surface area")
        Volavn = Variable("\\mathcal{V}_{avn}", 0.125, "ft^3",
                          "Avionics volume")
        W = Variable("W", "lbf", "Fuselage weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        m_fac = Variable("m_{fac}", 2.5, "-", "Fuselage weight margin factor")
        tmin = Variable("t_{min}", 0.03, "in", "minimum skin thickness")
        tskin = Variable("t_{skin}", "in", "skin thickness")
        hengine = Variable("h_{engine}", 6, "in", "engine height")

        ft = FuelTank()
        self.dynamic_model = FuselageP

        constraints = [
            mskin >= Sfuse*rhokevlar*tskin,
            tskin >= tmin,
            Sfuse >= np.pi*d*lfuse + np.pi*d**2,
            np.pi*(d/2)**2*lfuse >= ft["\\mathcal{V}"] + Volavn,
            lfuse >= ft["l"],
            d >= hengine,
            W/m_fac >= mskin*g + ft["W"],
            ]

        Model.__init__(self, None, [constraints, ft], **kwargs)

class FuselageP(Model):
    "fuselage drag model"
    def __init__(self, static, state, **kwargs):

        CDA = Variable("CDA", "-", "fuselage area drag normalized by wing area")
        Cf = Variable("C_f", "-", "fuselage skin friction coefficient")
        Re = Variable("Re", "-", "fuselage reynolds number")

        constraints = [
            Re == state["V"]*state["\\rho"]*static["l_{fuse}"]/state["\\mu"],
            Cf >= 0.455/Re**0.3
            ]

        Model.__init__(self, None, constraints, **kwargs)

class Mission(Model):
    "creates flight profile"
    def __init__(self, **kwargs):
        JHO = Aircraft()

        climb1 = Climb(10, JHO, alt=np.linspace(0, 15000, 11)[1:])
        cruise1 = Cruise(1, JHO)
        loiter1 = Loiter(5, JHO)
        cruise2 = Cruise(1, JHO)
        mission = [climb1, cruise1, loiter1, cruise2]

        mtow = Variable("MTOW", "lbf", "max-take off weight")

        constraints = [
            mtow >= mission[0]["W_{start}"][0],
            mtow >= JHO["W_{zfw}"] + JHO["W_{fuel-tot}"],
            JHO["W_{fuel-tot}"] >= sum(fs["W_{fuel-fs}"] for fs in mission),
            mission[-1]["W_{end}"][-1] >= JHO["W_{zfw}"],
            JHO["W_{cent}"] >= (JHO["W_{fuel-tot}"] +
                                sum(summing_vars(JHO.smeared_loads, "W")))
            ]

        for i, fs in enumerate(mission[1:]):
            constraints.extend([
                mission[i]["W_{end}"][-1] == fs["W_{start}"][0]
                ])

        Model.__init__(self, mtow, [JHO, mission, constraints], **kwargs)


if __name__ == "__main__":
    M = Mission()
    # JHO.debug(solver="mosek")
    sol = M.solve("mosek")
    print sol.table()
