""" gas_male_rubber.py """
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit import LinkedConstraintSet, ConstraintSet
from gpkit import SignomialsEnabled
from helpers import SummingConstraintSet
#from gpkit.tools import BoundedConstraintSet
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
PLOT = False
SIGNOMIALS = False

INCLUDE = ["l_{fuse}", "MTOW", "t_{loiter}", "S", "b", "AR", "W_{zfw}",
           "P_{shaft-maxMSL}", "S_{fuse}", "W_{cent}", "W_{fuel-tot}", "g",
           "V_{stall}"]

class Mission(Model):
    def __init__(self, h_station, wind, DF70, Nclimb, Nloiter, **kwargs):

        self.submodels = [
            TakeOff(1, [0.684], [0.1], False, wind, DF70),
            Climb(Nclimb, [0.502]*Nclimb, np.linspace(0, 15000, Nclimb+1)[1:],
                  False, wind, DF70, dh=15000),
            Cruise(1, [0.684], [15000], False, wind, DF70, R=180),
            Loiter(Nloiter, [0.647]*Nloiter, [h_station]*Nloiter, True, wind,
                   DF70),
            Cruise(1, [0.684], [15000], False, wind, DF70, R=200)
            ]

        mtow = Variable("MTOW", "lbf", "Max take off weight")
        W_zfw = Variable("W_{zfw}", "lbf", "Zero fuel weight")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "Total fuel weight")
        m_fac = Variable("m_{fac}", 1.0, "-", "MTOW margin factor")

        constraints = [
            mtow/m_fac >= self.submodels[0]["W_{start}"],
            W_zfw <= self.submodels[-1]["W_{end}"],
            W_fueltot >= sum(fs["W_{fuel-fs}"] for fs in self.submodels)
            ]

        for i, fs in enumerate(self.submodels[1:]):
            constraints.extend([
                self.submodels[i]["W_{end}"] == fs["W_{start}"]
                ])

        lc = LinkedConstraintSet(
            [fs for fs in self.submodels, constraints],
            include_only=INCLUDE + ["W"])

        Model.__init__(self, None, lc, **kwargs)

class FlightSegment(Model):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, **kwargs):

        self.aero = Aerodynamics(N)
        self.fuel = Fuel(N)
        self.slf = SteadyLevelFlight(N, eta_p)
        self.engine = EnginePerformance(N, alt, onStation, DF70)
        self.atm = Atmosphere(N, alt)
        self.wind = Wind(N, alt, wind)

        self.N = N
        self.exclude = ["m_{fac}"]

        # self.Vstall = VectorVariable(N, "V_{stall}", "m/s", "stall speed")

        # self.constraints = [
        #     self.Vstall == ((self.fuel["W_{N}"]*2/self.slf["\\rho"]/
        #                      self.aero["S"]/1.5)**0.5)
        #     ]
        self.constraints = []

        self.submodels = [self.aero, self.fuel, self.slf, self.engine,
                          self.atm, self.wind]

class TakeOff(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        self.submodels.extend([breguetendurance])

        t_to = Variable("t_{TO}", 10, "minutes", "take off time")
        Vstall = Variable("V_{stall}", "m/s", "stall speed")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")

        self.constraints = [breguetendurance["t"] >= t_to,
                            self.slf["V"] >= 1.3*Vstall,
                            Vstall >= (self.fuel["W_{N}"][0]*2/rhosl/self.aero["S"]/1.5)**0.5
                           ]

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 exclude=self.exclude)

        Model.__init__(self, None, lc, **kwargs)

class Cruise(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, R=200, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        R = Variable("R", R, "nautical_miles", "Range to station")

        self.submodels.extend([breguetendurance])

        self.constraints.extend([R/N <= self.slf["V"]*breguetendurance["t"]])

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 exclude=self.exclude)

        Model.__init__(self, None, lc, **kwargs)

class Loiter(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        t_loiter = Variable("t_{loiter}", "days", "time loitering")

        self.constraints.extend([breguetendurance["t"] >= t_loiter/N])

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 exclude=self.exclude)

        Model.__init__(self, None, lc, **kwargs)

class Climb(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, dh=5000,
                 **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        deltah = Variable("\\delta_h", dh, "ft", "altitude difference")
        h_dot = VectorVariable(N, "h_{dot}", "ft/min", "Climb rate")
        h_dotmin = Variable("h_{dot-min}", 100, "ft/min",
                            "minimum climb rate")
        self.constraints.extend([
            h_dot*breguetendurance["t"] >= deltah/N,
            h_dot >= h_dotmin,
            self.slf["T"] >= (0.5*self.slf["\\rho"]*self.slf["V"]**2*
                              self.slf["C_D"]*self.slf["S"] +
                              self.slf["W_{N}"]*h_dot/self.slf["V"])
            ])

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 exclude=self.exclude)

        Model.__init__(self, None, lc, **kwargs)

class Atmosphere(Model):
    """
    Model to capture density changes with altitude
    """
    def __init__(self, N, alt, **kwargs):

        h = VectorVariable(N, "h", alt, "ft", "Altitude")
        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")

        T_atm = VectorVariable(N, "T_{atm}",
                               [T_sl.value - L_atm.value*v.value for v in h],
                               "K", "Air temperature")
        mu_atm = VectorVariable(N, "\\mu", "N*s/m^2", "Dynamic viscosity")
        mu_sl = Variable("\\mu_{sl}", 1.789*10**-5, "N*s/m^2",
                         "Dynamic viscosity at sea level")
        R_spec = Variable("R_{spec}", 287.058, "J/kg/K",
                          "Specific gas constant of air")
        h_ref = Variable("h_{ref}", 15000, "ft", "Reference altitude")
        rho = VectorVariable(N, "\\rho", "kg/m^3", "Air density")

        # Atmospheric variation with altitude (valid from 0-7km of altitude)
        constraints = [rho == p_sl*T_atm**(5.257-1)/R_spec/(T_sl**5.257),
                       (mu_atm/mu_sl)**0.1 == 0.991*(h/h_ref)**(-0.00529)
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Fuel(Model):
    """
    Fuel weight model
    """
    def __init__(self, N, **kwargs):

        #----------------------------------------------------
        # Fuel weight model

        W_start = Variable("W_{start}", "lbf",
                           "weight at beginning of flight segment")
        W_nplus1 = VectorVariable(N, "W_{N+1}", "lbf", "vector-end weight")
        W_fuel = VectorVariable(N, "W_{fuel}", "lbf",
                                "Segment-fuel weight")
        W_fuelfs = Variable("W_{fuel-fs}", "lbf",
                            "flight segment fuel weight")
        W_end = Variable("W_{end}", "lbf",
                         "weight at beginning of flight segment")
        W_n = VectorVariable(N, "W_{N}", "lbf", "vector-begin weight")

        # end of first segment weight + first segment fuel weight must be
        # greater  than MTOW.  Each end of segment weight must be greater
        # than the next end of segment weight + the next segment fuel weight.
        # The last end segment weight must be greater than the zero fuel
        # weight
        constraints = [W_start >= W_nplus1[0] + W_fuel[0],
                       W_fuelfs >= W_fuel.sum(),
                       W_nplus1[-1] >= W_end,
                       W_n[0] == W_start]

        if N == 1:
            pass
        else:
            constraints.extend([
                W_nplus1[:-1] >= W_nplus1[1:] + W_fuel[1:],
                W_n[1:] == W_nplus1[:-1],
                ])

        Model.__init__(self, None, constraints, **kwargs)

class SteadyLevelFlight(Model):
    """
    Captures steady level flight mode
    """
    def __init__(self, N, eta_p, **kwargs):

        CD = VectorVariable(N, "C_D", "-", "Drag coefficient")
        CL = VectorVariable(N, "C_L", "-", "Lift coefficient")
        V = VectorVariable(N, "V", "m/s", "Cruise speed")
        S = Variable("S", "ft^2", "wing area")
        eta_prop = VectorVariable(N, "\\eta_{prop}", eta_p, "-",
                                  "Propulsive efficiency")
        P_shaft = VectorVariable(N, "P_{shaft}", "hp", "Shaft power")
        T = VectorVariable(N, "T", "lbf", "Thrust")

        rho = VectorVariable(N, "\\rho", "kg/m^3", "Air density")
        W_nplus1 = VectorVariable(N, "W_{N+1}", "lbf", "vector-end weight")
        W_n = VectorVariable(N, "W_{N}", "lbf", "vector-begin weight")

        # Climb model
        # Currently climb rate is being optimized to reduce fuel consumption.
        # In future, could implement min climb rate.

        constraints = [
            P_shaft == T*V/eta_prop,
            T >= 0.5*rho*V**2*CD*S,
            0.5*rho*CL*S*V**2 == (W_nplus1*W_n)**0.5,
            ]
        # Propulsive efficiency variation with different flight segments,
        # will change depending on propeller characteristics

        Model.__init__(self, None, constraints, **kwargs)

class EnginePerformance(Model):
    """
    Engine performance and weight model for small engine
    """
    def __init__(self, N, alt, onStation, DF70, **kwargs):

        h = VectorVariable(N, "h", alt, "ft", "Altitude")
        h_ref = Variable("h_{ref}", 15000, "ft", "Reference altitude")
        P_shaft = VectorVariable(N, "P_{shaft}", "hp", "Shaft power")
        bsfc = VectorVariable(N, "BSFC", "lb/hr/hp",
                              "Brake specific fuel consumption")
        rpm = VectorVariable(N, "RPM", "rpm", "Engine operating RPM")
        P_avn = Variable("P_{avn}", 40, "watts", "Avionics power")
        P_pay = Variable("P_{pay}", 10, "watts", "Payload power")
        P_shafttot = VectorVariable(N, "P_{shaft-tot}", "hp",
                                    "Total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        lfac = [1 - 0.906**(1/0.15)*(v.value/h_ref.value)**0.92
                for v in h]
        h_loss = VectorVariable(N, "h_{loss}", lfac, "-",
                                "Max shaft power loss factor")
        P_shaftmax = VectorVariable(N, "P_{shaft-max}",
                                    "hp", "Max shaft power at altitude")
        m_fac = Variable("m_{fac}", 1.0, "-", "BSFC margin factor")

        if DF70:
            P_shaftmaxmsl = Variable("P_{shaft-maxMSL}", 5.17, "hp",
                                     "Max shaft power at MSL")
            rpm_max = Variable("RPM_{max}", 7698, "rpm", "Maximum RPM")
            bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr",
                                "Minimum BSFC")

            constraints = [
                (bsfc/m_fac/bsfc_min)**35.7 >=
                (2.29*(rpm/rpm_max)**8.02 + 0.00114*(rpm/rpm_max)**-38.3),
                (P_shafttot/P_shaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.294,
                ]
        else:
            bsfc_min = Variable("BSFC_{min}", 0.32, "kg/kW/hr",
                                "Minimum BSFC")
            rpm_max = Variable("RPM_{max}", 9000, "rpm", "Maximum RPM")
            P_shaftmaxmsl = Variable("P_{shaft-maxMSL}", "hp",
                                     "Max shaft power at MSL")

            constraints = [
                (bsfc/m_fac/bsfc_min)**0.129 >=
                (0.972*(rpm/rpm_max)**-0.141 + 0.0268*(rpm/rpm_max)**9.62),
                (P_shafttot/P_shaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.292,
                ]

        constraints.extend([P_shaftmax/P_shaftmaxmsl == h_loss,
                            P_shaftmax >= P_shafttot,
                            rpm <= rpm_max,
                           ])

        if onStation:
            constraints.extend([
                P_shafttot >= P_shaft + (P_avn + P_pay)/eta_alternator])
        else:
            constraints.extend([P_shafttot >= P_shaft + P_avn/eta_alternator])


        Model.__init__(self, None, constraints, **kwargs)

class BreguetEndurance(Model):
    """
    Discritized Breguet Range model
    """
    def __init__(self, N, **kwargs):
        z_bre = VectorVariable(N, "z_{bre}", "-", "Breguet coefficient")
        t = VectorVariable(N, "t", "days", "Time per flight segment")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        P_shafttot = VectorVariable(N, "P_{shaft-tot}", "hp",
                                    "Total power, avionics included")
        bsfc = VectorVariable(N, "BSFC", "lb/hr/hp",
                              "Brake specific fuel consumption")
        W_nplus1 = VectorVariable(N, "W_{N+1}", "lbf", "vector-end weight")
        W_fuel = VectorVariable(N, "W_{fuel}", "lbf",
                                "Segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        W_n = VectorVariable(N, "W_{N}", "lbf", "vector-begin weight")

        constraints = [
            z_bre >= P_shafttot*t*bsfc*g/(W_nplus1*W_n)**0.5,
            # TCS([z_bre >= P_shafttot*t*bsfc*g/W_nplus1]),
            f_fueloil*W_fuel/W_nplus1 >= te_exp_minus1(z_bre, 3)
            ]

        Model.__init__(self, None, constraints, **kwargs)

class Aerodynamics(Model):
    """
    Aero model assuming jh01 airfoil, designed by Mark Drela
    """
    def __init__(self, N, **kwargs):

        CLmax = Variable("C_{L-max}", 1.5, "-", "Maximum lift coefficient")
        e = Variable("e", 0.9, "-", "Spanwise efficiency")
        AR = Variable("AR", "-", "Aspect ratio")
        b = Variable("b", "ft", "Span")
        Re = VectorVariable(N, "Re", "-", "Reynolds number")

        # fuselage drag
        Kfuse = Variable("K_{fuse}", 1.1, "-", "Fuselage form factor")
        S_fuse = Variable("S_{fuse}", "ft^2", "Fuselage surface area")
        Cffuse = VectorVariable(N, "C_{f-fuse}", "-",
                                "Fuselage skin friction coefficient")
        CDfuse = VectorVariable(N, "C_{D-fuse}", "-", "fueslage drag")
        l_fuse = Variable("l_{fuse}", "ft", "Fuselage length")
        Refuse = VectorVariable(N, "Re_{fuse}", "-",
                                "Fuselage Reynolds number")
        Re_ref = Variable("Re_{ref}", 3e5, "-", "Reference Re for cdp")
        cdp = VectorVariable(N, "c_{dp}", "-", "wing profile drag coeff")

        CD = VectorVariable(N, "C_D", "-", "Drag coefficient")
        CL = VectorVariable(N, "C_L", "-", "Lift coefficient")
        V = VectorVariable(N, "V", "m/s", "Cruise speed")
        S = Variable("S", "ft^2", "wing area")
        rho = VectorVariable(N, "\\rho", "kg/m^3", "Air density")
        mu_atm = VectorVariable(N, "\\mu", "N*s/m^2", "Dynamic viscosity")
        m_fac = Variable("m_{fac}", 1.0, "-", "cdp margin factor")

        constraints = [
            CD >= CDfuse*2 + cdp*1.3 + CL**2/(pi*e*AR),
            #jh01
            cdp/m_fac >= ((0.0075 + 0.002*CL**2 +
                           0.00033*CL**10)*(Re/Re_ref)**-0.4),
            #sd7032
            # cdp >= ((0.006 + 0.005*CL**2 + 0.00012*CL**10)*(Re/Re_ref)**-0.3),
            b**2 == S*AR,
            CL <= CLmax,
            Re == rho*V/mu_atm*(S/AR)**0.5,
            CDfuse >= Kfuse*S_fuse*Cffuse/S,
            Refuse == rho*V/mu_atm*l_fuse,
            Cffuse >= 0.455/Refuse**0.3,
            ]

        Model.__init__(self, None, constraints, **kwargs)

class LandingGear(Model):
    """
    Landing gear, with fixed dimensions.  Used for preliminary study
    """
    def __init__(self, NSeg, **kwargs):

        A_rearland = Variable("A_{rear-land}", 12, "in^2",
                              "rear landing gear frontal area")
        A_frontland = Variable("A_{front-land}", 18, "in^2",
                               "front landing gear frontal area")
        CDland = Variable("C_{D-land}", 0.2, "-",
                          "drag coefficient landing gear")
        CDAland = Variable("CDA_{land}", "-",
                           "normalized drag coefficient landing gear")

        CD = VectorVariable(NSeg, "C_D", "-", "Drag coefficient")
        CL = VectorVariable(NSeg, "C_L", "-", "Lift coefficient")
        CDfuse = VectorVariable(NSeg, "C_{D-fuse}", "-", "fueslage drag")
        S = Variable("S", "ft^2", "wing area")
        e = Variable("e", 0.95, "-", "Spanwise efficiency")
        AR = Variable("AR", "-", "Aspect ratio")
        cdp = VectorVariable(NSeg, "c_{dp}", "-", "wing profile drag coeff")

        constraints = [
            CD >= CDfuse*2 + cdp*1.3 + CL**2/(pi*e*AR) + CDAland,
            CDAland >= (2*CDland*A_rearland + CDland*A_frontland)/S
            ]

        Model.__init__(self, None, constraints, **kwargs)

class Beam(Model):
    def __init__(self, N, q, untapered, **kwargs):

        qbar = VectorVariable(N, "\\bar{q}", q, "-", "normalized loading")
        Sbar = VectorVariable(N, "\\bar{S}", "-", "normalized shear")
        Sbar_tip = Variable("\\bar{S}_{tip}", 1e-10, "-", "Tip loading")
        Mbar = VectorVariable(N, "\\bar{M}", "-", "normalized moment")
        Mbar_tip = Variable("\\bar{M}_{tip}", 1e-10, "-", "Tip moment")
        th = VectorVariable(N, "\\theta", "-", "deflection slope")
        th_root = Variable("\\theta_{root}", 1e-10, "-", "Base angle")
        dbar = VectorVariable(N, "\\bar{\\delta}", "-",
                              "normalized displacement")
        dbar_root = Variable("\\bar{\\delta}_{root}", 1e-10, "-",
                             "Base deflection")
        dx = Variable("dx", "-", "normalized length of element")
        if untapered:
            EIbar = Variable("\\bar{EI}", "-",
                             "normalized YM and moment of inertia")
        else:
            EIbar = VectorVariable(N-1, "\\bar{EI}", "-",
                                   "normalized YM and moment of inertia")

        constraints = [
            Sbar[:-1] >= Sbar[1:] + 0.5*dx*(qbar[:-1] + qbar[1:]),
            Sbar[-1] >= Sbar_tip,
            Mbar[:-1] >= Mbar[1:] + 0.5*dx*(Sbar[:-1] + Sbar[1:]),
            Mbar[-1] >= Mbar_tip,
            th[0] >= th_root,
            th[1:] >= th[:-1] + 0.5*dx*(Mbar[1:] + Mbar[:-1])/EIbar,
            dbar[0] >= dbar_root,
            dbar[1:] >= dbar[:-1] + 0.5*dx*(th[1:] + th[:-1]),
            1 == (N-1)*dx
            ]

        Model.__init__(self, None, constraints, **kwargs)

def c_bar(lam, N):
    eta = np.linspace(0, 1, N)
    c = 2/(1+lam)*(1+(lam-1)*eta)
    return c

class Spar(Model):
    def __init__(self, N=5, untapered=False, **kwargs):
        # NOTE: untapered does not solve with current values.  Relax w_lim
        # and t_loiter to solve.

        # phyiscal properties
        rho_cfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        rho_fg = Variable("\\rho_{fg}", 0.75, "g/cm^3", "density of fiberglass")
        E = Variable("E", 2e7, "psi", "Youngs modulus of CF")
        sigma_cfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")
        sigma_fg = Variable("\\sigma_{fg}", 175e6, "Pa",
                            "Fiberglass max stress")

        # Structural lengths
        cb = c_bar(0.5, N)
        cbavg = (cb[:-1] + cb[1:])/2
        cbar = VectorVariable(N-1, "\\bar{c}", cbavg, "-",
                              "normalized chord at mid element")
        if untapered:
            n = 1
        else:
            n = N-1
        t = VectorVariable(n, "t", "in", "spar cap thickness")
        hin = VectorVariable(n, "h_{in}", "in", "inner spar height")
        w = VectorVariable(n, "w", "in", "spar width")
        ts = VectorVariable(n, "t_s", "in", "shear casing thickness")
        I = VectorVariable(n, "I", "m^4", "spar x moment of inertia")
        dm = VectorVariable(n, "dm", "kg", "segment spar mass")
        m = Variable("m", "kg", "spar mass")

        S = Variable("S", "ft^2", "wing area")
        tau = Variable("\\tau", 0.115, "-", "Airfoil thickness ratio")
        b = Variable("b", "ft", "Span")

        N_max = Variable("N_{max}", 5, "-", "Load factor")
        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")

        kappa = Variable("\\kappa", 0.2, "-", "Max tip deflection ratio")
        w_lim = Variable("w_{lim}", 0.14, "-", "spar width to chord ratio")

        beam = Beam(N, cb, untapered)
        self.submodels = [beam]

        constraints = [
            dm >= rho_cfrp*w*t*b/n + rho_fg*b/2/n*2*ts*(w+hin),
            m >= dm.sum(),
            w <= w_lim*S/b*cbar,
            S/b*cbar*tau >= hin + 2*t + 2*ts,
            beam["\\bar{S}"][:-1]*W_cent*N_max/b*(b/2)/4/hin/ts <= sigma_fg,
            beam["\\bar{\\delta}"][-1] <= kappa,
            sigma_cfrp >= ((beam["\\bar{M}"][:-1] + beam["\\bar{M}"][1:])/
                           2*b*W_cent*N_max/4*(hin+t)/I),
            beam["\\bar{EI}"] <= E*I/N_max/W_cent*b/(b/2)**3
            ]

        if SIGNOMIALS:
            with SignomialsEnabled():
                constraints.extend([I <= w*t**3/6 + 2*w*t*(hin/2+t/2)**2])
        else:
            constraints.extend([I <= 2*w*t*(hin/2)**2])

        lc = LinkedConstraintSet([beam, constraints])

        Model.__init__(self, None, lc, **kwargs)

class Wing(Model):
    """
    Structural wing model.  Simple beam.
    """
    def __init__(self, **kwargs):

        # Structural parameters
        rho_skin = Variable("\\rho_{skin}", 0.1, "g/cm^2",
                            "Wing skin density")

        # wing parameters
        S = Variable("S", "ft^2", "wing area")
        #find better number
        wingloading = Variable("W/S", "lbf/ft^2", "Wing loading")

        # Structural evaluation parameters

        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        mtow = Variable("MTOW", "lbf", "Max take off weight")
        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        W = Variable("W", "lbf", "Total wing structural weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        m_fac = Variable("m_{fac}", 1.0, "-", "Wing weight margin factor")

        self.spar = Spar(5)
        self.submodels = [self.spar]

        constraints = [m_skin >= rho_skin*S*2,
                       wingloading == mtow/S,
                       W/m_fac >= m_skin*g + 2*self.spar["m"]*g,
                      ]

        lc = LinkedConstraintSet([self.spar, constraints])

        Model.__init__(self, None, lc, **kwargs)

class Fuselage(Model):
    """
    Sizes fuselage based off of volume constraints.  Assumes elliptical shape
    """
    def __init__(self, **kwargs):

        # Constants
        rho_fuel = Variable("\\rho_{fuel}", 6.01, "lbf/gallon",
                            "density of 100LL")

        # Non-dimensional variables
        k1fuse = Variable('k_{1-fuse}', 4.3246, '-', 'fuselage form factor 1')
        k2fuse = Variable('k-{2-fuse}', 7.124, '-', 'fuselage form factor 2')
        w_cent = Variable('w_{cent}', 'ft', 'center fuselage width')
        fr = Variable('fr', 6.5, '-', 'fineness ratio fuselage')

        # Volumes
        Vol_fuel = Variable("Vol_{fuel}", "m**3", "Fuel Volume")
        Vol_fuse = Variable("Vol_{fuse}", "m**3", "Fuselage volume")

        m_fuse = Variable("m_{fuse}", "kg", "Fuselage mass")
        rho_skin = Variable("\\rho_{skin}", 0.1, "g/cm^2",
                            "Wing skin density")
        S_fuse = Variable("S_{fuse}", "ft^2", "Fuselage surface area")
        l_fuse = Variable("l_{fuse}", "ft", "Fuselage length")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "Total fuel weight")
        l_cent = Variable("l_{cent}", "ft", "Center fuselage length")
        Vol_avionics = Variable("Vol_{avionics}", 0.125, "ft^3",
                                "Avionics volume")
        Vol_pay = Variable("Vol_{pay}", 1, "ft^3", "Payload volume")
        m_rib = Variable("m_{rib}", 1.36, "kg", "Rib mass")
        m_fuse = Variable("m_{fuse}", "kg", "Fuselage mass")
        W = Variable("W", "lbf", "Fuselage weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        m_fac = Variable("m_{fac}", 1.0, "-", "Fuselage weight margin factor")

        constraints = [m_fuse >= S_fuse*rho_skin,
                       l_cent == fr*w_cent,
                       l_fuse >= l_cent*1.1,
                       (l_fuse/k1fuse)**3 == Vol_fuse,
                       (S_fuse/k2fuse)**3 == Vol_fuse**2,
                       Vol_fuse >= l_cent*w_cent**2,
                       Vol_fuel >= W_fueltot/rho_fuel,
                       l_cent*w_cent**2 >= Vol_fuel+Vol_avionics+Vol_pay,
                       W/m_fac >= m_fuse*g + m_rib*g
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Wind(Model):
    """
    Model for wind speed
    wind = True, wind speed has specific value
    """
    def __init__(self, N, alt, wind, **kwargs):

        V = VectorVariable(N, "V", "m/s", "Cruise speed")
        h = VectorVariable(N, "h", alt, "ft", "Altitude")

        if wind:

            V_wind = Variable("V_{wind}", 25, "m/s", "Wind speed")
            constraints = [V >= V_wind]

        else:

            V_wind = VectorVariable(N, "V_{wind}", "m/s", "Wind speed")
            h_ref = Variable("h_{ref}", 15000, "ft", "Reference altitude")
            V_ref = Variable("V_{ref}", 25, "m/s", "Reference wind speed")

            constraints = [(V_wind/V_ref) >= 0.6462*(h/h_ref) + 0.3538,
                           V >= V_wind,
                          ]

        Model.__init__(self, None, constraints, **kwargs)

class EngineWeight(Model):
    def __init__(self, DF70, **kwargs):

        W = Variable("W", "lbf", "Installed/Total engine weight")
        m_fac = Variable("m_{fac}", 1.0, "-", "Engine weight margin factor")

        if DF70:
            W_df70 = Variable("W_{DF70}", 7.1, "lbf",
                              "Installed/Total DF70 engine weight")
            P_shaftmaxmsl = Variable("P_{shaft-maxMSL}", 5.17, "hp",
                                     "Max shaft power at MSL")
            constraints = [W/m_fac >= W_df70]

        else:
            P_shaftref = Variable("P_{shaft-ref}", 2.295, "hp",
                                  "Reference shaft power")
            W_engref = Variable("W_{eng-ref}", 4.4107, "lbf",
                                "Reference engine weight")
            W_eng = Variable("W_{eng}", "lbf", "engine weight")
            P_shaftmaxmsl = Variable("P_{shaft-maxMSL}", "hp",
                                     "Max shaft power at MSL")

            constraints = [
                W_eng/W_engref >= 0.5538*(P_shaftmaxmsl/P_shaftref)**1.075,
                W/m_fac >= 2.572*W_eng**0.922*units("lbf")**0.078]

        Model.__init__(self, None, constraints, **kwargs)

class TailBoom(Model):
    def __init__(self, **kwargs):

        F = Variable("F", "N", "point force from tail")
        L = Variable("L", 5, "ft", "tail boom length")
        E = Variable("E", 150e9, "N/m^2", "young's modulus carbon fiber")
        kfac = Variable("(1-k/2)", 0.6, "-", "(1-k/2) tail boom inertia value")
        I0 = Variable("I_0", "m^4", "tail boom moment of inertia")
        d0 = Variable("d_0", "m", "tail boom diameter")
        t0 = Variable("t_0", "mm", "tail boom thickness")
        tmin = Variable("t_{min}", 0.25, "mm", "minimum tail boom thickness")
        rho_cfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        W = Variable("W", "lbf", "tail boom weight")
        th = Variable("\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.5, "-",
                         "max tail boom deflection angle")

        constraints = [I0 <= pi*t0*d0**3/8.0,
                       W >= pi*g*rho_cfrp*d0*L*t0*kfac,
                       t0 >= tmin,
                       th <= thmax,
                       th >= F*L**2/E/I0*kfac
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Tail(Model):
    def __init__(self, case, **kwargs):
        W_vtail = Variable("W_{v-tail}", 3.4999, "lbf", "V-Tail weight")
        m_fac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "Tail weight")
        Sh = Variable("S_h", "ft**2", "horizontal tail area")
        CLh = Variable("C_{L_h}", 1.0, "-", "maximum lift of horizontal tail")
        CLw = Variable("C_{L_w}", 1.5, "-", "maximum lift of win")
        Vh = Variable("V_h", "-", "horizontal tail volume coefficient")
        S = Variable("S", "ft^2", "wing area")
        b = Variable("b", "ft", "Span")
        Cmw = Variable("C_{m_w}", 0.121, "-", "negative wing moment coefficent")
        xcg = Variable("x_{cg}-x_{ac}", 4, "in", "distance from AC to CG")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Vstall = Variable("V_{stall}", "m/s", "stall speed")

        tb = TailBoom()
        self.submodels = [tb]

        constraints = [
            W/m_fac >= W_vtail + tb["W"],
            tb["F"] >= 0.5*rhosl*Vstall**2*Sh*CLh,
            Cmw + xcg*b/S*CLw <= Vh*CLh,
            Vh <= Sh*tb["L"]/S**2*b
            ]

        lc = LinkedConstraintSet([tb, constraints], exclude=["W"])

        Model.__init__(self, None, lc, **kwargs)

class Avionics(Model):
    def __init__(self, **kwargs):
        W_fc = Variable("W_{fc}", 4, "lbf", "Flight controller weight")
        W_batt = Variable("W_{batt}", 4, "lbf", "Battery weight")
        W = Variable("W", "lbf", "Avionics Weight")
        m_fac = Variable("m_{fac}", 1.0, "-", "Avionics weight margin factor")

        constraints = [W/m_fac >= W_fc + W_batt]

        Model.__init__(self, None, constraints, **kwargs)

class Weight(ConstraintSet):
    """
    Weight brakedown of aircraft
    """
    def __init__(self, center_loads, zf_loads, **kwargs):

        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "Total fuel weight")
        W_fueltank = Variable('W_{fuel-tank}', 4, 'lbf', 'Fuel tank weight')
        W_skid = Variable("W_{skid}", 4, "lbf", "Skid weight")

        # gobal vars

        W_pay = Variable("W_{pay}", 10, "lbf", "Payload weight")
        W_zfw = Variable("W_{zfw}", "lbf", "Zero fuel weight")

        constraints = [
            SummingConstraintSet(W_cent, "W", center_loads,
                                 [W_fueltot, W_pay, W_skid, W_fueltank]),
            SummingConstraintSet(W_zfw, "W", zf_loads,
                                 [W_pay, W_skid, W_fueltank])
            ]

        ConstraintSet.__init__(self, constraints, **kwargs)

class SystemRequirements(Model):
    def __init__(self, **kwargs):

        mtow = Variable("MTOW", "lbf", "Max take off weight")
        t_loiter = Variable("t_{loiter}", "days", "time loitering")

        s = Variable("s", "-", "system margin factor")
        mtowreq = Variable("MTOW_{req}", 150, "lbf", "Max take off weight")
        t_loiterreq = Variable("t_{loiter-req}", 5, "days", "time loitering")

        constraints = [mtow <= s*mtowreq,
                       t_loiter*s >= t_loiterreq]

        Model.__init__(self, None, constraints, **kwargs)

class GasMALE(Model):
    """
    This model has a rubber engine and as many non fixed parameters as
    possible.  Model should be combed for variables that are incorrectly
    fixed.
    """
    def __init__(self, h_station=15000, wind=False, DF70=False, Nclimb=10,
                 Nloiter=5, margin=False, **kwargs):

        mission = Mission(h_station, wind, DF70, Nclimb, Nloiter)
        engineweight = EngineWeight(DF70)
        wing = Wing()
        tail = Tail(mission.submodels[0])
        avionics = Avionics()
        fuselage = Fuselage()
        center_loads = [fuselage, avionics, engineweight]
        zf_loads = center_loads + [tail, wing]
        weight = Weight(center_loads, zf_loads)

        self.submodels = zf_loads + [weight, mission]

        constraints = []

        if margin:
            sq = SystemRequirements()
            self.submodels.append(sq)

        lc = LinkedConstraintSet([self.submodels, constraints],
                                 include_only=INCLUDE)


        objective = 1/mission["t_{loiter}"]

        Model.__init__(self, objective, lc, **kwargs)

if __name__ == "__main__":
    M = GasMALE(DF70=True)
    M.substitutions.update({"t_{loiter}": 6})
    M.cost = M["MTOW"]
    if SIGNOMIALS:
        sol = M.localsolve("mosek")
    else:
        sol = M.solve("mosek")

    # M.substitutions.update({"t_{loiter}": 5.8})
    # M.substitutions.update({"P_{pay}": 100})
    # if SIGNOMIALS:
    #     sol = M.localsolve("mosek")
    # else:
    #     sol = M.solve("mosek")

