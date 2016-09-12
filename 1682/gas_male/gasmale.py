""" gas_male_rubber.py """
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit import LinkedConstraintSet
#from gpkit.tools import BoundedConstraintSet
from gpkit.tools import te_exp_minus1
PLOT = False

class Mission(Model):
    def __init__(self, h_station, wind, DF70, discrit, **kwargs):

        if discrit:
            self.submodels = [
                Climb(5, [0.502]*5, np.linspace(1, 5000, 5), False, wind,
                      5000, DF70),
                Cruise(5, [0.684]*5, [5000]*5, False, wind, 180, DF70),
                Climb(5, [0.567]*5, np.linspace(5000, h_station, 5), False,
                      wind, 10000, DF70),
                Loiter(20, [0.647]*20, [h_station]*20, True, wind, DF70),
                Cruise(5, [0.684]*5, [5000]*5, False, wind, 200, DF70)
                ]
        else:
            self.submodels = [
                Climb(1, [0.502], [5000], False, wind, 5000, DF70),
                Cruise(1, [0.684], [5000], False, wind, 180, DF70),
                Climb(1, [0.567], [h_station], False, wind, 10000, DF70),
                Loiter(5, [0.647]*5, [h_station]*5, True, wind, DF70),
                Cruise(1, [0.684], [5000], False, wind, 200, DF70)
                ]

        mtow = Variable("MTOW", "lbf", "max take off weight")
        W_zfw = Variable("W_{zfw}", "lbf", "zero fuel weight")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "total fuel weight")

        constraints = [
            mtow >= self.submodels[0]["W_{start}"],
            W_zfw <= self.submodels[-1]["W_{finish}"],
            W_fueltot >= sum(fs["W_{fuel-fs}"] for fs in self.submodels)
            ]

        for i, fs in enumerate(self.submodels[1:]):
            constraints.extend([
                self.submodels[i]["W_{finish}"] == fs["W_{start}"]
                ])

        lc = LinkedConstraintSet(
            [fs for fs in self.submodels, constraints], exclude={
                "t_{loiter}", "\\delta_h", "h_{dot}", "h", "T_{atm}",
                "\\mu", "\\rho", "W_{start}", "W_{end}", "W_{fuel}",
                "W_{fuel-fs}", "W_{finish}", "W_{begin}", "C_L", "C_D",
                "V", "\\eta_{prop}", "P_{shaft}", "T", "BSFC", "RPM",
                "P_{shaft-max}", "L_factor", "P_{shaft-tot}", "t", "R",
                "Re", "C_{f-fuse}", "C_{D-fuse}", "Re_{fuse}", "c_{dp}",
                "V_{wind}", "z_{bre}", "h_{loss}", "K_{fuse}"})

        Model.__init__(self, None, lc, **kwargs)

class FlightSegment(Model):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, **kwargs):

        self.aero = Aerodynamics(N)
        self.fuel = Fuel(N)
        self.slf = SteadyLevelFlight(N, eta_p)
        self.engine = Engine(N, alt, onStation, DF70)
        self.atm = Atmosphere(N, alt)
        self.wind = Wind(N, alt, wind)

        self.submodels = [self.aero, self.fuel, self.slf, self.engine,
                          self.atm, self.wind]

class Cruise(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, R, DF70, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetrange = BreguetRange(N, R)

        self.submodels.extend([breguetrange])

        lc = LinkedConstraintSet([self.submodels])

        Model.__init__(self, None, lc, **kwargs)

class Loiter(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        t_loiter = Variable("t_{loiter}", "days", "time loitering")

        constraints = [breguetendurance["t"] >= t_loiter/N]

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, constraints])

        Model.__init__(self, None, lc, **kwargs)

class Climb(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, dh, DF70, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        deltah = Variable("\\delta_h", dh, "ft", "altitude difference")
        h_dot = VectorVariable(N, "h_{dot}", "ft/min", "climb rate")
        h_dotmin = Variable("h_{dot-min}", 100, "ft/min",
                            "minimum climb rate")

        constraints = [
            h_dot*breguetendurance["t"] >= deltah/N,
            h_dot >= h_dotmin,
            self.slf["T"] >= (0.5*self.slf["\\rho"]*self.slf["V"]**2*
                              self.slf["C_D"]*self.slf["S"] +
                              self.slf["W_{begin}"]*h_dot/self.slf["V"])
            ]

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, constraints])

        Model.__init__(self, None, lc, **kwargs)

class Atmosphere(Model):
    """
    Model to capture density changes with altitude
    """
    def __init__(self, N, alt, **kwargs):

        h = VectorVariable(N, "h", alt, "ft", "altitude")
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
        h_ref = Variable("h_{ref}", 15000, "ft", "ref altitude")
        rho = VectorVariable(N, "\\rho", "kg/m^3", "air density")

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
        W_end = VectorVariable(N, "W_{end}", "lbf", "segment-end weight")
        W_fuel = VectorVariable(N, "W_{fuel}", "lbf",
                                "segment-fuel weight")
        W_fuelfs = Variable("W_{fuel-fs}", "lbf",
                            "flight segment fuel weight")
        W_finish = Variable("W_{finish}", "lbf",
                            "weight at beginning of flight segment")
        W_begin = VectorVariable(N, "W_{begin}", "lbf", "segment-begin weight")

        # end of first segment weight + first segment fuel weight must be
        # greater  than MTOW.  Each end of segment weight must be greater
        # than the next end of segment weight + the next segment fuel weight.
        # The last end segment weight must be greater than the zero fuel
        # weight
        constraints = [W_start >= W_end[0] + W_fuel[0],
                       W_fuelfs >= W_fuel.sum(),
                       W_end[-1] >= W_finish,
                       W_begin[0] >= W_start]

        if N == 1:
            pass
        else:
            constraints.extend([W_end[:-1] >= W_end[1:] + W_fuel[1:],
                                W_begin[1:] == W_end[:-1],
                               ])

        Model.__init__(self, None, constraints, **kwargs)

class SteadyLevelFlight(Model):
    """
    Captures steady level flight mode
    """
    def __init__(self, N, eta_p, **kwargs):

        CD = VectorVariable(N, "C_D", "-", "Drag coefficient")
        CL = VectorVariable(N, "C_L", "-", "Lift coefficient")
        V = VectorVariable(N, "V", "m/s", "cruise speed")
        S = Variable("S", "ft^2", "wing area")
        eta_prop = VectorVariable(N, "\\eta_{prop}", eta_p, "-",
                                  "propulsive efficiency")
        P_shaft = VectorVariable(N, "P_{shaft}", "hp", "Shaft power")
        T = VectorVariable(N, "T", "lbf", "Thrust")

        rho = VectorVariable(N, "\\rho", "kg/m^3", "air density")
        W_end = VectorVariable(N, "W_{end}", "lbf", "segment-end weight")
        W_begin = VectorVariable(N, "W_{begin}", "lbf", "segment-begin weight")

        # Climb model
        # Currently climb rate is being optimized to reduce fuel consumption.
        # In future, could implement min climb rate.

        constraints = [
            P_shaft == T*V/eta_prop,
            T >= 0.5*rho*V**2*CD*S,
            0.5*rho*CL*S*V**2 == (W_end*W_begin)**0.5,
            ]
        # Propulsive efficiency variation with different flight segments,
        # will change depending on propeller characteristics

        Model.__init__(self, None, constraints, **kwargs)

class Engine(Model):
    """
    Engine performance and weight model for small engine
    """
    def __init__(self, N, alt, onStation, DF70, **kwargs):

        h = VectorVariable(N, "h", alt, "ft", "altitude")
        h_ref = Variable("h_{ref}", 15000, "ft", "ref altitude")
        P_shaft = VectorVariable(N, "P_{shaft}", "hp", "Shaft power")
        bsfc = VectorVariable(N, "BSFC", "lb/hr/hp",
                              "brake specific fuel consumption")
        rpm = VectorVariable(N, "RPM", "rpm", "Engine operating RPM")
        P_avn = Variable("P_{avn}", 40, "watts", "avionics power")
        P_pay = Variable("P_{pay}", 10, "watts", "payload power")
        P_shafttot = VectorVariable(N, "P_{shaft-tot}", "hp",
                                    "total power, avionics included")
        eta_alternator = Variable("\\eta_{alternator}", 0.8, "-",
                                  "alternator efficiency")
        lfac = [1 - 0.906**(1/0.15)*(v.value/h_ref.value)**0.92
                for v in h]
        h_loss = VectorVariable(N, "h_{loss}", lfac, "-",
                                "Max shaft power loss factor")
        P_shaftmax = VectorVariable(N, "P_{shaft-max}",
                                    "hp", "Max shaft power at altitude")

        if DF70:
            W_engtot = Variable("W_{eng-tot}", 7.1, "lbf",
                                "Installed engine weight")
            P_shaftmaxmsl = Variable("P_{shaft-maxMSL}", 5.17, "hp",
                                     "Max shaft power at MSL")
            rpm_max = Variable("RPM_{max}", 7698, "rpm", "Maximum RPM")
            bsfc_min = Variable("BSFC_{min}", 0.3162, "kg/kW/hr",
                                "Minimum BSFC")

            constraints = [
                (bsfc/bsfc_min)**35.7 >= (2.29*(rpm/rpm_max)**8.02 +
                                          0.00114*(rpm/rpm_max)**-38.3),
                (P_shafttot/P_shaftmax)**0.1 == 0.999*(rpm/rpm_max)**0.294,
                ]
        else:
            P_shaftref = Variable("P_{shaft-ref}", 2.295, "hp",
                                  "reference shaft power")
            W_engref = Variable("W_{eng-ref}", 4.4107, "lbf",
                                "Reference engine weight")
            W_eng = Variable("W_{eng}", "lbf", "engine weight")
            W_engtot = Variable("W_{eng-tot}", "lbf",
                                "Installed engine weight")
            bsfc_min = Variable("BSFC_{min}", 0.32, "kg/kW/hr",
                                "Minimum BSFC")
            P_shaftmaxmsl = Variable("P_{shaft-maxMSL}", "hp",
                                     "Max shaft power at MSL")
            rpm_max = Variable("RPM_{max}", 9000, "rpm", "Maximum RPM")

            constraints = [
                W_eng/W_engref >= 0.5538*(P_shaftmaxmsl/P_shaftref)**1.075,
                W_engtot >= 2.572*W_eng**0.922*units("lbf")**0.078,
                (bsfc/bsfc_min)**0.129 >= (2*.486*(rpm/rpm_max)**-0.141 +
                                           0.0268*(rpm/rpm_max)**9.62),
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
        z_bre = VectorVariable(N, "z_{bre}", "-", "breguet coefficient")
        t = VectorVariable(N, "t", "days", "time per flight segment")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        P_shafttot = VectorVariable(N, "P_{shaft-tot}", "hp",
                                    "total power, avionics included")
        bsfc = VectorVariable(N, "BSFC", "lb/hr/hp",
                              "brake specific fuel consumption")
        W_end = VectorVariable(N, "W_{end}", "lbf", "segment-end weight")
        W_fuel = VectorVariable(N, "W_{fuel}", "lbf",
                                "segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")

        constraints = [
            z_bre >= P_shafttot*t*bsfc*g/W_end,
            f_fueloil*W_fuel/W_end >= te_exp_minus1(z_bre, 3)
            ]

        Model.__init__(self, None, constraints, **kwargs)

class BreguetRange(Model):
    """
    Discritized Breguet Range model
    """
    def __init__(self, N, R, **kwargs):
        z_bre = VectorVariable(N, "z_{bre}", "-", "breguet coefficient")
        t = VectorVariable(N, "t", "days", "time per flight segment")
        R = Variable("R", R, "nautical_miles", "range to station")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        P_shafttot = VectorVariable(N, "P_{shaft-tot}", "hp",
                                    "total power, avionics included")
        bsfc = VectorVariable(N, "BSFC", "lb/hr/hp",
                              "brake specific fuel consumption")
        W_end = VectorVariable(N, "W_{end}", "lbf", "segment-end weight")
        V = VectorVariable(N, "V", "m/s", "cruise speed")
        W_fuel = VectorVariable(N, "W_{fuel}", "lbf",
                                "segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")

        constraints = [
            z_bre >= P_shafttot*t*bsfc*g/W_end,
            R/N <= V*t,
            f_fueloil*W_fuel/W_end >= te_exp_minus1(z_bre, 3)
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
        l_fuse = Variable("l_{fuse}", "ft", "fuselage length")
        Refuse = VectorVariable(N, "Re_{fuse}", "-",
                                "fuselage Reynolds number")
        Re_ref = Variable("Re_{ref}", 3e5, "-", "Reference Re for cdp")
        cdp = VectorVariable(N, "c_{dp}", "-", "wing profile drag coeff")

        CD = VectorVariable(N, "C_D", "-", "Drag coefficient")
        CL = VectorVariable(N, "C_L", "-", "Lift coefficient")
        V = VectorVariable(N, "V", "m/s", "cruise speed")
        S = Variable("S", "ft^2", "wing area")
        rho = VectorVariable(N, "\\rho", "kg/m^3", "air density")
        mu_atm = VectorVariable(N, "\\mu", "N*s/m^2", "Dynamic viscosity")

        constraints = [
            CD >= CDfuse*2 + cdp*1.3 + CL**2/(pi*e*AR),
            #jh01
            cdp >= ((0.0075 + 0.002*CL**2 + 0.00033*CL**10)*(Re/Re_ref)**-0.4),
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

class Weight(Model):
    """
    Weight brakedown of aircraft
    """
    def __init__(self, DF70, **kwargs):

        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        W_fuse = Variable("W_{fuse}", "lbf", "fuselage weight")
        W_wing = Variable("W_{wing}", "lbf", "Total wing structural weight")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "total fuel weight")
        W_fueltank = Variable('W_{fuel-tank}', 4, 'lbf', 'fuel tank weight')
        W_skid = Variable("W_{skid}", 4, "lbf", "skid weight")
        m_tail = Variable("m_{tail}", 1.587, "kg", "tail mass")
        m_rib = Variable("m_{rib}", 1.36, "kg", "rib mass")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")

        # gobal vars
        m_fuse = Variable("m_{fuse}", "kg", "fuselage mass")
        m_cap = Variable("m_{cap}", "kg", "Cap mass")
        m_skin = Variable("m_{skin}", "kg", "Skin mass")

        W_pay = Variable("W_{pay}", 10, "lbf", "Payload weight")
        W_avionics = Variable("W_{avionics}", 8, "lbf", "Avionics weight")
        W_zfw = Variable("W_{zfw}", "lbf", "Zero fuel weight")

        if DF70:
            W_engtot = Variable("W_{eng-tot}", 7.1, "lbf",
                                "Installed engine weight")
        else:
            W_engtot = Variable("W_{eng-tot}", "lbf",
                                "Installed engine weight")

        constraints = [
            W_wing >= m_skin*g + 1.2*m_cap*g,
            W_fuse >= m_fuse*g + m_rib*g,
            W_cent >= (W_fueltot + W_pay + W_engtot + W_fuse + W_avionics +
                       W_skid + W_fueltank),
            W_zfw >= (W_pay + W_engtot + W_fuse + W_wing + m_tail*g +
                      W_avionics + W_skid + W_fueltank)
            ]

        Model.__init__(self, None, constraints, **kwargs)

class Structures(Model):
    """
    Structural wing model.  Simple beam.
    """
    def __init__(self, **kwargs):

        # Structural parameters
        rho_skin = Variable("\\rho_{skin}", 0.1, "g/cm^2",
                            "Wing Skin Density")
        rho_cap = Variable("\\rho_{cap}", 1.76, "g/cm^3", "Density of CF cap")
        E_cap = Variable("E_{cap}", 2e7, "psi", "Youngs modulus of CF cap")
        sigma_cap = Variable("\\sigma_{cap}", 475e6, "Pa", "Cap stress")

        # Structural lengths
        h_spar = Variable("h_{spar}", "m", "Spar height")
        t_cap = Variable("t_{cap}", 0.028, "in", "Spar cap thickness")
        #arbitrarily placed based on available cf
        w_cap = Variable("w_{cap}", "in", "Spar cap width")
        c = Variable("c", "ft", "Wing chord")
        #assumes straight, untapered wing

        # Structural ratios
        tau = Variable("\\tau", 0.115, "-", "Airfoil thickness ratio")
        #find better number
        LoverA = Variable("LoverA", "lbf/ft^2", "Wing loading")

        # Structural areas
        A_capcent = Variable("A_{capcent}", "m**2", "Cap area at center")
        #currently assumes constant area

        # Structural volumes
        Vol_cap = Variable("Vol_{cap}", "m**3", "Cap volume")

        # Structural evaluation parameters
        M_cent = Variable("M_cent", "N*m", "Center bending moment")
        F = Variable("F", "N", "Load on wings")
        N_max = Variable("N_{max}", 5, "-", "Load factor")
        #load rating for max number of g"s
        P_cap = Variable("P_{cap}", "N", "Cap load")
        delta_tip = Variable("\\delta_{tip}", "ft", "Tip deflection")
        delta_tip_max = Variable("\\delta_{tip-max}", 0.2, "-",
                                 "max tip deflection ratio")

        S = Variable("S", "ft^2", "wing area")
        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        b = Variable("b", "ft", "Span")
        m_cap = Variable("m_{cap}", "kg", "Cap mass")
        mtow = Variable("MTOW", "lbf", "max take off weight")

        constraints = [m_skin >= rho_skin*S*2,
                       F >= W_cent*N_max,
                       c == S/b,
                       M_cent >= b*F/8,
                       P_cap >= M_cent/h_spar,
                       A_capcent >= P_cap/sigma_cap,
                       Vol_cap >= A_capcent*b/3,
                       m_cap == rho_cap*Vol_cap,
                       h_spar <= tau*c,
                       w_cap == A_capcent/t_cap,
                       LoverA == mtow/S,
                       delta_tip == b**2*sigma_cap/(4*E_cap*h_spar),
                       delta_tip/b <= delta_tip_max]

        Model.__init__(self, None, constraints, **kwargs)

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
        Vol_fuse = Variable("Vol_{fuse}", "m**3", "fuselage volume")

        m_fuse = Variable("m_{fuse}", "kg", "fuselage mass")
        rho_skin = Variable("\\rho_{skin}", 0.1, "g/cm^2",
                            "Wing Skin Density")
        S_fuse = Variable("S_{fuse}", "ft^2", "Fuselage surface area")
        l_fuse = Variable("l_{fuse}", "ft", "fuselage length")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "total fuel weight")
        l_cent = Variable("l_{cent}", "ft", "center fuselage length")
        Vol_avionics = Variable("Vol_{avionics}", 0.125, "ft^3",
                                "Avionics volume")
        Vol_pay = Variable("Vol_{pay}", 1, "ft^3", "Payload volume")

        constraints = [m_fuse >= S_fuse*rho_skin,
                       l_cent == fr*w_cent,
                       l_fuse >= l_cent*1.1,
                       (l_fuse/k1fuse)**3 == Vol_fuse,
                       (S_fuse/k2fuse)**3 == Vol_fuse**2,
                       Vol_fuse >= l_cent*w_cent**2,
                       Vol_fuel >= W_fueltot/rho_fuel,
                       l_cent*w_cent**2 >= Vol_fuel+Vol_avionics+Vol_pay
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Wind(Model):
    """
    Model for wind speed
    wind = True, wind speed has specific value
    """
    def __init__(self, N, alt, wind, **kwargs):

        V = VectorVariable(N, "V", "m/s", "cruise speed")
        h = VectorVariable(N, "h", alt, "ft", "altitude")

        if wind:

            V_wind = Variable("V_{wind}", 25, "m/s", "wind speed")
            constraints = [V >= V_wind]

        else:

            V_wind = VectorVariable(N, "V_{wind}", "m/s", "wind speed")
            h_ref = Variable("h_{ref}", 15000, "ft", "ref altitude")
            V_ref = Variable("V_{ref}", 25, "m/s", "wind speed")

            constraints = [(V_wind/V_ref) >= 0.6462*(h/h_ref) + 0.3538,
                           V >= V_wind,
                          ]

        Model.__init__(self, None, constraints, **kwargs)

class GasMALE(Model):
    """
    This model has a rubber engine and as many non fixed parameters as
    possible.  Model should be combed for variables that are incorrectly
    fixed.
    """
    def __init__(self, h_station=15000, wind=False, DF70=False,
                 discrit=False, **kwargs):

        mission = Mission(h_station, wind, DF70, discrit)
        weight = Weight(DF70)
        fuselage = Fuselage()
        structures = Structures()

        self.submodels = [mission, weight, fuselage, structures]

        constraints = []

        lc = LinkedConstraintSet([self.submodels, constraints], exclude={
            "t_{loiter}", "\\delta_h", "h_{dot}", "h", "T_{atm}",
            "\\mu", "\\rho", "W_{start}", "W_{end}", "W_{fuel}",
            "W_{fuel-fs}", "W_{finish}", "W_{begin}", "C_L", "C_D",
            "V", "\\eta_{prop}", "P_{shaft}", "T", "BSFC", "RPM",
            "P_{shaft-max}", "L_factor", "P_{shaft-tot}", "t", "R",
            "Re", "C_{f-fuse}", "C_{D-fuse}", "Re_{fuse}", "c_{dp}",
            "V_{wind}", "z_{bre}", "h_{loss}", "K_{fuse}"})

        objective = 1/mission["t_{loiter}"]

        Model.__init__(self, objective, lc, **kwargs)

if __name__ == "__main__":
    M = GasMALE()
    sol = M.solve("mosek")
