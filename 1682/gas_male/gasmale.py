""" gas_male_rubber.py """
from numpy import pi
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit import LinkedConstraintSet, ConstraintSet
from helpers import SummingConstraintSet
#from gpkit.tools import BoundedConstraintSet
from gpkit.tools import te_exp_minus1
PLOT = False

INCLUDE = ["l_{fuse}", "MTOW", "t_{loiter}", "S", "b", "AR", "P_{shaft-maxMSL}",
           "S_{fuse}", "W_{cent}", "W_{zfw}", "W_{fuel-tot}", "g"]

class Mission(Model):
    def __init__(self, h_station, wind, DF70, N, Nloiter, **kwargs):

        self.submodels = [
            Climb(N, [0.502]*N, [5000]*N, False, wind, DF70, dh=5000),
            Cruise(N, [0.684]*N, [5000]*N, False, wind, DF70, R=180),
            Climb(N, [0.567]*N, [h_station]*N, False, wind, DF70, dh=10000),
            Loiter(Nloiter, [0.647]*Nloiter, [h_station]*Nloiter, True, wind,
                   DF70),
            Cruise(N, [0.684]*N, [5000]*N, False, wind, DF70, R=200)
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

        self.submodels = [self.aero, self.fuel, self.slf, self.engine,
                          self.atm, self.wind]

class Cruise(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, R=200, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetrange = BreguetRange(N, R)

        self.submodels.extend([breguetrange])

        self.constraints = []

        lc = LinkedConstraintSet([self.submodels], exclude=self.exclude)

        Model.__init__(self, None, lc, **kwargs)

class Loiter(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70)

        breguetendurance = BreguetEndurance(N)

        t_loiter = Variable("t_{loiter}", "days", "time loitering")

        constraints = [breguetendurance["t"] >= t_loiter/N]

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, constraints],
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
        constraints = [
            h_dot*breguetendurance["t"] >= deltah/N,
            h_dot >= h_dotmin,
            self.slf["T"] >= (0.5*self.slf["\\rho"]*self.slf["V"]**2*
                              self.slf["C_D"]*self.slf["S"] +
                              self.slf["W_{N}"]*h_dot/self.slf["V"])
            ]

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, constraints],
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
                       W_n[0] >= W_start]

        if N == 1:
            pass
        else:
            constraints.extend([W_nplus1[:-1] >= W_nplus1[1:] + W_fuel[1:],
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

        constraints = [
            z_bre >= P_shafttot*t*bsfc*g/W_nplus1,
            f_fueloil*W_fuel/W_nplus1 >= te_exp_minus1(z_bre, 3)
            ]

        Model.__init__(self, None, constraints, **kwargs)

class BreguetRange(Model):
    """
    Discritized Breguet Range model
    """
    def __init__(self, N, R, **kwargs):
        z_bre = VectorVariable(N, "z_{bre}", "-", "Breguet coefficient")
        t = VectorVariable(N, "t", "days", "Time per flight segment")
        R = Variable("R", R, "nautical_miles", "Range to station")
        f_fueloil = Variable("f_{(fuel/oil)}", 0.98, "-", "Fuel-oil fraction")
        P_shafttot = VectorVariable(N, "P_{shaft-tot}", "hp",
                                    "Total power, avionics included")
        bsfc = VectorVariable(N, "BSFC", "lb/hr/hp",
                              "Brake specific fuel consumption")
        W_nplus1 = VectorVariable(N, "W_{N+1}", "lbf", "vector-end weight")
        V = VectorVariable(N, "V", "m/s", "Cruise speed")
        W_fuel = VectorVariable(N, "W_{fuel}", "lbf",
                                "Segment-fuel weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")

        constraints = [
            z_bre >= P_shafttot*t*bsfc*g/W_nplus1,
            R/N <= V*t,
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
        m_fac = Variable("m_{fac}", 1.0, "-", "c_{dp} margin factor")

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


class Wing(Model):
    """
    Structural wing model.  Simple beam.
    """
    def __init__(self, **kwargs):

        # Structural parameters
        rho_skin = Variable("\\rho_{skin}", 0.1, "g/cm^2",
                            "Wing skin density")
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
                                 "Max tip deflection ratio")

        S = Variable("S", "ft^2", "wing area")
        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        b = Variable("b", "ft", "Span")
        m_cap = Variable("m_{cap}", "kg", "Cap mass")
        mtow = Variable("MTOW", "lbf", "Max take off weight")
        m_cap = Variable("m_{cap}", "kg", "Cap mass")
        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        W = Variable("W", "lbf", "Total wing structural weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        m_fac = Variable("m_{fac}", 1.0, "-", "Wing weight margin factor")

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
                       delta_tip/b <= delta_tip_max,
                       W/m_fac >= m_skin*g + 1.2*m_cap*g
                      ]

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

class Tail(Model):
    def __init__(self, **kwargs):
        W_vtail = Variable("W_{v-tail}", 3.4999, "lbf", "V-Tail weight")
        m_fac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "Tail weight")

        constraints = [W/m_fac >= W_vtail]

        Model.__init__(self, None, constraints, **kwargs)

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
    def __init__(self, h_station=15000, wind=False, DF70=False, N=1,
                 Nloiter=5, margin=False, **kwargs):

        mission = Mission(h_station, wind, DF70, N, Nloiter)
        engineweight = EngineWeight(DF70)
        wing = Wing()
        tail = Tail()
        avionics = Avionics()
        fuselage = Fuselage()
        center_loads = [wing, fuselage, avionics, engineweight]
        zf_loads = center_loads + [tail]
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
    sol = M.solve("mosek")

