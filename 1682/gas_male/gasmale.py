""" gas_male_rubber.py """
from numpy import pi
from operator import mul
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit import LinkedConstraintSet, ConstraintSet, VarKey
from gpkit import SignomialsEnabled
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
from helpers import SummingConstraintSet
#from gpkit.tools import BoundedConstraintSet
plt.rcParams.update({'font.size':19})
PLOT = False
SIGNOMIALS = False

INCLUDE = ["l_{fuse}", "MTOW", "t_{loiter}", "S", "b", "AR", "W_{zfw}",
           "P_{shaft-maxMSL}", "S_{fuse}", "W_{cent}", "W_{fuel-tot}", "g",
           "V_{stall}"]

class Mission(Model):
    def __init__(self, h_station, wind, DF70, Nclimb, Nloiter, dragcomps,
                 **kwargs):

        self.submodels = [
            TakeOff(1, [0.684], [0.1], False, wind, DF70, dragcomps),
            Climb(Nclimb, [0.502]*Nclimb, np.linspace(0, h_station,
                                                      Nclimb+1)[1:],
                  False, wind, DF70, dragcomps, dh=h_station),
            Cruise(1, [0.684], [h_station], False, wind, DF70, dragcomps,
                   R=180),
            Loiter(Nloiter, [0.647]*Nloiter, [h_station]*Nloiter, True, wind,
                   DF70, dragcomps),
            Cruise(1, [0.684], [h_station], False, wind, DF70, dragcomps, R=200)
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

        cs = ConstraintSet([fs for fs in self.submodels, constraints])

        linked = {}
        for c in dragcomps:
            for name in ["l_{ref}", "S_{ref}"]:
                vks = cs.varkeys[name]
                for vk in vks:
                    descr = dict(vk.descr)
                    if c.name in descr["models"]:
                        descr.pop("value", None)
                        descr["models"] = [c.name]
                        descr["modelnums"] = [c.num]
                        newvk = VarKey(**descr)
                        linked.update({vk: newvk})
        cs.subinplace(linked)

        lc = LinkedConstraintSet(cs, include_only=INCLUDE)


        Model.__init__(self, None, lc, **kwargs)

class FlightSegment(Model):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, dragcomps,
                 **kwargs):

        self.aero = Aerodynamics(N, dragcomps)
        self.fuel = Fuel(N)
        self.slf = SteadyLevelFlight(N, eta_p)
        self.engine = EnginePerformance(N, alt, onStation, DF70)
        self.atm = Atmosphere(N, alt)
        self.wind = Wind(N, alt, wind)

        self.N = N
        self.include = ["V", "\\rho", "\\mu", "BSFC", "W_{N}", "W_{N+1}",
                        "P_{shaft}", "P_{shaft-tot}", "C_D", "C_L", "S",
                        "W_{fuel}", "h"]

        self.constraints = []

        self.submodels = [self.aero, self.fuel, self.slf, self.engine,
                          self.atm, self.wind]

class TakeOff(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, dragcomps,
                 **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70,
                               dragcomps)

        breguetendurance = BreguetEndurance(N)

        self.submodels.extend([breguetendurance])

        t_to = Variable("t_{TO}", 10, "minutes", "take off time")
        Vstall = Variable("V_{stall}", "m/s", "stall speed")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")

        self.constraints = [
            breguetendurance["t"] >= t_to,
            self.slf["V"] >= 1.3*Vstall,
            Vstall >= (self.fuel["W_{N}"][0]*2/rhosl/self.aero["S"]/1.5)**0.5
            ]

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 include_only=self.include)

        Model.__init__(self, None, lc, **kwargs)

class Cruise(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, dragcomps, R=200,
                 **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70,
                               dragcomps)

        breguetendurance = BreguetEndurance(N)

        R = Variable("R", R, "nautical_miles", "Range to station")

        self.submodels.extend([breguetendurance])

        self.constraints.extend([R/N <= self.slf["V"]*breguetendurance["t"]])

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 include_only=self.include)

        Model.__init__(self, None, lc, **kwargs)

class Loiter(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, dragcomps,
                 **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70,
                               dragcomps)

        breguetendurance = BreguetEndurance(N)

        t_loiter = Variable("t_{loiter}", "days", "time loitering")

        self.constraints.extend([breguetendurance["t"] >= t_loiter/N])

        self.submodels.extend([breguetendurance])

        lc = LinkedConstraintSet([self.submodels, self.constraints],
                                 include_only=self.include)

        Model.__init__(self, None, lc, **kwargs)

class Climb(FlightSegment):
    def __init__(self, N, eta_p, alt, onStation, wind, DF70, dragcomps,
                 dh=15000,
                 **kwargs):
        FlightSegment.__init__(self, N, eta_p, alt, onStation, wind, DF70,
                               dragcomps)

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
                                 include_only=self.include)

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
            # TCS([z_bre >= P_shafttot*t*bsfc*g/(W_nplus1*W_n)**0.5]),
            # TCS([z_bre >= P_shafttot*t*bsfc*g/W_nplus1]),
            f_fueloil*W_fuel/W_nplus1 >= te_exp_minus1(z_bre, 3)
            ]

        Model.__init__(self, None, constraints, **kwargs)

class ComponentDrag(Model):
    def __init__(self, N, comp, **kwargs):

        CDA = VectorVariable(N, "CDA", "-",
                             "%s area drag normalized by wing area" % comp.name)
        Cf = VectorVariable(N, "C_f", "-",
                            "%s skin friction coefficient" % comp.name)
        Re = VectorVariable(N, "Re", "-", "%s reynolds number" % comp.name)
        S = Variable("S", "ft^2", "wing area")
        rho = VectorVariable(N, "\\rho", "kg/m^3", "Air density")
        mu_atm = VectorVariable(N, "\\mu", "N*s/m^2", "Dynamic viscosity")
        V = VectorVariable(N, "V", "m/s", "Cruise speed")

        constraints = [CDA >= Cf*comp["S_{ref}"]/S,
                       Re == V*rho*comp["l_{ref}"]/mu_atm,
                       Cf >= 0.455/Re**0.3
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Aerodynamics(Model):
    """
    Aero model assuming jh01 airfoil, designed by Mark Drela
    """
    def __init__(self, N, dragcomps=None, jh01=True, **kwargs):

        CLmax = Variable("C_{L-max}", 1.39, "-", "Maximum lift coefficient")
        e = Variable("e", 0.9, "-", "Spanwise efficiency")
        AR = Variable("AR", "-", "Aspect ratio")
        b = Variable("b", "ft", "Span")
        Re = VectorVariable(N, "Re", "-", "Reynolds number")

        Re_ref = Variable("Re_{ref}", 3e5, "-", "Reference Re for cdp")
        cdp = VectorVariable(N, "c_{dp}", "-", "wing profile drag coeff")

        CD = VectorVariable(N, "C_D", "-", "Drag coefficient")
        CL = VectorVariable(N, "C_L", "-", "Lift coefficient")
        V = VectorVariable(N, "V", "m/s", "Cruise speed")
        S = Variable("S", "ft^2", "wing area")
        rho = VectorVariable(N, "\\rho", "kg/m^3", "Air density")
        mu_atm = VectorVariable(N, "\\mu", "N*s/m^2", "Dynamic viscosity")
        m_fac = Variable("m_{fac}", 1.7, "-", "CDA0 margin factor")
        CDA0 = VectorVariable(N, "CDA_0", "-", "sum of component drag")

        constraints = [CD >= CDA0 + cdp + CL**2/(pi*e*AR),
                       b**2 == S*AR,
                       CL <= CLmax,
                       Re == rho*V/mu_atm*(S/AR)**0.5,
                      ]

        if jh01:
            constraints.extend([
                (cdp/(Re/Re_ref)**-0.4)**0.00544 >= (0.33*CL**-0.0809 +
                                                     0.645*CL**0.045 +
                                                     7.35e-5*CL**12)])
        else:
            #sd7032
            constraints.extend([cdp >= ((0.006 + 0.005*CL**2 +
                                         0.00012*CL**10)*(Re/Re_ref)**-0.3)])

        if dragcomps:
            dragbuild = [ComponentDrag(N, c) for c in dragcomps]
            constraints.extend([
                CDA0/m_fac >= sum(db["CDA"] for db in dragbuild)])

            includes = ["\\rho", "\\mu", "S", "V"]
            for c, db in zip(dragcomps, dragbuild):
                linked = {}
                for vk in db.varkeys:
                    descr = dict(vk.descr)
                    if "ComponentDrag" == descr["models"][0]:
                        if descr["name"] not in includes:
                            descr["models"][0] = c.name
                            descr["modelnums"][0] = c.num
                            newvk = VarKey(**descr)
                            linked.update({vk: newvk})
                db.subinplace(linked)

            cs = LinkedConstraintSet([constraints, dragbuild],
                                     include_only=includes)

        else:
            cd0 = Variable("c_{d_0}", 0.002, "-", "non wing drag")
            constraints.extend([CDA0/m_fac >= cd0])
            cs = ConstraintSet([constraints])

        Model.__init__(self, None, cs, **kwargs)

class Beam(Model):
    def __init__(self, N, q, **kwargs):

        self.q = q
        self.N = N

        qbar = VectorVariable(N, "\\bar{q}", self.q, "-", "normalized loading")
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

    # def process_solution(self, sol):
    #     load = sol("W_{cent}")/sol("b")*self.q
    #     dx = sol("b")/2/(self.N-1)
    #     S = [0]*self.N
    #     for i in range(1, self.N):
    #         S[self.N-i-1] = S[self.N-i] + 0.5*dx*(load[self.N-i] +
    #                                               load[self.N-i-1])
    #     M = [0]*self.N
    #     for i in range(1, self.N):
    #         M[self.N-i-1] = M[self.N-i] + 0.5*dx*(S[self.N-i] + S[self.N-i-1])
    #     th = [0]*self.N
    #     for i in range(self.N-1):
    #         th[i+1] = (th[i] + 0.5*dx*(M[i] + M[i+1])/
    #                    sol("E_CapSpar, Wing, GasMALE")/sol("I")[i])
    #     d = [0]*self.N
    #     for i in range(self.N-1):
    #         d[i+1] = d[i] + 0.5*dx*(th[i] + th[i+1])
    #     load = load.to("N/m").magnitude
    #     for i in range(self.N-1):
    #         S[i] = S[i].to("N").magnitude
    #         M[i] = M[i].to("N*m").magnitude
    #         th[i+1] = th[i+1].to("dimensionless").magnitude
    #         d[i+1] = d[i+1].to("ft").magnitude

    #     fig, axis = plt.subplots(5)
    #     loading = [load, S, M, th, d]
    #     lunits = ["N/m", "N", "N*m", "-", "ft"]
    #     label = ["Loading", "Shear", "Moment", "Angle", "Deflection"]
    #     for ax, y, u, l in zip(axis, loading, lunits, label):
    #         ax.plot(dx.magnitude*np.linspace(0, 4, 5), y)
    #         ax.set_ylabel("%s [%s]" % (l, u), rotation=0)
    #         ax.set_yticks([min(y), max(y)])
    #         box = ax.get_position()
    #         ax.set_position([box.x0*1.5, box.y0, box.width*0.85, box.height])
    #         if not l == "Deflection":
    #             ax.set_xticklabels([])
    #             continue
    #         ax.set_xlabel("y [%s]" % u)
    #         ax.set_xticks(dx.magnitude*np.linspace(0, 4, 5))
    #     fig.savefig("sparforces.pdf", bbox_inches="tight")

    #     fig, axis = plt.subplots(3)
    #     d = len(sol("I"))
    #     dx = dx.magnitude*np.array([0, 1, 1, 2, 2, 3, 3, 4])
    #     I = np.insert(sol("I").magnitude, np.arange(d), sol("I").magnitude)
    #     w = np.insert(sol("w").magnitude, np.arange(d), sol("w").magnitude)
    #     t = np.insert(sol("\\vec{t_CapSpar, Wing, GasMALE}").magnitude,
    #                   np.arange(d),
    #                   sol("\\vec{t_CapSpar, Wing, GasMALE}").magnitude)
    #     dims = [I, w, t]
    #     lunits = ["m^4", "in", "in"]
    #     label = ["Inertia", "Width", "Thickness"]
    #     for ax, y, u, l in zip(axis, dims, lunits, label):
    #         ax.plot(dx, y)
    #         ax.set_ylabel("%s [%s]" % (l, u), rotation=0)
    #         ax.set_yticks([min(y), max(y)])
    #         box = ax.get_position()
    #         ax.set_position([box.x0*1.5, box.y0, box.width*0.85, box.height])
    #         if not l == "Thickness":
    #             ax.set_xticklabels([])
    #             continue
    #         ax.set_xlabel("y [ft]")
    #         ax.set_xticks(dx)
    #     fig.savefig("spardims.pdf", bbox_inches="tight")


def c_bar(lam, N):
    eta = np.linspace(0, 1, N)
    # c = np.array([1.0, 1.0, 0.83, 0.66, 0.5])
    c = 2/(1+lam)*(1+(lam-1)*eta)
    return c

class TubeSpar(Model):
    def __init__(self, N=5, **kwargs):

        # phyiscal properties
        rho_cfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        E = Variable("E", 2e7, "psi", "Youngs modulus of CF")
        sigma_cfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")

        # Structural lengths
        cb = c_bar(0.5, N)
        cbavg = (cb[:-1] + cb[1:])/2
        cbar = VectorVariable(N-1, "\\bar{c}", cbavg, "-",
                              "normalized chord at mid element")
        d = VectorVariable(N-1, "d", "in", "spar diameter")
        I = VectorVariable(N-1, "I", "m^4", "spar x moment of inertia")
        A = VectorVariable(N-1, "A", "in**2", "spar cross sectional area")
        dm = VectorVariable(N-1, "dm", "kg", "segment spar mass")
        m = Variable("m", "kg", "spar mass")

        S = Variable("S", "ft^2", "wing area")
        tau = Variable("\\tau", 0.115, "-", "Airfoil thickness ratio")
        b = Variable("b", "ft", "Span")

        N_max = Variable("N_{max}", 5, "-", "Load factor")
        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        kappa = Variable("\\kappa", 0.2, "-", "Max tip deflection ratio")

        beam = Beam(N, cb)
        self.submodels = [beam]

        constraints = [
            dm >= rho_cfrp*A*b/(N-1),
            m >= dm.sum(),
            S/b*cbar*tau >= d,
            4*I**2/A**2/(d/2)**2 + A/pi <= (d/2)**2,
            beam["\\bar{\\delta}"][-1] <= kappa,
            sigma_cfrp >= ((beam["\\bar{M}"][:-1] + beam["\\bar{M}"][1:])/
                           2*b*W_cent*N_max/4*(d/2)/I),
            beam["\\bar{EI}"] <= E*I/N_max/W_cent*b/(b/2)**3
            ]

        lc = LinkedConstraintSet([beam, constraints])

        Model.__init__(self, None, lc, **kwargs)

class CapSpar(Model):
    def __init__(self, N=5, **kwargs):

        # phyiscal properties
        rho_cfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        E = Variable("E", 2e7, "psi", "Youngs modulus of CFRP")
        sigma_cfrp = Variable("\\sigma_{CFRP}", 475e6, "Pa", "CFRP max stress")

        # Structural lengths
        cb = c_bar(0.5, N)
        cbavg = (cb[:-1] + cb[1:])/2
        cbar = VectorVariable(N-1, "\\bar{c}", cbavg, "-",
                              "normalized chord at mid element")
        t = VectorVariable(N-1, "t", "in", "spar cap thickness")
        hin = VectorVariable(N-1, "h_{in}", "in", "inner spar height")
        w = VectorVariable(N-1, "w", "in", "spar width")
        tshear = VectorVariable(N-1, "t_{shear}", "in",
                                "shear casing thickness")
        I = VectorVariable(N-1, "I", "m^4", "spar x moment of inertia")
        dm = VectorVariable(N-1, "dm", "kg", "segment spar mass")
        m = Variable("m", "kg", "spar mass")

        S = Variable("S", "ft^2", "wing area")
        tau = Variable("\\tau", 0.115, "-", "Airfoil thickness ratio")
        b = Variable("b", "ft", "Span")

        N_max = Variable("N_{max}", 5, "-", "Load factor")
        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")

        kappa = Variable("\\kappa", 0.2, "-", "Max tip deflection ratio")
        w_lim = Variable("w_{lim}", "-", "spar width to chord ratio")

        beam = Beam(N, cb)
        self.submodels = [beam]

        constraints = [
            dm >= rho_cfrp*w*t*b/(N-1) + rho_cfrp*b/2/(N-1)*2*tshear*(w+hin),
            m >= dm.sum(),
            w_lim <= 2*units("in")/(S/b*1.3),
            w <= w_lim*S/b*cbar,
            S/b*cbar*tau >= hin + 2*t + 2*tshear,
            sigma_cfrp >= (beam["\\bar{S}"][:-1]*W_cent*N_max/2/tshear/
                           (S/b*cbar*tau)),
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

        N = 5

        rho_cfrp = Variable("\\rho_{CFRP}", 1.4, "g/cm^3", "density of CFRP")
        rho_foam = Variable("\\rho_{foam}", 0.036, "g/cm^3", "foam density")
        # wing parameters
        S = Variable("S", "ft^2", "wing area")
        wingloading = Variable("W/S", "lbf/ft^2", "Wing loading")

        # Structural evaluation parameters
        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        mtow = Variable("MTOW", "lbf", "Max take off weight")
        m_skin = Variable("m_{skin}", "kg", "Skin mass")
        W = Variable("W", "lbf", "Total wing structural weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        m_fac = Variable("m_{fac}", 1.2, "-", "Wing weight margin factor")
        taucfrp = Variable("\\tau_{CFRP}", 570, "MPa", "torsional stress limit")
        ts = Variable("t_s", "mm", "wing skin thickness")
        Vne = Variable("V_{NE}", 45, "m/s", "never exceed vehicle speed")
        Cmw = Variable("C_{m_w}", 0.121, "-", "negative wing moment coefficent")
        Jtbar = Variable("\\bar{J/t}", 0.01114, "1/mm",
                         "torsional moment of inertia")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        b = Variable("b", "ft", "Span")
        tsmin = Variable("t_{s-min}", 0.012, "in",
                         "minimum wing skin thickness")
        Abar = Variable("\\bar{A}_{jh01}", 0.0753449, "-",
                        "jh01 non dimensional area")
        mfoam = VectorVariable(N-1, "m_{foam}", "kg", "interior mass of wing")

        self.spar = CapSpar(N)
        # self.spar = TubeSpar(5)
        self.submodels = [self.spar]
        # loads normalized by N*W_cent/b

        constraints = [
            m_skin >= rho_cfrp*S*2*ts,
            mfoam >= rho_foam*Abar*(S/b*self.spar["\\bar{c}"])**2*(b/2/4),
            ts >= tsmin,
            wingloading == mtow/S,
            W/m_fac >= m_skin*g + 2*self.spar["m"]*g + 2*mfoam.sum()*g,
            taucfrp >= (1/Jtbar/(S/b*self.spar["\\bar{c}"][0])**2/
                        ts*Cmw*S*rhosl*Vne**2)
            ]

        lc = LinkedConstraintSet([self.spar, constraints], include_only=INCLUDE)

        Model.__init__(self, None, lc, **kwargs)

class FuelTank(Model):
    """
    Returns the weight of the fuel tank.  Assumes a cylinder shape with some
    fineness ratio
    """
    def __init__(self, **kwargs):

        d = Variable("d", "ft", "fuel tank diameter")
        phi = Variable("\\phi", 6, "-", "fuel tank fineness ratio")
        l = Variable("l", "ft", "fuel tank length")
        Stank = Variable("S_{tank}", "ft^2", "fuel tank surface area")
        W = Variable("W", "lbf", "fuel tank weight")
        W_fueltot = Variable("W_{fuel-tot}", "lbf", "Total fuel weight")
        m_fac = Variable("m_{fac}", 1.1, "-", "fuel volume margin factor")
        rho_fuel = Variable("\\rho_{fuel}", 6.01, "lbf/gallon",
                            "density of 100LL")
        rhotank = Variable("\\rho_{fuel-tank}", 0.089, "g/cm^2",
                           "density of plastic fuel tank")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        Voltank = Variable("\\mathcal{V}", "ft^3", "fuel tank volume")

        constraints = [W >= Stank*rhotank*g,
                       Stank/4/phi >= Voltank/l,
                       Voltank/m_fac >= W_fueltot/rho_fuel,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Fuselage(Model):
    """
    Sizes fuselage based off of volume constraints.  Assumes elliptical shape
    """
    def __init__(self, **kwargs):

        # Constants
        d = Variable("d", "ft", "fuselage diameter")
        l = Variable("l", "ft", "fuselage length")

        mskin = Variable("m_{skin}", "kg", "fuselage skin mass")
        rhokevlar = Variable("\\rho_{kevlar}", 1.3629, "g/cm**3",
                             "kevlar density")
        Sfuse = Variable("S_{fuse}", "ft^2", "Fuselage surface area")
        Volavn = Variable("\\mathcal{V}_{avn}", 0.125, "ft^3",
                          "Avionics volume")
        Vol_pay = Variable("\\mathcal{V}_{pay}", 1.0, "ft^3", "Payload volume")
        W = Variable("W", "lbf", "Fuselage weight")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        m_fac = Variable("m_{fac}", 2.1, "-", "Fuselage weight margin factor")
        S_ref = Variable("S_{ref}", "ft**2", "fuselage reference area")
        l_ref = Variable("l_{ref}", "ft", "fuselage reference length")
        rhofoam = Variable("\\rho_{foam400}", 29, "kg/m^3",
                           "density of Foamular 400")
        sigmafoam = Variable("\\sigma_{foam400}", 0.414, "MPa",
                             "shear stress of Foamular 400")
        tfoam = Variable("t_{foam}", "in", "structural foam thickness")
        wbolt = Variable("w_{bolt}", 1, "in", "supporting bolt width")
        Nbolt = Variable("N_{bolt}", 6, "-",
                         "number of bolts from wing to fuselage")
        W_cent = Variable("W_{cent}", "lbf", "Center aircraft weight")
        N_max = Variable("N_{max}", 5, "-", "Load factor")
        mfoam = Variable("m_{foam}", "kg", "mass of structural foam")
        tmin = Variable("t_{min}", 0.03, "in", "minimum skin thickness")
        tskin = Variable("t_{skin}", "in", "skin thickness")
        hengine = Variable("h_{engine}", 6, "in", "engine height")

        ft = FuelTank()

        constraints = [
            mskin >= Sfuse*rhokevlar*tskin,
            tskin >= tmin,
            Sfuse >= pi*d*l + pi*d**2,
            pi*(d/2)**2*l >= ft["\\mathcal{V}"] + l*tfoam*d/2 + Volavn,
            l >= ft["l"],
            d >= hengine,
            W/m_fac >= mskin*g + mfoam*g + ft["W"],
            sigmafoam >= W_cent*N_max/Nbolt/tfoam/wbolt,
            mfoam >= rhofoam*l*tfoam*d/2,
            l_ref == l,
            S_ref == Sfuse
            ]

        lc = LinkedConstraintSet([ft, constraints], include_only=["g"])

        Model.__init__(self, None, lc, **kwargs)

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

        F = VectorVariable(2, "F", "N", "horizontal and vertical tail force")
        T = Variable("T", "N*m", "vertical tail moment")
        L = Variable("L", "ft", "tail boom length")
        E = Variable("E", 150e9, "N/m^2", "young's modulus carbon fiber")
        k = Variable("k", 0.8, "-", "tail boom inertia value")
        kfac = Variable("(1-k/2)", 1-k.value/2, "-", "(1-k/2)")
        I0 = Variable("I_0", "m^4", "tail boom moment of inertia")
        J = Variable("J", "m^4", "tail boom polar moment of inertia")
        d0 = Variable("d_0", "ft", "tail boom diameter")
        t0 = Variable("t_0", "mm", "tail boom thickness")
        tmin = Variable("t_{min}", 0.25, "mm", "minimum tail boom thickness")
        rho_cfrp = Variable("\\rho_{CFRP}", 1.6, "g/cm^3", "density of CFRP")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        W = Variable("W", "lbf", "tail boom weight")
        th = VectorVariable(2, "\\theta", "-", "tail boom deflection angle")
        thmax = Variable("\\theta_{max}", 0.3, "-",
                         "max tail boom deflection angle")
        taucfrp = Variable("\\tau_{CFRP}", 210, "MPa", "torsional stress limit")
        l_ref = Variable("l_{ref}", "ft", "tail boom reference length")
        S_ref = Variable("S_{ref}", "ft**2", "tail boom reference area")
        Vne = Variable("V_{NE}", 40, "m/s", "never exceed vehicle speed")
        rhosl = Variable("\\rho_{sl}", 1.225, "kg/m^3",
                         "air density at sea level")
        Fne = Variable("F_{NE}", "-", "tail boom flexibility factor")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")
        Sh = Variable("S_h", "ft**2", "horizontal tail area")
        Sv = Variable("S_v", "ft**2", "vertical tail surface area")
        bv = Variable("b_v", "ft", "vertical tail span")
        Lmax = Variable("L_{max}", 5.5, "ft", "maximum tail boom length")
        CLmax = Variable("C_{L-max}", 1.39, "-", "Maximum lift coefficient")
        CLhtmax = Variable("C_{L-max_h}", "-",
                           "max lift coefficient of horizontal tail")
        mw = Variable("m_w", 2.0*pi/(1+2.0/23), "-",
                      "assumed span wise effectiveness")

        constraints = [I0 <= pi*t0*d0**3/8.0,
                       W >= pi*g*rho_cfrp*d0*L*t0*kfac,
                       t0 >= tmin,
                       th <= thmax,
                       # L <= Lmax,
                       th >= F*L**2/E/I0*(1+k)/2,
                       CLhtmax/mh >= CLmax/mw,
                       Fne >= 1 + mh*0.5*Vne**2*rhosl*Sh*L**2/E/I0*kfac,
                       F[0] >= 0.5*rhosl*Vne**2*Sh*CLhtmax,
                       F[1] >= 0.5*rhosl*Vne**2*Sv*1.1,
                       T >= 0.5*rhosl*Vne**2*Sv*1.1*bv,
                       taucfrp >= T*d0/2/J,
                       J <= pi/8.0*d0**3*t0,
                       l_ref == L,
                       S_ref == L*pi*d0,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class HorizontalTail(Model):
    def __init__(self, **kwargs):
        Sh = Variable("S_h", "ft**2", "horizontal tail area")
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
        S = Variable("S", "ft^2", "wing area")
        b = Variable("b", "ft", "Span")
        g = Variable("g", 9.81, "m/s^2", "Gravitational acceleration")
        L = Variable("L", "ft", "tail boom length")
        CLhmin = Variable("(C_{L_h})_{min}", 0.75, "-",
                          "max downlift coefficient")
        CLmax = Variable("C_{L-max}", 1.39, "-", "Maximum lift coefficient")
        CM = Variable("C_M", 0.14, "-", "wing moment coefficient")
        l_ref = Variable("l_{ref}", "ft", "horizontal tail reference length")
        S_ref = Variable("S_{ref}", "ft**2", "horizontal tail reference area")

        SMcorr = Variable("SM_{corr}", 0.35, "-", "corrected static margin")
        Fne = Variable("F_{NE}", "-", "tail boom flexibility factor")
        deda = Variable("d\\epsilon/d\\alpha", "-", "wing downwash derivative")
        mw = Variable("m_w", 2.0*pi/(1+2.0/23), "-",
                      "assumed span wise effectiveness")
        mh = Variable("m_h", "-", "horizontal tail span effectiveness")
        cth = Variable("c_{t_h}", "ft", "horizontal tail tip chord")
        lamhfac = Variable("\\lambda_h/(\\lambda_h+1)", 1.0/(1.0+1), "-",
                           "horizontal tail taper ratio factor")

        # signomial helper variables
        sph1 = Variable("sph1", "-", "first term involving $V_h$")
        sph2 = Variable("sph2", "-", "second term involving $V_h$")

        constraints = [
            Vh <= Sh*L/S**2*b,
            bh**2 == ARh*Sh,
            sph1*(mw*Fne/mh/Vh) + deda <= 1,
            sph2 <= Vh*CLhmin/CLmax,
            (sph1 + sph2).mono_lower_bound(
                {"sph1": .48, "sph2": .52}) >= SMcorr + CM/CLmax,
            deda >= mw*S/b/4/pi/L,
            mh*(1+2/ARh) <= 2*pi,
            W >= rhofoam*Sh**2/bh*Abar + g*rhoskin*Sh,
            cth == 2*Sh/bh*lamhfac,
            l_ref == Sh/bh,
            S_ref == Sh
            ]

        Model.__init__(self, None, constraints, **kwargs)

class VerticalTail(Model):
    def __init__(self, **kwargs):
        W = Variable("W", "lbf", "one vertical tail weight")
        Sv = Variable("S_v", "ft**2", "total vertical tail surface area")
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
        S = Variable("S", "ft^2", "wing area")
        b = Variable("b", "ft", "Span")
        L = Variable("L", "ft", "tail boom length")
        l_ref = Variable("l_{ref}", "ft", "vertical tail reference length")
        S_ref = Variable("S_{ref}", "ft**2", "vertical tail reference area")
        ctv = Variable("c_{t_v}", "ft", "vertical tail tip chord")
        lamvfac = Variable("\\lambda_v/(\\lambda_v+1)", 1.0/(1.0+1), "-",
                           "vertical tail taper ratio factor")
        lantenna = Variable("l_{antenna}", 13.4, "in", "antenna length")
        wantenna = Variable("w_{antenna}", 10.2, "in", "antenna width")
        tanlam = Variable("tan(\\Lambda)", np.tan(10*np.pi/180))

        constraints = [Vv <= Sv*L/S/b,
                       bv**2 == ARv*Sv,
                       W >= rhofoam*Sv**2/bv*Abar + g*rhoskin*Sv,
                       ctv == 2*Sv/bv*lamvfac,
                       ctv >= wantenna + bv*tanlam,
                       bv >= lantenna,
                       S_ref == Sv,
                       l_ref == Sv/bv,
                      ]

        Model.__init__(self, None, constraints, **kwargs)

class Empennage(Model):
    def __init__(self, **kwargs):
        m_fac = Variable("m_{fac}", 1.0, "-", "Tail weight margin factor")
        W = Variable("W", "lbf", "empennage weight")

        tb = TailBoom()
        ht = HorizontalTail()
        vt = VerticalTail()
        self.submodels = [tb, ht, vt]

        constraints = [
            W/m_fac >= tb["W"] + ht["W"] + vt["W"],
            ]

        lc = LinkedConstraintSet(
            [self.submodels, constraints],
            include_only=["L", "F_{NE}", "m_h", "S_h", "S", "b", "b_v", "b_h",
                          "S_v", "m_w"]
            )

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
                                 [W_fueltot, W_pay, W_skid]),
            SummingConstraintSet(W_zfw, "W", zf_loads,
                                 [W_pay, W_skid])
            ]

        ConstraintSet.__init__(self, constraints, **kwargs)


def return_cg(weights, xlocs):
    assert len(weights) == len(xlocs)
    xcg = sum(list(map(mul, weights, xlocs)))/sum(weights)
    return xcg

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

        engineweight = EngineWeight(DF70)
        wing = Wing()
        empennage = Empennage()
        avionics = Avionics()
        fuselage = Fuselage()
        dragcomps = [fuselage] + empennage.submodels
        mission = Mission(h_station, wind, DF70, Nclimb, Nloiter, dragcomps)
        center_loads = [fuselage, avionics, engineweight]
        zf_loads = center_loads + [empennage, wing]
        weight = Weight(center_loads, zf_loads)

        self.submodels = zf_loads + [weight, mission]

        constraints = []

        if margin:
            sq = SystemRequirements()
            self.submodels.append(sq)

        cs = ConstraintSet([self.submodels, constraints])

        linked = {}
        for c in dragcomps:
            for name in ["l_{ref}", "S_{ref}"]:
                vks = cs.varkeys[name]
                for vk in vks:
                    descr = dict(vk.descr)
                    if c.name in descr["models"]:
                        descr.pop("value", None)
                        descr["models"] = [c.name]
                        descr["modelnums"] = [c.num]
                        newvk = VarKey(**descr)
                        linked.update({vk: newvk})
        cs.subinplace(linked)

        lc = LinkedConstraintSet([cs, constraints],
                                 include_only=INCLUDE)


        objective = 1/mission["t_{loiter}"]

        Model.__init__(self, objective, lc, **kwargs)

    # def process_solution(self, sol):
    #     xwing = 0.5*(sol("l_Fuselage, GasMALE")+sol("d"))*0.93 # 0.93 for 10lb
    #     self.xwing = xwing
    #     xnp = xwing + (sol("m_h")/sol("m_w")*(1.0-4.0/(sol("AR")+2.0))*
    #                    sol("V_h"))*sol("S")/sol("b")
    #     weights = [sol("W_Fuselage, GasMALE"),
    #                sol("W_{pay}")*2.5,
    #                sol("W_EngineWeight, GasMALE"),
    #                sol("W_{fuel-tot}"),
    #                sol("W_Wing, GasMALE"),
    #                sol("W_TailBoom, Empennage, GasMALE"),
    #                sol("W_HorizontalTail, Empennage, GasMALE"),
    #                sol("W_VerticalTail, Empennage, GasMALE")
    #               ]

    #     xlocs = [0.5*(sol("l_Fuselage, GasMALE") + sol("d")),
    #              0.25*sol("d"),
    #              sol("l_Fuselage, GasMALE") + 0.75*sol("d"),
    #              0.5*(sol("l_Fuselage, GasMALE") + sol("d")),
    #              xwing,
    #              xwing + 0.5*sol("L"),
    #              xwing + sol("L"),
    #              xwing + sol("L")
    #             ]

    #     xcg = return_cg(weights, xlocs)
    #     SM = (xnp - xcg)/(sol("S")/sol("b"))

    #     self.xnp = xnp
    #     self.xcg = xcg
    #     self.SM = SM

    #     # wpay = sol("W_{pay}")*np.linspace(0.5, 2.5, 15)
    #     # ind = weights.index(sol("W_{pay}"))
    #     # cgs = []
    #     # SMs = []
    #     # xws = []
    #     # for w in wpay:
    #     #     weights[ind] = w
    #     #     cg = return_cg(weights, xlocs)
    #     #     cgs.extend([cg])
    #     #     SMs.extend([((xnp-cg)/sol("S")*sol("b")).magnitude])
    #     #     xws.extend([(SM*sol("S")/sol("b")+cg-xnp+xwing).magnitude])

    #     # weights[ind] = sol("W_{pay}")

    #     # fig, ax = plt.subplots(2)
    #     # ax[0].plot(wpay, SMs)
    #     # ax[1].plot(wpay, xws)
    #     # ax[0].set_ylabel("Static Margin")
    #     # ax[1].set_ylabel("Wing location [ft]")
    #     # ax[0].set_xticks([])
    #     # ax[1].set_xlabel("Payload weight [lbf]")
    #     # ax[0].set_title("Static Margin Full Fuel")
    #     # fig.savefig("smvswpay.pdf")

    #     indf = weights.index(sol("W_{fuel-tot}"))
    #     del weights[indf]
    #     del xlocs[indf]
    #     xcg = return_cg(weights, xlocs)
    #     SM = (xnp - xcg)/(sol("S")/sol("b"))
    #     self.xcgempty = xcg
    #     self.SMempty = SM
    #     # cgs = []
    #     # SMs = []
    #     # xws = []
    #     # for w in wpay:
    #     #     weights[ind] = w
    #     #     cg = return_cg(weights, xlocs)
    #     #     cgs.extend([cg])
    #     #     SMs.extend([((xnp-cg)/sol("S")*sol("b")).magnitude])
    #     #     xws.extend([(SM*sol("S")/sol("b")+cg-xnp+xwing).magnitude])

    #     # fig, ax = plt.subplots(2)
    #     # ax[0].plot(wpay, SMs)
    #     # ax[1].plot(wpay, xws)
    #     # ax[0].set_ylabel("Static Margin")
    #     # ax[1].set_ylabel("Wing location [ft]")
    #     # ax[0].set_xticks([])
    #     # ax[1].set_xlabel("Payload weight [lbf]")
    #     # ax[0].set_title("Static Margin Empty Fuel")
    #     # fig.savefig("smvswpayempty.pdf")

    #     # indh = weights.index(sol("W_HorizontalTail, Empennage, GasMALE"))
    #     # indv = weights.index(sol("W_VerticalTail, Empennage, GasMALE"))
    #     # del weights[indv]
    #     # del weights[indh]
    #     # del xlocs[indv]
    #     # del xlocs[indh]
    #     # weights[ind] = sol("W_{pay}")
    #     # wtail = []
    #     # for w in wpay:
    #     #     weights[ind] = w
    #     #     wtail.append(((self.xcg*sum(weights)-sum(list(map(mul, weights, xlocs))))/(xwing+sol("L"))).magnitude)

    #     # fig, ax = plt.subplots()
    #     # ax.plot(wpay, wtail)
    #     # ax.set_xlabel("Payload weight [lbf]")
    #     # ax.set_ylabel("Tail weight [lbf]")
    #     # fig.savefig("wtailvswpay.pdf")

    # def get_cg(self):
    #     """ return cg location at full and empty """
    #     return self.xcg, self.xcgempty

    # def get_np(self):
    #     """ return neutral point """
    #     return self.xnp

    # def get_SM(self):
    #     """ return static margin at full and empty """
    #     return self.SM, self.SMempty

    # def get_wingloc(self):
    #     """ return wing location"""
    #     return self.xwing

if __name__ == "__main__":
    M = GasMALE(DF70=True)
    M.substitutions.update({"t_{loiter}": 6})
    # M.substitutions.update({"b": 24})
    M.cost = M["MTOW"]
    if SIGNOMIALS:
        sol = M.localsolve("mosek")
    else:
        sol = M.solve("mosek")
