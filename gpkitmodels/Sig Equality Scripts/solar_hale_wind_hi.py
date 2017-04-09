from gpkit import Variable, Model, units, SignomialsEnabled, SignomialEquality
from gpkit.nomials import SignomialEqualityTriv, SignomialEqualityTrivTrust, SignomialEqualityLinTrust, SignomialEqualityLin
from gpkit.constraints.set import ConstraintSet
from gpkit import LinkedConstraintSet
import numpy as np
import matplotlib.pyplot as plt

CD = Variable('C_D', '-', 'Drag coefficient')
CL = Variable('C_L', '-', 'Lift coefficient')
P_shaft = Variable('P_{shaft}', 'W', 'Shaft power')
S = Variable('S', 'm^2', 'Wing reference area')
S_solar = Variable('S_{solar}', 'm^2', 'solar cell area')
V = Variable('V', 'm/s', 'Cruise velocity')
W = Variable('W', 'lbf', 'Aircraft weight')
# rho = Variable(r'\rho', 0.1895, 'kg/m^3', 'Density of air')
rho = Variable(r'\rho', 'kg/m^3', 'Density of air')
E_batt = Variable('E_{batt}', 'J', 'Battery energy')
h = Variable("h", "ft", "Altitude")
g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
AR = Variable('AR', 27, '-', "aspect ratio")

# Steady Level Flight

class SteadyLevelFlight(Model):
    def __init__(self, **kwargs):
        eta_prop = Variable(r'\eta_{prop}', 0.80, '-', 'Propulsive efficiency')

        constraints = [P_shaft >= V*W*CD/CL/eta_prop,   # eta*P = D*V
                       W == 0.5*rho*V**2*CL*S]
        Model.__init__(self, None, constraints, **kwargs)

# Aerodynamics

class Aero(Model):
    def __init__(self, **kwargs):

        Cd0 = Variable('C_{d0}', 0.002, '-', "non-wing drag coefficient")
        cdp = Variable("c_{dp}", "-", "wing profile drag coeff")
        CLmax = Variable('C_{L-max}', 1.5, '-', 'maximum lift coefficient')
        e = Variable('e', 0.9, '-', "spanwise efficiency")
        b = Variable('b', 'ft', 'span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', "dynamic viscosity")
        Re = Variable("Re", '-', "Reynolds number")
        Re_ref = Variable("Re_{ref}", 3e5, "-", "Reference Re for cdp")
        Cf = Variable("C_f", "-", "wing skin friction coefficient")
        Kwing = Variable("K_{wing}", 1.3, "-", "wing form factor")

        constraints = [
            CD >= Cd0 + cdp + CL**2/(np.pi*e*AR),
            cdp >= ((0.006 + 0.005*CL**2 + 0.00012*CL**10)*(Re/Re_ref)**-0.3),
            b**2 == S*AR,
            CL <= CLmax,
            Re == rho*V/mu*(S/AR)**0.5,
            ]
        Model.__init__(self, None, constraints, **kwargs)

# Weights

class Weight(Model):
    def __init__(self, **kwargs):

        W_batt = Variable('W_{batt}', 'lbf', 'Battery weight')
        W_airframe = Variable('W_{airframe}', 'lbf', 'Airframe weight')
        W_solar = Variable('W_{solar}', 'lbf', 'Solar panel weight')
        W_pay = Variable('W_{pay}', 10, 'lbf', 'Payload weight')
        W_avionics = Variable('W_{avionics}', 8, 'lbf', 'avionics weight')
        rho_solar = Variable(r'\rho_{solar}', 0.3, 'kg/m^2',
                             'Solar cell area density')
        f_airframe = Variable('f_{airframe}', 0.30, '-',
                              'Airframe weight fraction') # incl propulsion
        h_batt = Variable('h_{batt}', 350, 'W*hr/kg', 'Battery energy density')

	coeff_airframe = Variable('coeff_af', 0.44, 'lbf/m^3.1', 'coefficient for airframe')

        constraints = [#W_airframe >= coeff_airframe*S**1.55*AR**1.30,
		       W_airframe >= W*f_airframe,
                       W_batt >= E_batt/h_batt*g,
                       W_solar >= rho_solar*g*S_solar,
                       W >= W_pay + W_solar + W_airframe + W_batt + W_avionics]
        Model.__init__(self, None, constraints, **kwargs)

# Power

class Power(Model):
    def __init__(self, **kwargs):

        ES_irr = Variable('(E/S)_{irr}', 2.9, 'kW*hr/m^2',
                          'Average daytime solar energy')
        P_oper = Variable('P_{oper}', 'W', 'Aircraft operating power')
        P_acc = Variable('P_{acc}', 25, 'W', 'Accessory power draw')
        eta_solar = Variable(r'\eta_{solar}', 0.2, '-',
                             'Solar cell efficiency')
        eta_charge = Variable(r'\eta_{charge}', 0.95, '-',
                              'Battery charging efficiency')
        eta_discharge = Variable(r'\eta_{discharge}', 0.95, '-',
                                 'Battery discharging efficiency')
        t_day = Variable('t_{day}', 8, 'hr', 'Daylight span')
        t_night = Variable('t_{night}', 16, 'hr', 'Night span')

        constraints = [ES_irr*eta_solar*S_solar >= P_oper*t_day + E_batt/eta_charge,
                       S_solar <= S,
                       P_oper >= P_shaft + P_acc,
                       E_batt >= P_oper*t_night/eta_discharge]
        Model.__init__(self, None, constraints, **kwargs)

# Atmosphere

class Atmosphere(Model):
    def __init__(self, **kwargs):

        p_sl = Variable("p_{sl}", 101325, "Pa", "Pressure at sea level")
        T_sl = Variable("T_{sl}", 288.15, "K", "Temperature at sea level")
        L_atm = Variable("L_{atm}", 0.0065, "K/m", "Temperature lapse rate")
        T_atm = Variable("T_{atm}", "K", "air temperature")
        M_atm = Variable("M_{atm}", 0.0289644, "kg/mol",
                         "Molar mass of dry air")
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K",
                         "air specific heating value")
        p_atm = Variable("P_{atm}", "Pa", "air pressure")
        TH = (g*M_atm/R_atm/L_atm).value

        with SignomialsEnabled():
            constraints = [
                h <= 100000*units.m,
                h >= 1E-6*units.m,
                #h <= 20000*units.m,  # Model valid to top of troposphere

                # Pressure-altitude relation
                (p_atm/p_sl)**(1/TH) == T_atm/T_sl,

                # Ideal gas law
                rho == p_atm/(R_atm/M_atm*T_atm),

                # T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                SignomialEquality(T_sl, T_atm + L_atm*h)]
                # http://en.wikipedia.org/wiki/Density_of_air#Altitude
                #rho <= p_sl*T_atm**(TH-1)*M_atm/R_atm/(T_sl**TH)]
	Model.__init__(self, None, constraints, **kwargs)

class StationKeeping(Model):
    def __init__(self, **kwargs):

        h_min = Variable('h_{min}', 11000, 'm', 'minimum altitude for GP fit')
        h_max = Variable('h_{max}', 20000, 'm', 'maximum altitude for GP fit')
        V_wind = Variable('V_{wind}', 26, 'm/s', 'wind speed')

        h_0 = Variable('h_0',1,'m','h0 for nondimensionalization')
        V_0 = Variable('V_0',1,'m/s','V0 for nondimensionalization')

	with SignomialsEnabled():
	    constraints = [h >= h_min,
                       (2.63 * (h/h_0)**-0.162
                        + 2.68 * (h/h_0)**-0.156
                        + 2.37 * (h/h_0)**-0.178)/(V/V_0)**(0.131) <= 1,
                       h <= h_max] # units?
        Model.__init__(self, None, constraints, **kwargs)

class SolarHALE(Model):
    """High altitude long endurance solar UAV"""
    def __init__(self, **kwargs):
        """Setup method should return objective, list of constraints"""

        slf = SteadyLevelFlight()
        power = Power()
        weight = Weight()
        sk = StationKeeping()
        aero = Aero()
        atmosphere = Atmosphere()
        self.submodels = [slf, power, weight, sk, aero, atmosphere]
        # self.submodels = [slf, power, weight, sk, aero]

        constraints = []
        lc = LinkedConstraintSet([self.submodels, constraints])

        objective = W

        Model.__init__(self, objective, lc, **kwargs)

if __name__ == "__main__":
    M = SolarHALE()
    sol = M.localsolve("mosek",x0={h:60000*units.ft},iteration_limit=500)
    print sol.table()
