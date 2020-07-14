import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, SignomialEquality, VarKey, units, Vectorize

# Importing atmospheric model
from gpkitmodels.SP.atmosphere.atmosphere import Atmosphere

# SimPleAC with mission design and flight segments, and lapse rate and BSFC model (3.4.2)

class SimPleAC(Model):
    def setup(self):
        self.engine = Engine()
        self.wing = Wing()
        self.fuse = Fuselage()
        self.components = [self.engine, self.wing, self.fuse]

        # Environmental constants
        g         = Variable("g", 9.81, "m/s^2", "gravitational acceleration")
        rho_f     = Variable("\\rho_f", 817, "kg/m^3", "density of fuel")

        # Free Variables
        W         = Variable("W", "N", "maximum takeoff weight")
        W_f       = Variable("W_f", "N", "maximum fuel weight")
        V_f       = Variable("V_f", "m^3", "maximum fuel volume")
        V_f_avail = Variable("V_{f_{avail}}", "m^3", "fuel volume available")

        constraints = []

        # Fuel volume model
        with SignomialsEnabled():
            constraints += [V_f == W_f / g / rho_f,
                    V_f_avail <= self.wing['V_{f_{wing}}'] + self.fuse['V_{f_{fuse}}'], #[SP]
                    V_f_avail >= V_f]

        return constraints, self.components

    def dynamic(self,state):
        return SimPleACP(self,state)

class SimPleACP(Model):
    def setup(self,aircraft,state):
        self.aircraft = aircraft
        self.engineP  = aircraft.engine.dynamic(state)
        self.wingP    = aircraft.wing.dynamic(state)
        self.fuseP    = aircraft.fuse.dynamic(state)
        self.Pmodels  = [self.engineP, self.wingP, self.fuseP]

        # Free variables
        C_D       = Variable("C_D", "-", "drag coefficient")
        D         = Variable("D", "N", "total drag force")
        LoD       = Variable('L/D','-','lift-to-drag ratio')
        V         = Variable("V", "m/s", "cruising speed")

        constraints = []

        constraints += [self.engineP['T'] * V <= self.aircraft.engine['\\eta_{prop}'] * self.engineP['P_{shaft}'],
                    C_D >= self.fuseP['C_{D_{fuse}}'] + self.wingP['C_{D_{wpar}}'] + self.wingP['C_{D_{ind}}'],
                    D >= 0.5 * state['\\rho'] * self.aircraft['S'] * C_D * V ** 2,
                    self.wingP['Re'] == (state['\\rho'] / state['\\mu']) * V * (self.aircraft['S'] / self.aircraft['A']) ** 0.5,
                    self.fuseP['Re_{fuse}'] == state['\\rho']*V*self.aircraft.fuse['l_{fuse}']/state['\\mu'],
                    LoD == self.wingP['C_L'] / C_D]

        return constraints, self.Pmodels

class Fuselage(Model):
    def setup(self):
        # Free Variables
        S         = Variable('S_{fuse}', 'm^2', 'fuselage surface area')
        l         = Variable('l_{fuse}', 'm', 'fuselage length', fix = True)
        r         = Variable('r_{fuse}', 'm', 'fuselage minor radius')
        f         = Variable('f_{fuse}', '-', 'fuselage fineness ratio', fix = True)
        k         = Variable('k_{fuse}', '-', 'fuselage form factor')
        # Free variables (fixed for performance eval.)
        V         = Variable('V_{fuse}', 'm^3', 'total volume in the fuselage', fix = True)
        V_f_fuse  = Variable('V_{f_{fuse}}', 'm^3', 'fuel volume in the fuselage')
        W_fuse    = Variable('W_{fuse}', 'N', 'fuselage weight')
        p = 1.6075

        constraints = [f == l/r/2,
                       f <= 6,
                       k >= 1 + 60/f**3 + f/400,
                        3*(S/np.pi)**p >= 2*(l*2*r)**p + (2*r)**(2*p),
                        V == 4./6.*np.pi*r**2*l,
                        V_f_fuse >= 1*10**-10*units('m^3'),
                       ]

        return constraints

    def dynamic(self,state):
        return FuselageP(self,state)

class FuselageP(Model):
    def setup(self,fuselage,state):
        # Constants
        Cfref      = Variable('C_{f_{fuse,ref}}', 0.455, '-', 'fuselage reference skin friction coefficient', pr=10.)
        # Free Variables
        Re         = Variable('Re_{fuse}', '-', 'fuselage Reynolds number')
        Cf         = Variable('C_{f_{fuse}}', '-', 'fuselage skin friction coefficient')
        Cd         = Variable('C_{D_{fuse}}', '-', 'fuselage drag coefficient')

        constraints = [Cf >= Cfref/Re**0.3,
                       Cd >= fuselage['k_{fuse}']*Cf,
                       ]
        return constraints


class Wing(Model):
    def setup(self):
        # Non-dimensional constants
        C_Lmax     = Variable("C_{L,max}", 1.6, "-", "lift coefficient at stall", pr=5.)
        e          = Variable("e", 0.92, "-", "Oswald efficiency factor", pr=3.)
        N_ult      = Variable("N_{ult}", 3, "-", "ultimate load factor", pr=15.)
        tau        = Variable("\\tau", "-", "airfoil thickness to chord ratio", fix = True)
        tau_ref    = Variable("\\tau_{ref}", 0.12, "-", "reference airfoil thickness to chord ratio")

        # Dimensional constants
        W_w_coeff1 = Variable("W_{w_{coeff1}}", 2e-5, "1/m",
                              "wing weight coefficient 1", pr= 30.) #orig  12e-5
        W_w_coeff2 = Variable("W_{w_{coeff2}}", 60., "Pa",
                              "wing weight coefficient 2", pr=10.)

        # Free Variables (fixed for performance eval.)
        A         = Variable("A", "-", "aspect ratio",fix = True)
        S         = Variable("S", "m^2", "total wing area", fix = True)
        W_w       = Variable("W_w", "N", "wing weight")
        W_w_strc  = Variable('W_{w_{strc}}','N','wing structural weight', fix = True)
        W_w_surf  = Variable('W_{w_{surf}}','N','wing skin weight', fix = True)
        V_f_wing  = Variable("V_{f_{wing}}",'m^3','fuel volume in the wing', fix = True)

        constraints = []

        # Structural model
        constraints += [W_w_surf >= W_w_coeff2 * S,
                        W_w >= W_w_surf + W_w_strc]

      # Wing fuel and form factor model
        constraints += [V_f_wing**2 <= 0.0009*S**3/A*tau**2, # linear with b and tau, quadratic with chord
                        tau >= 0.08, tau <= 0.23,
                        ]

        # Form factor model

        return constraints

    def dynamic(self,state):
        return WingP(self,state)

class WingP(Model):
    def setup(self,wing,state):
        self.wing = wing
        # Free Variables
        C_D_ind   = Variable('C_{D_{ind}}', '-', "wing induced drag coefficient")
        C_D_wpar  = Variable('C_{D_{wpar}}', '-', "wing profile drag coefficient")
        C_L       = Variable("C_L", "-", "wing lift coefficient")
        Re        = Variable("Re", "-", "Reynolds number")
        Re_ref    = Variable("Re_{ref}", 1500000, "-", "reference Reynolds number")

        constraints = []

        # Drag model
        w = C_D_wpar
        u_1 = C_L
        u_2 = Re/Re_ref
        u_3 = self.wing['\\tau']/self.wing['\\tau_{ref}']
        nc = w**0.00488697 >= 0.000347324 * (u_1)**6.64787 * (u_2)**-0.00842527 * (u_3)**-0.406817 + \
                            0.974515 * (u_1)**-0.00206058 * (u_2)**-0.00117649 * (u_3)**-0.000597604 + \
                            0.000211504 * (u_1)**1.35483 * (u_2)**-0.252459 * (u_3)**3.91243
        nc.name = 'drag'
        constraints += [C_D_ind == C_L ** 2 / (np.pi * self.wing['A'] * self.wing['e']),
                        nc]

        return constraints

class Engine(Model):
    def setup(self):
        # Dimensional constants
        BSFC_ref    = Variable("BSFC_{ref}", 0.32, "lbf/(hp*hr)", "reference brake specific fuel consumption")
        eta_prop    = Variable("\\eta_{prop}", 0.8, '-',"propeller efficiency")
        P_shaft_ref = Variable("P_{shaft,ref}", 10, "hp", "reference MSL maximum shaft power")
        W_e_ref     = Variable("W_{e,ref}", 10, "lbf","reference engine weight")
        h_ref       = Variable("h_{ref}", 15000,'ft','engine lapse reference altitude')

        # Free variables
        P_shaft_max = Variable("P_{shaft,max}","kW","MSL maximum shaft power")
        W_e         = Variable("W_e", "N", "engine weight", fix = True)

        constraints = [(W_e/W_e_ref) == 1.27847 * (P_shaft_max/P_shaft_ref)**0.772392]
        return constraints

    def dynamic(self,state):
        return EngineP(self,state)

class EngineP(Model):
    def setup(self,engine,state):
        self.engine = engine
        # Dimensional constants

        # Free variables
        BSFC        = Variable("BSFC", "lbf/(hp*hr)", "brake specific fuel consumption")
        P_shaft     = Variable("P_{shaft}", "kW","shaft power")
        P_shaft_alt = Variable("P_{shaft,alt}", "kW", 'maximum shaft power at altitude')
        Thrust      = Variable("T", "N", "propeller thrust")

        L           = Variable("L","-","power lapse percentage")

        constraints = []

        with SignomialsEnabled():
            constraints += [P_shaft <= P_shaft_alt,
                        L == (0.937 * (state['h']/self.engine['h_{ref}'])**0.0922)**10,
                        SignomialEquality(1, L + P_shaft_alt / self.engine['P_{shaft,max}']),
                        (BSFC/self.engine['BSFC_{ref}'])**(0.1) >= 0.984*(P_shaft/P_shaft_alt)**-0.0346,
                        BSFC/self.engine['BSFC_{ref}'] >= 1.,
                        ]
        return constraints


class Mission(Model):
    def setup(self,aircraft,Nsegments):
        self.aircraft = aircraft
        W_f_m   = Variable('W_{f_m}','N','total mission fuel')
        t_m     = Variable('t_m','hr','total mission time')

        with Vectorize(Nsegments):
            Wavg    = Variable('W_{avg}','N','segment average weight')
            Wstart  = Variable('W_{start}', 'N', 'weight at the beginning of flight segment')
            Wend    = Variable('W_{end}', 'N', 'weight at the end of flight segment')
            h       = Variable('h','m','final segment flight altitude')
            havg    = Variable('h_{avg}','m','average segment flight altitude')
            dhdt    = Variable('\\frac{dh}{dt}','m/hr','climb rate')
            W_f_s   = Variable('W_{f_s}','N', 'segment fuel burn')
            t_s     = Variable('t_s','hr','time spent in flight segment')
            R_s     = Variable('R_s','km','range flown in segment')
            state   = Atmosphere()
            self.aircraftP = self.aircraft.dynamic(state)

        # Mission variables
        hcruise    = Variable('h_{cruise_m}', 'm', 'minimum cruise altitude')
        Range      = Variable("Range_m", "km", "aircraft range")
        W_p        = Variable("W_{p_m}", "N", "payload weight", pr=20.)
        rho_p      = Variable("\\rho_{p_m}", "kg/m^3", "payload density", pr = 10.)
        V_min      = Variable("V_{min_m}", "m/s", "takeoff speed", pr=20.)
        TOfac      = Variable('T/O factor_m', '-','takeoff thrust factor')
        cost_index = Variable("C_m", '1/hr','hourly cost index')

        constraints = []

        # Setting up the mission
        with SignomialsEnabled():
            constraints += [havg == state['h'], # Linking states
                        h[1:Nsegments-1] >= hcruise,  # Adding minimum cruise altitude

                        # Weights at beginning and end of mission
                        Wstart[0] >= W_p + self.aircraft.wing['W_w'] + self.aircraft.engine['W_e'] + self.aircraft.fuse['W_{fuse}'] + W_f_m,
                        Wend[Nsegments-1] >= W_p + self.aircraft.wing['W_w'] + self.aircraft.engine['W_e'] + self.aircraft.fuse['W_{fuse}'],

                        # Lift, and linking segment start and end weights
                        Wavg <= 0.5 * state['\\rho'] * self.aircraft['S'] * self.aircraftP.wingP['C_L'] * self.aircraftP['V'] ** 2,
                        Wstart >= Wend + W_f_s, # Making sure fuel gets burnt!
                        Wstart[1:Nsegments] == Wend[:Nsegments-1],
                        Wavg == Wstart ** 0.5 * Wend ** 0.5,

                        # Altitude changes
                        h[0] == t_s[0]*dhdt[0], # Starting altitude
                        dhdt >= 1.*units('m/hr'),
                        havg[0] == 0.5*h[0],
                        havg[1:Nsegments] == (h[1:Nsegments]*h[0:Nsegments-1])**(0.5),
                        SignomialEquality(h[1:Nsegments],h[:Nsegments-1] + t_s[1:Nsegments]*dhdt[1:Nsegments]),

                        # Thrust and fuel burn
                        W_f_s >= self.aircraftP.engineP['BSFC'] * self.aircraftP.engineP['P_{shaft}'] * t_s,
                        self.aircraftP.engineP['T'] * self.aircraftP['V'] >= self.aircraftP['D'] * self.aircraftP['V'] + Wavg * dhdt,

                        # Max MSL thrust at least 2*climb thrust
                        self.aircraft.engine['P_{shaft,max}'] >= TOfac*self.aircraftP.engineP['P_{shaft}'][0],

                        # Flight time
                        t_s == R_s/self.aircraftP['V'],

                        # Aggregating segment variables
                        self.aircraft['W_f'] >= W_f_m,
                        R_s == Range/Nsegments, # Dividing into equal range segments
                        W_f_m >= sum(W_f_s),
                        t_m >= sum(t_s)
                        ]

        # Maximum takeoff weight
        constraints += [self.aircraft['W'] >= W_p + self.aircraft.wing['W_w'] + self.aircraft['W_f'] +
                        self.aircraft.engine['W_e'] + self.aircraft.fuse['W_{fuse}']]

        # Stall constraint
        constraints += [self.aircraft['W'] <= 0.5 * state['\\rho'] *
                            self.aircraft['S'] * self.aircraft['C_{L,max}'] * V_min ** 2]

        # Wing weight model
        constraints += [self.aircraft.wing['W_{w_{strc}}']**2. >=
                        self.aircraft.wing['W_{w_{coeff1}}']**2. / self.aircraft.wing['\\tau']**2. *
                        (self.aircraft.wing['N_{ult}']**2. * self.aircraft.wing['A'] ** 3. *
                        ((W_p + self.aircraft.fuse['W_{fuse}'] +
                          self.aircraft['W_e'] + self.aircraft.fuse['V_{f_{fuse}}']*self.aircraft['g']*self.aircraft['\\rho_f']) *
                         self.aircraft['W'] * self.aircraft.wing['S']))]

        # Fuselage volume and weight
        constraints += [self.aircraft.fuse['V_{fuse}'] >=
                        self.aircraft.fuse['V_{f_{fuse}}'] + W_p/(rho_p*self.aircraft['g']),
                        self.aircraft.fuse['W_{fuse}'] == self.aircraft.fuse['S_{fuse}']*self.aircraft.wing['W_{w_{coeff2}}'],
                        ]

        # Upper bounding variables
        constraints += [t_m <= 100000*units('hr'),
            W_f_m <= 1e10*units('N'),
            cost_index >= 1e-10*units('1/hr')]

        return constraints, state, self.aircraft, self.aircraftP

def test():
    m = Mission(SimPleAC(),4)
    m.substitutions.update({
        'h_{cruise_m}'   :5000*units('m'),
        'Range_m'        :3000*units('km'),
        'W_{p_m}'        :3000*units('N'),
        '\\rho_{p_m}'    :1500*units('kg/m^3'),
        'C_m'            :120*units('1/hr'),
        'V_{min_m}'      :35*units('m/s'),
        'T/O factor_m'   :2,
    })
    m.cost = m['W_{f_m}']*units('1/N') + m['C_m']*m['t_m']
    sol = m.localsolve(verbosity=0)

if __name__ == "__main__":
    test()
