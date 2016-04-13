# SOLAR HIGH ALTITUDE LONG ENDURANCE AIRCRAFT

This paper is a presentation of the infeasibility of a solar powered aircraft for the given requirements presented in table 1. This study will prove that the amount of sunlight accessible to the aircraft for the 45 degree latitude is not enough to meet the 90% availablility requirement. The study was done using a convex optimization tool called, GPkit.  The optimization model was constructed using various models to describe the physics required to achieve the specified requirements. The specific requirements to which this plane was optimized are found in table 1. 

\begin{table}[H]
\begin{center}
\label{t:battery}
\begin{tabular}{ | c |  c | }
    \hline
    Requirement & Specification \\ \hline
    Payload & 10 [lbs] 10 [Watts] \\\hline
    Endurance & 5-40 days \\\hline
    Availability & 90\% \\\hline
    Latitude Coverage & 45-60 degrees \\\hline
\end{tabular}
\caption{Table outlining the requirements given to us by the customer.}
\end{center}
\end{table}

# Assumptions

\begin{enumerate}
    \item We assumed steady level flight. 
    \item We assumed a standard aerodynamic model using profile drag and induced drag. 
    \item We assumed a constant fractional airframe weight ($W_{airframe} = W_{plane}f_{airframe}$), where $f_{airframe} = 20\%$.  The fractional weight includes the weight of the wing, fuselage, engine (with parts), tail, and any other structural elements. 
    \item The avionics weight was assumed to be 6 lbs and draws 15 watts of power and is based on required avionics for flight control, ground communication, and satellite communication.
    \item The solar cell area density is based off of the best available solar cell density. 
    \item We assumed that the payload weighs 10 lbs and draws 10 watts of power. 
    \item We assumed that we are flying at 45 degree latitude and during the winter solstice. By being able to fly during the winter solstice we guarantee that we can fly during any other day of the year.  This assumptions translates to 8 hours of charging time and an angle of incidence of about 66 degrees. 
    \item We assume the objective function to be weight becasue aircraft cost typically scale with weight.
\end{enumerate}

# Station Keeping Requirements

To achieve the minimum footprint of 100km, the aircraft must fly at 15,000 ft. In order to avoid diminished solar irradiance due to weather, and to avoid air traffic conflicts, the aircraft must fly at 50,000 ft or above. 

\begin{figure}[h!]
	\label{f:windspeeds}
	\begin{center}
	\includegraphics[scale = 0.6]{windspeeds}
	\caption{This figure shows the wind speeds necessary to station keep at the 90\%, 95\% and the 99\% winds for a given altitude.  As can be seen, the winds at 50,000 ft are higher than those at 15,000 ft.} 
	\end{center}
\end{figure}


```python
#inPDF: skip

from gpkit import Variable, Model, units
from gpkit.constraints.set import ConstraintSet
from gpkit import LinkedConstraintSet
import numpy as np
import matplotlib.pyplot as plt

PLOT = False 
        
CD = Variable('C_D', '-', 'Drag coefficient')
CL = Variable('C_L', '-', 'Lift coefficient')
P_shaft = Variable('P_{shaft}', 'W', 'Shaft power')
S = Variable('S', 'm^2', 'Wing reference area')
V = Variable('V', 'm/s', 'Cruise velocity')
W = Variable('W', 'lbf', 'Aircraft weight')
rho = Variable(r'\rho', 'kg/m^3', 'Density of air')
E_batt = Variable('E_{batt}', 'J', 'Battery energy')
h = Variable("h", "ft", "Altitude")
g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

# Steady Level Flight

class SteadyLevelFlight(Model):
    def __init__(self, **kwargs):
        eta_prop = Variable(r'\eta_{prop}', 0.7, '-', 'Propulsive efficiency')

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
        AR = Variable('AR', 27, '-', "aspect ratio")
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
        W_pay = Variable('W_{pay}', 4, 'lbf', 'Payload weight')
        W_avionics = Variable('W_{avionics}', 4, 'lbf', 'avionics weight')
        rho_solar = Variable(r'\rho_{solar}', 1.2, 'kg/m^2',
                             'Solar cell area density')
        f_airframe = Variable('f_{airframe}', 0.20, '-',
                              'Airframe weight fraction')
        h_batt = Variable('h_{batt}', 250, 'W*hr/kg', 'Battery energy density')

        constraints = [W_airframe >= W*f_airframe,
                       W_batt >= E_batt/h_batt*g,
                       W_solar >= rho_solar*g*S,
                       W >= W_pay + W_solar + W_airframe + W_batt + W_avionics]
        Model.__init__(self, None, constraints, **kwargs)

# Power 

class Power(Model):
    def __init__(self, **kwargs):

        PS_irr = Variable('(P/S)_{irr}', 1000*0.5, 'W/m^2',
                          'Average daytime solar irradiance')
        P_oper = Variable('P_{oper}', 'W', 'Aircraft operating power')
        P_charge = Variable('P_{charge}', 'W', 'Battery charging power')
        P_acc = Variable('P_{acc}', 25, 'W', 'Accessory power draw')
        eta_solar = Variable(r'\eta_{solar}', 0.2, '-',
                             'Solar cell efficiency')
        eta_charge = Variable(r'\eta_{charge}', 0.95, '-',
                              'Battery charging efficiency')
        eta_discharge = Variable(r'\eta_{discharge}', 0.95, '-',
                                 'Battery discharging efficiency')
        t_day = Variable('t_{day}', 'hr', 'Daylight span')
        t_night = Variable('t_{night}', 16, 'hr', 'Night span')
        th = Variable(r'\theta', 0.35, '-', 
                      'cosine of incidence angle at 45 lat')

        constraints = [PS_irr*eta_solar*S >= P_oper + P_charge,
                       P_oper >= P_shaft + P_acc,
                       P_charge >= E_batt/(t_day*eta_charge*th),
                       t_day + t_night <= 24*units.hr,
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
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K", "air specific heating value")
        TH = (g*M_atm/R_atm/L_atm).value

        constraints = [
            #h <= 20000*units.m,  # Model valid to top of troposphere
            T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
            # http://en.wikipedia.org/wiki/Density_of_air#Altitude
            rho <= p_sl*T_atm**(TH-1)*M_atm/R_atm/(T_sl**TH)]
        Model.__init__(self, None, constraints, **kwargs)

class StationKeeping(Model):
    def __init__(self, **kwargs):
        
        h_min = Variable('h_{min}', 15000, 'ft', 'minimum altitude')
        V_wind = Variable('V_{wind}', 10, 'm/s', 'wind speed')
        
        constraints = [h >= h_min,
                       V >= V_wind]
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



        constraints = []
        lc = LinkedConstraintSet([self.submodels, constraints])

        objective = W

        Model.__init__(self, objective, lc, **kwargs)

if __name__ == "__main__":
    M = SolarHALE()
    sol = M.solve("mosek")

    if PLOT:

        # number of data points
        N = 30

        M.substitutions.update({ 'h_{batt}': ('sweep', [250,350,500])})

        # plotting at 45 degree latitude

        M.substitutions.update({'V_{wind}': ('sweep', np.linspace(5,40,N))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
        
        W = sol('W')
        b = sol('b')
        V_wind = sol('V_{wind}')
        ind = np.nonzero(V_wind.magnitude==5.)
        ind = ind[0]

        plt.close()
        plt.rcParams['lines.linewidth']=2
        plt.rcParams['font.size']=13
        plt.rcParams
        curve1, = plt.plot(V_wind[ind[0]:ind[1]], b[ind[0]:ind[1]], 
                           label='batt: 250 [Whr/kg]') 
        curve2, = plt.plot(V_wind[ind[1]:ind[2]], b[ind[1]:ind[2]], 
                           label = 'batt: 350 [Whr/kg]')
        curve3, = plt.plot(V_wind[ind[2]:V_wind.size], b[ind[2]:V_wind.size], 
                           label='batt: 500 [Whr/kg]')
        curve4, = plt.plot([25,25], [0,200], '--', color='r', label='90% wind speed')
        plt.ylabel('wing span [ft]')
        plt.xlabel('wind speed on station [m/s]')
        plt.legend(handles = [curve1, curve2, curve3, curve4], loc =4), 
        plt.grid()
        plt.axis([5,40,0,200])
        plt.savefig('bvsV_wind.pdf')
        
        M.substitutions.update({'V_{wind}': 10})
        M.substitutions.update({'h': ('sweep', np.linspace(15000,50000,N))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
        
        W = sol('W')
        b = sol('b')
        h = sol('h')
        ind = np.nonzero(h.magnitude==15000.)
        ind = ind[0]

        plt.close()
        plt.rcParams['lines.linewidth']=2
        plt.rcParams['font.size']=13
        plt.rcParams
        curve1, = plt.plot(h[ind[0]:ind[1]], b[ind[0]:ind[1]], 
                           label='batt: 250 [Whr/kg]') 
        curve2, = plt.plot(h[ind[1]:ind[2]], b[ind[1]:ind[2]], 
                           label = 'batt: 350 [Whr/kg]')
        curve3, = plt.plot(h[ind[2]:h.size],  b[ind[2]:h.size], 
                           label='batt: 500 [Whr/kg]')
        plt.ylabel('wing span [ft]')
        plt.xlabel('wind speed on station [m/s]')
        plt.legend(handles = [curve1, curve2, curve3], loc =2), 
        plt.grid()
        plt.axis([15000,50000,0,200])
        plt.savefig('bvsh.pdf')
        
        # plotting at 30 degree latitude
        M.substitutions.update({r'\theta': 0.6, 't_{night}': 14})

        M.substitutions.update({'h':15000})
        M.substitutions.update({'V_{wind}': ('sweep', np.linspace(5,40,N))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
        
        W = sol('W')
        b = sol('b')
        V_wind = sol('V_{wind}')
        ind = np.nonzero(V_wind.magnitude==5.)
        ind = ind[0]

        plt.close()
        plt.rcParams['lines.linewidth']=2
        plt.rcParams['font.size']=13
        plt.rcParams
        curve1, = plt.plot(V_wind[ind[0]:ind[1]], b[ind[0]:ind[1]], 
                           label='batt: 250 [Whr/kg]') 
        curve2, = plt.plot(V_wind[ind[1]:ind[2]], b[ind[1]:ind[2]], 
                           label = 'batt: 350 [Whr/kg]')
        curve3, = plt.plot(V_wind[ind[2]:V_wind.size],  b[ind[2]:V_wind.size], 
                           label='batt: 500 [Whr/kg]')
        curve4, = plt.plot([25,25], [0,200], '--', color='r', 
                            label='90% wind speed')
        plt.ylabel('wing span [ft]')
        plt.xlabel('wind speed on station [m/s]')
        plt.legend(handles = [curve1, curve2, curve3, curve4], loc =2), 
        plt.grid()
        plt.axis([5,40,0,200])
        plt.savefig('bvsV_wind30.pdf')
        
        M.substitutions.update({'V_{wind}': 10})
        M.substitutions.update({'h': ('sweep', np.linspace(15000,50000,N))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
        
        W = sol('W')
        b = sol('b')
        h = sol('h')
        ind = np.nonzero(h.magnitude==15000.)
        ind = ind[0]

        plt.close()
        plt.rcParams['lines.linewidth']=2
        plt.rcParams['font.size']=13
        plt.rcParams
        curve1, = plt.plot(h[ind[0]:ind[1]], b[ind[0]:ind[1]], 
                           label='batt: 250 [Whr/kg]') 
        curve2, = plt.plot(h[ind[1]:ind[2]], b[ind[1]:ind[2]], 
                           label = 'batt: 350 [Whr/kg]')
        curve3, = plt.plot(h[ind[2]:h.size],  b[ind[2]:h.size], 
                           label='batt: 500 [Whr/kg]')
        plt.ylabel('wing span [ft]')
        plt.xlabel('wind speed on station [m/s]')
        plt.legend(handles = [curve1, curve2, curve3], loc =2), 
        plt.grid()
        plt.axis([15000,50000,0,200])
        plt.savefig('bvsh30.pdf')
```

# Results

Using the above assumption the solar aircraft was sized for three different battery types and at two different latitudes, 45 degrees and 30 degrees latitude. Table 2 shows the different battery efficienies that were used. 

\begin{table}[H]
\begin{center}
\label{t:battery}
\begin{tabular}{ | c |  c | }
    \hline
    Battery Efficiency & Battery Type \\ \hline
    250 $[Whr/kg]$ & Best known lithium ion battery \\\hline
    350 $[Whr/kg]$ & Best known lithium sulphur battery \\
    & (used on the Zephyr) \\\hline
    500 $[Whr/kg]$ & Theoretical limit \\
    & (only achieved in labs) \\\hline
\end{tabular}
\caption{Table showing the best battery efficiencies for lithium ion, lithium sulphur and the theoretical limit.}
\end{center}
\end{table}

## Size vs Wind Speed

As explained in figure 1, wind speeds are important for this analysis because the wind affects the size of the engine and therefore battery needed to sustain flight. As shown in figure 2, in order to reach the 90th percentile winds (i.e. 90% availability), the wing span must be in excess of 200 ft.  Please note that this result was done assuming that we are flying at 45 degrees latitude and during the winter solstice at 15,000 ft. These conditions will guarantee year-round availability.  Figure 3 is shows the same plot but for 30 degree latitude. (Each point on the graph is the minimum weight aircraft that will fly at that latitude during the winter solstice.)


\begin{figure}[H]
	\label{f:bvsV_wind}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsV_wind}
	\caption{This figure shows how at the 45 degree latitude, at the winter solstice, in order to overcome the 90\% windspeeds at 15,000 ft the aircraft needs to be larger than is structurally feasible.}
	\end{center}
\end{figure}

\begin{figure}[H]
	\label{f:bvsV_wind30}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsV_wind30}
	\caption{This figure shows how at the 30 degree latitude, at the winter solstice, in order to overcome the 90\% wind speeds at 15,000 ft the aircraft needs to be larger than is structurally feasible.}
	\end{center}
\end{figure}

## Size vs Altitude

Reaching an altitude of near 50,000 ft is desirable to: 1) reach lower winds, and 2) escape heavy jet stream.  The difficulty with going higher is mainly that the air is thinner and requires faster speeds to produce sufficient lift.  Faster speeds mean larger engines, which means larger batteries and solar cell area.  Figure 4 shows the trade off in size, shown in wing span, vs altitude.  To reach an altitude of 50,000 ft or higher the plane must become larger than is feasibly buildable. Figure 5 shows the same trade off but for the 30 degree latitude case. Please note that this study assumes that wind is not a constraint.  So even if it were possible to reach higher altitudes, the higher winds at that altitude would still make this design infeasible. 

\begin{figure}[H]
	\label{f:bvsh}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsh}
	\caption{This figure shows how at the 45 degree latitude, during the winter solstice, in order to reach altitudes of 50,000 ft the plane must be larger than is structurally feasible. (Note: this study assumes that wind is not a constraint.)}
	\end{center}
\end{figure}

\begin{figure}[H]
	\label{f:bvsh30}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsh30}
	\caption{This figure shows how at the 30 degree latitude, during the winter solstice, the wing span size is reasonable.  However, this assumes that wind is not a constraint, therefore further investigation needs to be done to show how wind speeds would affect the design at the 30 degree latitude case.}
	\end{center}
\end{figure}

# Conclusion

Based on the sizing results done in this study we believe that a solar aircraft cannot achieve the threshold requirements outlined by the customer.  Flying at the 45 parallel on the winter solstice gives insufficient sunlight at low angles to the solar array. Combined with the need to fly in windy conditions, the aircraft size grows to unbuildable dimensions. If the customer is willing to sacrifice availability by flying in slower winds at lower latitudes and limited periods of the year, a solar option could be further explored.   
