# SOLAR HIGH ALTITUDE LONG ENDURANCE AIRCRAFT

This paper is a presentation of the feasibility of a solar powered aircraft for the given requirements presented in table 1. 

\begin{table}[H]
\begin{center}
\label{t:battery}
\begin{tabular}{ | c |  c | }
    \hline
    \textbf{Requirement} & \textbf{Specification} \\ \hline
    Payload & 4-10 [lbs] 10 [Watts] \\\hline
    Endurance & 5-35 days \\\hline
    Availability & 95\% \\\hline
    Latitude Coverage & 45-60 degrees \\\hline
    Diameter Footprint & 100-300 [km] \\\hline
\end{tabular}
\caption{Table outlining the requirements given to us by the customer.}
\end{center}
\end{table}

This study will show that at the 45th degree parallel, the amount of power required to fly fast enough to match the wind speeds is greater than the power available through solar energy.

The study was done using a convex optimization tool called, GPkit.  The optimization model was constructed using various models which describe the physics required to achieve the specified requirements. The specific requirements to which this plane was optimized are found in table 1. 


# Assumptions in Optimization Model

\begin{enumerate}
    \item Steady level flight. 
    \item Standard aerodynamic model using profile drag and induced drag. 
    \item Constant fractional airframe weight. ($W_{airframe} = W_{plane}f_{airframe}$, where $f_{airframe} = 20\%$)  The fractional weight includes the weight of the wing, fuselage, engine (with parts), tail, and any other structural elements. 
    \item The avionics weight was calculated to be 6 lbs and draws 15 watts of power and is based on an analysis of avionics required for flight control, ground communication, and satellite communication.
    \item The solar cell area density and efficiency is based off of the best available solar cell density. ($\rho_{solar-cell} = 1.2$ [kg/m$^2$], $\eta_{solar-cell}$ = 20\%)
    \item Payload weight is 4 lbs and draws 10 watts of power. 
    \item Flight occurs at the worst case scenario, which at 45 degree latitude and during the winter solstice when the least amount of sunlight is available. By constraining the plane to fly during the winter solstice we guarantee that we can fly during any other day of the year.  This assumptions translates to 8 hours of charging time and an angle of incidence of about 66 degrees. (This is due to the 24 degree range off of the horizon.) 
    \item The optimization minimized weight of the aircraft becasue aircraft cost typically scales with weight.
\end{enumerate}

# Station Keeping Requirements

To achieve the minimum footprint of 100km, the aircraft must fly at 15,000 ft. In order meet the 300km footprint diameter requirement, to avoid diminished solar irradiance due to weather, and to avoid air traffic conflicts, the aircraft must fly at 50,000 ft or above.  Wind speed is critical to this design because as wind speed increases, more flight power is needed, and therefore more batteries and solar cells are needed. 

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
        plt.savefig('bvsV_wind4515.pdf')
        
        # plot at 45 degrees and 50000 ft
        M.substitutions.update({'h': 50000})
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
        curve4, = plt.plot([32,32], [0,200], '--', color='r', label='90% wind speed')
        plt.ylabel('wing span [ft]')
        plt.xlabel('wind speed on station [m/s]')
        plt.legend(handles = [curve1, curve2, curve3, curve4], loc =3), 
        plt.grid()
        plt.axis([5,40,0,200])
        plt.savefig('bvsV_wind4550.pdf')
        
        # plotting at 30 degree latitude and 15000 ft
        M.substitutions.update({r'\theta': 0.6, 't_{night}': 14})

        M.substitutions.update({'h':15000})
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
        plt.savefig('bvsV_wind3015.pdf')
        
        # plot at 30 degrees and 50000 ft
        M.substitutions.update({'h': 50000})
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
        curve4, = plt.plot([32,32], [0,200], '--', color='r', label='90% wind speed')
        plt.ylabel('wing span [ft]')
        plt.xlabel('wind speed on station [m/s]')
        plt.legend(handles = [curve1, curve2, curve3, curve4], loc =2), 
        plt.grid()
        plt.axis([5,40,0,200])
        plt.savefig('bvsV_wind3050.pdf')
```

# Results

The trade study was done by observing the size of the aircraft need to fly at the 45 degree latitude at both the 15,000 and 50,000 ft altitude case.  The 45 degree latitude meets the threshold requirement set by the customer.  The 15,000 ft altitude meets the threshold footprint diameter requirement and the 50,000 ft altitude meets the objective footprint diameter requirment.  An additional, easier case was oberserved at 30 degree latitude at both the 15,000 ft and 50,000 ft altitude cases for reference. 

The plane was also sized for each altitude and latitude case for three different state of the art battery types.  These are shown in table 2. 

\begin{table}[H]
\begin{center}
\label{t:battery}
\begin{tabular}{ | c |  c | }
    \hline
    \textbf{Battery Efficiency} & \textbf{Battery Type} \\ \hline
    250 $[Whr/kg]$ & Best known lithium ion battery \\\hline
    350 $[Whr/kg]$ & Best known lithium sulphur battery \\
    & (used on the Zephyr solar aircraft) \\\hline
    500 $[Whr/kg]$ & Highest observed value for lithium sulphur \\
    & (only achieved in labs) \\\hline
\end{tabular}
\caption{Table showing the best battery efficiencies for lithium ion, lithium sulphur and the theoretical limit. http://techportal.eere.energy.gov/technology.do/techID=1141}
\end{center}
\end{table}

## 45 Degree Latitude Sizing 

Figures 1 and 2 show how the aircraft would be sized at the 45 degree latitude case.  The plots show how wing span would grow based on the average wind speed during the mission.  The 90% wind speed is shown for reference.  As is observed at the 45 degree latitude case, in order for the plane to produce enough power to match the wind speeds of the 90% case the wing span would be well over 200 ft.  As the wind speed increases the plane requires more power to station keep and therefore more solar cells and more batteries to harvest and store that power.  Even for the 15,000 ft altitude case, the wind speeds between 15-18 m/s require about a 100 ft wing span to produce the power necessary regardless of the battry type.  Reaching the 90% wind speed limit is simply infeasible.  For the 50,000 ft altitude case the best available lithium ion battery is not even feasible and does not appear on the graph and reaching the 90% wind speeds is even more difficult than the 15,000 ft altitude case. (Please note that this assumes that flight occurs during the winter solstice.)

\begin{figure}[H]
	\label{f:bvsV_wind}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsV_wind4515}
	\caption{This figure shows how at the 45 degree latitude and at 15,000 ft, the wind span of the aircraft grows based on the average wind speed during the mission.}
	\end{center}
\end{figure}

\begin{figure}[H]
	\label{f:bvsV_wind30}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsV_wind4550}
	\caption{This figure shows how at the 45 degree latitude and at 50,000 ft, the wind span of the aircraft grows based on the average wind speed during the mission.}
	\end{center}
\end{figure}

## 30 Degree Latitude Sizing 

A sizing for the 30 degree latitude case, at both 15,000 and 50,000 ft, was done to oberseve how the aircraft would be sized under less constraining conditions.  Figures 3 and 4 show how the wing span grows with average wind speed during the mission.  While, for the lithium sulphur battery case, a 100 ft wing span span would achieve the 90% availability requirement, it can only do so at the 30th degree parallel and under optimistic assumptions.  Again, the 50,000 ft case is not feasible for the 90% wind speed case. 

\begin{figure}[H]
	\label{f:bvsV_wind}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsV_wind3015}
	\caption{This figure shows how at the 30 degree latitude and at 15,000 ft, the wind span of the aircraft grows based on the average wind speed during the mission.}
	\end{center}
\end{figure}

\begin{figure}[H]
	\label{f:bvsV_wind30}
	\begin{center}
	\includegraphics[scale = 0.5]{bvsV_wind3050}
	\caption{This figure shows how at the 30 degree latitude and at 50,000 ft, the wind span of the aircraft grows based on the average wind speed during the mission.}
	\end{center}
\end{figure}

# Conclusion

Based on the sizing results done in this study we believe that a solar aircraft cannot achieve the threshold requirements outlined by the customer.  The follwing reasons outline why we believe a gas powered aircraft to be a superior solution to solar aircraft for this set of requirements. 

\begin{enumerate}
    \item In order to reach the 90\% avaibility, based on the 90\% wind speeds, the wing span, and therefore weight, become too big to feasibly build.  The power needed to maintain the required availablity at the required latitude is much greater than the power able to be stored and collected by batteries and solar cells for a reasonably sized aircraft. 
    \item State of the art batteries and solar cell technology was assumed.  Even if the constraints were relaxed to a feasible design, the cost to build and maintain the aircraft would be very high. 
    \item This is a relatively unproven concept and the given set of requirements makes it even more prone to risk. 
    \item Gas power aircraft are not limited by latitude. 
    \item Gas powered engines are highly developed and efficient and powerful enough to station keep easily for the 95\% availibility requirements. 
    \item Gas powered aircraft are a highly proven concept with relatively low associated risk. 
\end{enumerate}

