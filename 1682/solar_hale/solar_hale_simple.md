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
from solar_hale import SolarHALE
import matplotlib.pyplot as plt
from plotting import poor_mans_contour
import numpy as np

PLOT = True

if __name__ == "__main__":
    M = SolarHALE()
    sol = M.solve("mosek")

    YLIM_b = [0, 200]

    if PLOT:
    
        fig, ax = poor_mans_contour(M, "V_{wind}", np.linspace(5, 40, 30), "h_{batt}", [250,350, 500], "b", YLIM_b, vref=25, vrefname="90% wind speed")
        fig.savefig('bvsV_wind4515.pdf')
        
        M.substitutions.update({'h': 50000})
        fig, ax = poor_mans_contour(M, "V_{wind}", np.linspace(5, 40, 30), "h_{batt}", [250,350, 500], "b", [0, 200], vref=25, vrefname="90% wind speed")
        fig.savefig('bvsV_wind4550.pdf')
        
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

