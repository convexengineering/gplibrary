"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units, LinkedConstraintSet
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
from atmosphere import Troposphere, Tropopause

"""
minimizes the aircraft total weight, must specify all weights except fuel weight, so in effect
we are minimizing the fuel weight

Rate of climb equation taken from John Anderson's Aircraft Performance and Design (eqn 5.85)
"""

#TODO

#-----------------------------------------------------------
#defining the number of segments
Nclimb1 = 3
Ntakeoff = 0
Nseg = Nclimb1
Nclimb = Nclimb1
#partial flight profile for use during development
iclimb1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))

#---------------------------------------------------------------------
#Variable definitions that apply to multiple classes
#physical constants
g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

#air properties
a = VectorVariable(Nseg, 'a', 'm/s', 'Speed of Sound')
rho = VectorVariable(Nseg, '\rho', 'kg/m^3', 'Air Density')
mu = VectorVariable(Nseg, '\mu', 'kg/m/s', 'Air Kinematic Viscosity')
T = VectorVariable(Nseg, 'T', 'K', 'Air Temperature')

#altitude
hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude [feet]')
dhft = VectorVariable(Nclimb, 'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
h = VectorVariable(Nseg, 'h', 'm', 'Altitude [meters]')

#time
tmin = VectorVariable(Nseg, 'tmin', 'min', 'Flight Time in Minutes')
thours = VectorVariable(Nseg, 'thr', 'hour', 'Flight Time in Hours')

#Range  
RngClimb = VectorVariable(Nclimb, 'RngClimb', 'miles', 'Segment Range During Climb')
ReqRng = Variable('ReqRng', 300, 'miles', 'Required Mission Range')

#aircraft weights
W_e = Variable('W_{e}', 'lbf', 'Empty Weight of Aircraft')
W_start = VectorVariable(Nseg, 'W_{start}', 'lbf', 'Segment Start Weight')
W_fuel = VectorVariable(Nseg, 'W_{fuel}', 'lbf', 'Segment Fuel Weight')
W_ftotal = Variable('W_{f_{total}}', 'lbf', 'Total Fuel Weight')
W_end = VectorVariable(Nseg, 'W_{end}', 'lbf', 'Segment End Weight')
W_payload = Variable('W_{payload}', 'lbf', 'Aircraft Payload Weight')
W_total = Variable('W_{total}', 'lbf', 'Total Aircraft Weight')

#aero
LD = VectorVariable(Nseg, '\\frac{L}{D}', '-', 'Lift to Drag')
LDmax = Variable('\\frac{L}{D}_{max}', 15, '-', 'Maximum Lift to Drag')

Cd0 = Variable('C_{d_0}', .025, '-', 'Aircraft Cd0')
K = Variable('K', .9, '-', 'K for Parametric Drag Model')

#aircraft geometry

S = Variable('S', 124.58, 'm^2', 'Wing Planform Area')

#velocitites and mach numbers
V = VectorVariable(Nseg, 'V', 'knots', 'Aircraft Flight Speed')
M = VectorVariable(Nseg, 'M', '-', 'Aircraft Mach Number')
Vstall = Variable('V_{stall}', 120, 'knots', 'Aircraft Stall Speed')

#Climb/Decent Rate
RC = VectorVariable(Nseg, 'RC', 'feet/min', 'Rate of Climb/Decent')
theta = VectorVariable(Nseg, '\\theta', '-', 'Aircraft Climb Angle')

#engine
TSFC = VectorVariable(Nseg, 'TSFC', 'lb/hr/lbf', 'Thrust Specific Fuel Consumption')

#currently sets the value of TSFC, just a place holder
c1 = Variable('c1', 2, 'lb/lbf/hr', 'Constant')

#thrust
thrust = Variable('thrust', 87000, 'N', 'Engine Thrust')

class CommericalMissionConstraints(Model):
    """
    class that is general constraints that apply across the mission
    """
    def __init__(self, **kwargs):        
        with gpkit.SignomialsEnabled():
        
            constraints = []
            
            #define variables local to this class
            alt10k = Variable('alt10k', 10000, 'feet', 'Altitude where 250kt Speed Limit Stops')
            alt1k = Variable('alt1k', 1500, 'feet', 'Altitude where Climb Profile Starts')
            
            constraints.extend([
                #speed of sound
                a  == (gamma * R * T)**.5,
                
                #convert m to ft
                hft  == h,
                
                #convert min to hours
                tmin  == thours ,
                
                #constraints on the various weights
                TCS([W_e + W_payload + W_ftotal <= W_total]),
                W_start[0]  == W_total,
                TCS([W_e + W_payload <= W_end[Nseg-1]]),
                TCS([W_ftotal >= sum(W_fuel)]),
                
                #altitude at end of climb segment 1...constraint comes from 250kt speed limit below 10,000'
                hft[Ntakeoff + Nclimb1 - 1] == alt10k,
                
                #range constraints
                TCS([sum(RngClimb) >= ReqRng]),
                ])

            #constrain the segment weights in a loop
            for i in range(1, Nseg):
                constraints.extend([
                    W_start[i] == W_end[i-1] 
                    ])
            for i in range(0,Nseg):
                constraints.extend([
                    TCS([W_start[i] >= W_end[i] + W_fuel[i]])
                    ])

            constraints.extend([
                #substitue these values later
                TSFC[iclimb1]  == c1,
                rho[iclimb1] == 1.225*units('kg/m^3'),
                W_e  == 400000*units('N'),
                W_payload == 400000*units('N'),
                T ==273*units('K'),
                ])
        Model.__init__(self, W_total, constraints, **kwargs)
        
#---------------------------------------
#takeoff

class Climb1(Model):
    """
    class to model the climb portion of a flight, applies to all climbs below
    10,000'
    """
    def __init__(self,**kwargs):
        #Climb #1 (sub 10K subject to 250KTS speed limit)
        climbspeed = Variable('climbspeed', 250, 'kts', 'Speed Limit Under 10,000 ft')

        constraints = []

        with gpkit.SignomialsEnabled():
            constraints.extend([            
                #set the velocity limits
                V[iclimb1] <= climbspeed,
                V[iclimb1] >= Vstall,
                
                #climb rate constraints
                TCS([RC[iclimb1] + 0.5 * (V[iclimb1]**3) * rho[iclimb1] * S / W_start[iclimb1] * Cd0 +
                W_start[iclimb1] / S * 2 * K / rho[iclimb1] / V[iclimb1] <= V[iclimb1] * thrust / W_start[iclimb1]]),
                
                #make the small angle approximation and compute theta
                theta[iclimb1]*V[iclimb1]  == RC[iclimb1],
               
                dhft[iclimb1]  == tmin[iclimb1] * RC[iclimb1],
                #compute the distance traveled for each segment

                #takes into account two terms of a cosine expansion
                TCS([RngClimb[iclimb1] + .5*thours[iclimb1]*V[iclimb1]*theta[iclimb1]**2 <= thours[iclimb1]*V[iclimb1]]),
                
                #compute fuel burn from TSFC
                W_fuel[iclimb1]  == g * TSFC[iclimb1] * thours[iclimb1] * thrust,
                ])
            
            for i in range(0, Nclimb):
                if i==0:
                    constraints.extend([
                        TCS([hft[i] <= 1500*units('ft')+dhft[i]])
                        ])
                else:
                     constraints.extend([
                        TCS([hft[i] <= hft[i-1]+dhft[i]])
                        ])
        Model.__init__(self, None, constraints, **kwargs)
        

class CommercialAircraft(Model):
    """
    class to link all models needed to simulate a commercial flight
    """
    def __init__(self, **kwargs):
        #define all the submodels
        cmc = CommericalMissionConstraints()
        climb1 = Climb1()
    

        self.submodels = [cmc, climb1]

        lc = LinkedConstraintSet([self.submodels])

        Model.__init__(self, cmc.cost, lc, **kwargs)
    
if __name__ == '__main__':
    m = CommercialAircraft()
    sol = m.localsolve(kktsolver = "ldl", verbosity = 4)
    
    plt.plot(np.cumsum(sol('tmin')), sol('hft'))
    plt.title('Altitude vs Time')
    plt.ylabel('Altitude [ft]')
    plt.xlabel('Time [min]')
    plt.show()
    
