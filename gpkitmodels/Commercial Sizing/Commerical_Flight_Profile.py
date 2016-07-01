"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units, LinkedConstraintSet, SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
from atmosphere import Atmosphere

"""
minimizes the aircraft total weight, must specify all weights except fuel weight, so in effect
we are minimizing the fuel weight
Rate of climb equation taken from John Anderson's Aircraft Performance and Design (eqn 5.85)
"""

#TODO
#link in with atmosphere
#link in with engine

#-----------------------------------------------------------
#defining the number of segments
Ntakeoff = 0
Nclimb1 = 2
Ntrans = 0
Nclimb2 = 0
Ncruise1 = 0
Ncruiseclimb = 0
Ncruise2 = 0
Ndecent = 0
Nlanding = 0

#segments below here that start with res --> fuel reserves
Nresclimb = 0
Nresdecent = 0
Nreslanding = 0
Nreshold = 0

#total number of flight segments in each category, possibly useful
Nclimb = Nclimb1 + Nclimb2 + Ncruiseclimb
Ncruise = Ncruise1 + Ncruise2
Nres = Nresclimb + Nresdecent + Nreslanding + Nreshold

Nseg = Nclimb + Ncruise + Nres + Ntakeoff + Ntrans + Ndecent + Nlanding

#create index lists for the vector variables

#partial flight profile for use during development
iclimb1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))
if Nclimb2 != 0:
    iclimb2 = map(int, np.linspace(iclimb1[len(iclimb1)-1] + 1, Nclimb2 + iclimb1[len(iclimb1)-1], Nclimb2))
else:
    iclimb2 =0
    
if Ncruise2 != 0:
    icruise2 = map(int, np.linspace(iclimb2[len(iclimb2)-1] + 1, Ncruise2 + iclimb2[len(iclimb2)-1], Ncruise2))
else:
    icruise2 = 0

izbre = map(int, np.linspace(0, Ncruise2 - 1, Ncruise2))

#---------------------------------------------------------------------
#Variable definitions that apply to multiple classes
#physical constants
g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

#air properties
a = VectorVariable(Nseg, 'a', 'm/s', 'Speed of Sound')
rho = VectorVariable(Nseg, '\rho', 'kg/m^3', 'Air Density')
p = VectorVariable(Nseg, 'p', 'Pa', 'Pressure')
mu = VectorVariable(Nseg, '\mu', 'kg/m/s', 'Air Kinematic Viscosity')
T = VectorVariable(Nseg, 'T', 'K', 'Air Temperature')

#altitude
hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude [feet]')
dhft = VectorVariable(Nclimb, 'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
h = VectorVariable(Nseg, 'h', 'm', 'Altitude [meters]')
htoc = Variable('h_{toc}', 'ft', 'Altitude at Top of Climb')
dhClimb1 = Variable('dh_{climb1}', 8500, 'feet', 'Total Altitude Change Required in Climb 1')

#time
tmin = VectorVariable(Nseg, 'tmin', 'min', 'Flight Time in Minutes')
thours = VectorVariable(Nseg, 'thr', 'hour', 'Flight Time in Hours')

#Range  
RngClimb = VectorVariable(Nclimb, 'RngClimb', 'miles', 'Segment Range During Climb')
if Ncruise !=0:
    RngCruise = VectorVariable(Ncruise2 + Ncruise1, 'RngCCruise', 'miles', 'Segment Range During Cruise')
ReqRngCruise = Variable('ReqRngCruise', 'miles', 'Required Cruise Range')
ReqRng = Variable('ReqRng', 'miles', 'Required Mission Range')

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
LDmax = Variable('\\frac{L}{D}_{max}', '-', 'Maximum Lift to Drag')

Cd0 = Variable('C_{d_0}', '-', 'Aircraft Cd0')
K = Variable('K', '-', 'K for Parametric Drag Model')

#aircraft geometry

S = Variable('S', 'm^2', 'Wing Planform Area')

#velocitites and mach numbers
V = VectorVariable(Nseg, 'V', 'knots', 'Aircraft Flight Speed')
M = VectorVariable(Nseg, 'M', '-', 'Aircraft Mach Number')
Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')

#Climb/Decent Rate
RC = VectorVariable(Nseg, 'RC', 'feet/min', 'Rate of Climb/Decent')
theta = VectorVariable(Nseg, '\\theta', '-', 'Aircraft Climb Angle')

#breguet parameter
z_bre = VectorVariable(Nseg, 'z_{bre}', '-', 'Breguet Parameter')

#engine
TSFC = VectorVariable(Nseg, 'TSFC', 'lb/hr/lbf', 'Thrust Specific Fuel Consumption')

#currently sets the value of TSFC, just a place holder
c1 = Variable('c1', 'lb/lbf/hr', 'Constant')

#thrust
thrust = Variable('thrust', 'N', 'Engine Thrust')

#temporary args
#pass i na 1 for testing climb segment 1

class CommericalMissionConstraints(Model):
    """
    class that is general constraints that apply across the mission
    """
    def __init__(self, test=0, **kwargs): 
        constraints = []
        
        #define variables local to this class
        alt10k = Variable('alt10k', 10000, 'feet', 'Altitude where 250kt Speed Limit Stops')
        alt1k = Variable('alt1k', 1500, 'feet', 'Altitude where Climb Profile Starts')
        
        constraints.extend([
            #speed of sound
            a  == (gamma * R * T)**.5,
            
            #convert m to ft
##            hft  == h,
            
            #convert min to hours
            tmin  == thours ,

##            W_fuel[iclimb1] == 1000*units('lbf'),
            
            #constraints on the various weights
            TCS([W_e + W_payload + W_ftotal <= W_total]),
            W_start[0]  == W_total,
            TCS([W_e + W_payload <= W_end[Nseg-1]]),
            TCS([W_ftotal >= sum(W_fuel)]),
            ])
        
        with gpkit.SignomialsEnabled():
            if test != 1:
                constraints.extend([
                    #range constraints
                    TCS([sum(RngClimb) + sum(RngCruise) >= ReqRng]),
                    TCS([ReqRngCruise   >= sum(RngCruise)]),
                    #altitude matching constraints
                    hft[icruise2]==hft[Nclimb-1],
                    ])
            if test ==1:
                 constraints.extend([
                    sum(RngClimb)>= ReqRng
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
        Model.__init__(self, W_ftotal, constraints, **kwargs)
        
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
        excessP = VectorVariable(Nclimb, 'Excess Power', 'W', 'Excess Power During Climb')
        D = VectorVariable(Nclimb, 'Drag', 'N', 'Drag')
        RCtest = VectorVariable(Nseg, 'RCtest', 'feet/min', 'Rate of Climb/Decent TEST')
        constraints = []

        constraints.extend([            
            #set the velocity limits
            V[iclimb1] <= climbspeed,
            V[iclimb1] >= Vstall,
            
            #climb rate constraints
##            RC[iclimb1] == 14*units('ft/min'),
            
            TCS([excessP[iclimb1]+V[iclimb1]*D <= V[iclimb1]*thrust]),

            TCS([D >= (.5*S*rho[iclimb1]*V[iclimb1]**2)*(Cd0 + K*(W_start[iclimb1]/(.5*S*rho[iclimb1]*V[iclimb1]**2))**2)]),

            RC[iclimb1] == excessP[iclimb1]/W_start[iclimb1],
            
            #make the small angle approximation and compute theta
            theta[iclimb1]*V[iclimb1]  == RC[iclimb1],
           
            dhft[iclimb1]  == tmin[iclimb1] * RC[iclimb1],
##            #compute the distance traveled for each segment

##            #takes into account two terms of a cosine expansion
            TCS([RngClimb[iclimb1] + .5*thours[iclimb1]*V[iclimb1]*theta[iclimb1]**2 <= thours[iclimb1]*V[iclimb1]]),
##            TCS([RngClimb[iclimb1]  == thours[iclimb1]*V[iclimb1]]),
            
##            #RngClimb[iclimb1]  == 150*units('miles'),
            
##            #compute fuel burn from TSFC
            W_fuel[iclimb1]  == g * TSFC[iclimb1] * thours[iclimb1] * thrust,

##            #compute the dh required for each climb 1 segment
            dhft[iclimb1] == dhClimb1/Nclimb1,

            #subsitute later
            TSFC[iclimb1]  == c1,
            rho[iclimb1] == 1.225*units('kg/m^3'),
            T[iclimb1] == 273 * units('K')
            ])
            
        Model.__init__(self, None, constraints, **kwargs)
        
#--------------------------------------
#transition between climb 1 and 2 if desired

            
class Climb2(Model):
    """
    class to model the climb portion above 10,000'
    """
    def __init__(self, **kwargs):
        constraints = []

        with gpkit.SignomialsEnabled():
        
            constraints.extend([            
                #set the velocity limits

                #needs to be replaced by an actual Vne and a mach number
                V[iclimb2] <= 500*units('kts'),
                V[iclimb2] >= Vstall,
                
                #climb rate constraints
                TCS([RC[iclimb2] + 0.5 * (V[iclimb2]**3) * rho[iclimb2] * S / W_start[iclimb2] * Cd0 +
                W_start[iclimb2] / S * 2 * K / rho[iclimb2] / V[iclimb2] <= V[iclimb2] * thrust / W_start[iclimb2]]),
                
                #make the small angle approximation and compute theta
                theta[iclimb2]*V[iclimb2]  == RC[iclimb2],
               
                dhft[iclimb2]  == tmin[iclimb2] * RC[iclimb2],
                #compute the distance traveled for each segment

                #takes into account two terms of a cosine expansion
                TCS([RngClimb[iclimb2] + .5*thours[iclimb2]*V[iclimb2]*theta[iclimb2]**2 <= thours[iclimb2]*V[iclimb2]]),

                #compute fuel burn from TSFC
                W_fuel[iclimb2]  == g * TSFC[iclimb2] * thours[iclimb2] * thrust,




                #substitute later
                TSFC[iclimb2]  == c1,
                ])
        Model.__init__(self, None, constraints, **kwargs)
        
#--------------------------------------
#cruise #1
#Breguet Range discretized to model the cruise


#---------------------------------------
#cruise climb/decent...might switch altitude in cruise

class Cruise2(Model):
    """
    class to model the second cruise portion of a flight (if the flight is
    long enough to mandate two cruise portions)
    Model is based off of a discretized Breguet Range equation
    """
    def __init__(self, **kwargs):
        constraints = []

        with gpkit.SignomialsEnabled():
            constraints.extend([
                #constrain the climb rate by holding altitude constant
                hft[icruise2]  == 35000*units('ft'),
                
                #taylor series expansion to get the weight term
                TCS([W_fuel[icruise2]/W_end[icruise2] >= te_exp_minus1(z_bre[izbre], nterm=3)]),
                
                #breguet range eqn
                TCS([RngCruise[izbre] <= z_bre[izbre]*LD[icruise2]*V[icruise2]/(TSFC[icruise2]*g)]),
                
                #time
                thours[icruise2]*V[icruise2]  == RngCruise[izbre],
                ])
            
            #constraint on the aircraft meeting the required range
            for i in range(min(izbre), max(izbre)+1):
                constraints.extend([
                     TCS([RngCruise[i]   >= (i+1) * ReqRngCruise/Nseg])
                    ])
                
            constraints.extend([
                #substitue these values later
                TSFC[icruise2]  == 1.4*units('lb/lbf/hr'),
                LD[icruise2]  == 10,
                V[icruise2]  == 420*units('kts'),
                TSFC[icruise2]  == 1.4*units('lb/lbf/hr'),
                LD[icruise2]  == 10,
                V[icruise2]  == 420*units('kts'),
                ])
        Model.__init__(self, None, constraints, **kwargs)
        
#---------------------------------------
#decent

#----------------------------------------
#landing

class CommercialAircraft(Model):
    """
    class to link all models needed to simulate a commercial flight
    """
    def __init__(self, **kwargs):
        #define all the submodels
        cmc = CommericalMissionConstraints(1)
        climb1 = Climb1()
##        climb2 = Climb2()
##        c2 = Cruise2()
        atm = Atmosphere()

        substitutions = {      
            'W_{e}': 40000*units('lbf'),
            'W_{payload}': 400000*units('lbf'),
            'V_{stall}': 120,
            '\\frac{L}{D}_{max}': 15,
            'ReqRng': 300,
            'C_{d_0}': .025,
            'K': 0.035,
            'S': 124.58,
            'thrust': 3640000,
            'c1': 2,
            }

        self.submodels = [cmc, climb1]

        lc = LinkedConstraintSet([self.submodels])

        Model.__init__(self, cmc.cost, lc, substitutions, **kwargs)
    
if __name__ == '__main__':
    m = CommercialAircraft()
    sol = m.localsolve(kktsolver = "ldl", verbosity = 4)
    
    plt.plot(np.cumsum(sol('tmin')), sol('hft'))
    plt.title('Altitude vs Time')
    plt.ylabel('Altitude [ft]')
    plt.xlabel('Time [min]')
    plt.show()
    
#full flight profile
##        itakeoff = map(int, np.linspace(0, Ntakeoff - 1, Ntakeoff))
##        iclimb1 = map(int, np.linspace(Ntakeoff, itakeoff[len(itakeoff)-1]+Nclimb1, Nclimb1))
##        itrans = map(int, np.linspace(iclimb1[len(iclimb1)-1] + 1, Ntrans + iclimb1[len(iclimb1)-1], Ntrans))
##        iclimb2 = map(int, np.linspace(itrans[len(itrans)-1] + 1, Nclimb2 + itrans[len(itrans)-1], Nclimb2))
##        icruise1 = map(int, np.linspace(iclimb2[len(iclimb2)-1] + 1, Ncruise1 + iclimb2[len(iclimb2)-1], Ncruise1))
##        icruiseclimb = map(int, np.linspace(icruise1[len(icruise1)-1] + 1, Ncruiseclimb + icruise1[len(icruise1)-1], Ncruiseclimb))
##        icruise2 = map(int, np.linspace(icruiseclimb[len(icruiseclimb)-1] + 1, Ncruise2 + icruiseclimb[len(icruiseclimb)-1], Ncruise2))
##        idecent = map(int, np.linspace(icruise2[len(icruise2)-1] + 1, Ndecent + icruise2[len(icruise2)-1], Ndecent))
##        ilanding = map(int, np.linspace(idecent[len(idecent)-1] + 1, Nlanding + idecent[len(idecent)-1], Nlanding))
##        iresclimb = map(int, np.linspace(ilanding[len(ilanding)-1] + 1, Nresclimb + ilanding[len(ilanding)-1], Nresclimb))
##        iresdecent = map(int, np.linspace(iresclimb[len(iresclimb)-1] + 1, Nresdecent + iresclimb[len(iresclimb)-1], Nresdecent))
##        ireslanding = map(int, np.linspace(iresdecent[len(iresdecent)-1] + 1, Nreslanding + iresdecent[len(iresdecent)-1], Nreslanding))
##        ireshold = map(int, np.linspace(ireslanding[len(ireslanding)-1] + 1, Nreshold + ireslanding[len(ireslanding)-1], Nresdecent))
