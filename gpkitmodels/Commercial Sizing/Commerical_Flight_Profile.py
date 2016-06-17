"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units, LinkedConstraintSet
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
from atmosphere import Troposphere, Tropopause

#TODO
#integrate atmosphere model
#integrate engine model

GAS_CONSTANT = 287  # [J/(kg*K)]
GRAVITATIONAL_ACCEL = 9.81 # [m/s^2]

#define the number of discretizations for each flight segment

Ntakeoff = 0
Nclimb1 = 3
Ntrans = 0
Nclimb2 = 3
Ncruise1 = 0
Ncruiseclimb = 0
Ncruise2 = 3
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

#partial flight profile for use during development
iclimb1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))
iclimb2 = map(int, np.linspace(iclimb1[len(iclimb1)-1] + 1, Nclimb2 + iclimb1[len(iclimb1)-1], Nclimb2))
icruise2 = map(int, np.linspace(iclimb2[len(iclimb2)-1] + 1, Ncruise2 + iclimb2[len(iclimb2)-1], Ncruise2))

izbre = map(int, np.linspace(0, Ncruise2 - 1, Ncruise2))

class FlightProfile(Model):
    """
    minimizes the aircraft total weight, must specify all weights except fuel weight, so in effect
    we are minimizing the fuel weight

    Rate of climb equation taken from John Anderson's Aircraft Performance and Design (eqn 5.85)
    """
    def __init__(self, **kwargs):
        #---------------------------------------
        #Variable definitions
        
        #physical constants
        g   = Variable('g', GRAVITATIONAL_ACCEL, 'm/s^2', 'Gravitational acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R   = Variable('R', GAS_CONSTANT, 'J/(kg*K)', 'Specific gas constant (air)')

        #altitude
        hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude [feet]')
        dhft = VectorVariable(Nclimb, 'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        h = VectorVariable(Nseg, 'h', 'm', 'Altitude [meters]')

        #air properties
        a = VectorVariable(Nseg, 'a', 'm/s', 'Speed of Sound')
        mu  = VectorVariable(Nseg, '\\mu', 'kg/(m*s)', 'Dynamic viscosity')
        p   = VectorVariable(Nseg, 'p', 'Pa', 'Pressure')
        rho = VectorVariable(Nseg, '\\rho', 'kg/m^3', 'Density')
        T   = VectorVariable(Nseg, 'T', 'K', 'Temperature')

        #time
        tmin = VectorVariable(Nseg, 'tmin', 'min', 'Flight Time in Minutes')
        thours = VectorVariable(Nseg, 'thr', 'hour', 'Flight Time in Hours')

        #Range  
        RngClimb = VectorVariable(Nclimb, 'RngClimb', 'miles', 'Segment Range During Climb')
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

        #CHECK VALUE OF CD0, K
        
        Cd0 = Variable('C_{d_0}', '-', 'Aircraft Cd0')
        K = Variable('K', '-', 'K for Parametric Drag Model')

        #aircraft geometry

        #CHECK THE VALUE OF S
        
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

        with gpkit.SignomialsEnabled():
        
            constraints = []

            #---------------------------------------
            #basic constraint definitions, apply to the entire flight
            alt10k = Variable('alt10k', 10000, 'feet', 'Altitude where 250kt Speed Limit Stops')
            alt1k = Variable('alt1k', 1500, 'feet', 'Altitude where Climb Profile Starts')
            
            constraints.extend([
                #speed of sound
                a  == (gamma * R * T)**.5,
                T == 273*units('K'),
                
                #convert m to ft
                hft == h,
                
                #convert min to hours
                tmin == thours,
                
                #constraints on the various weights
                W_e + W_payload + W_ftotal <= W_total,
                W_start[0]  == W_total,
                W_e + W_payload <= W_end[Nseg-1],
                W_ftotal   >= sum(W_fuel),
                
                #altitude at end of climb segment 1...constraint comes from 250kt speed limit below 10,000'
                hft[Ntakeoff + Nclimb1 - 1] <= alt10k,
                
                #range constraints
                sum(RngClimb) + sum(RngCruise) >= ReqRng,
                ReqRngCruise   <= sum(RngCruise),

                #altitude matching constraints
                hft[icruise2]==hft[Nclimb-1]
                ])

            #constrain the segment weights in a loop
            for i in range(1, Nseg):
                constraints.extend([
                    W_start[i] == W_end[i-1] 
                    ])
            for i in range(0,Nseg):
                constraints.extend([
                    W_start[i] >= W_end[i] + W_fuel[i]
                    ])
            
            #---------------------------------------
            #takeoff

            #--------------------------------------
            #Climb #1 (sub 10K subject to 250KTS speed limit)
            climbspeed = Variable('climbspeed', 250, 'kts', 'Speed Limit Under 10,000 ft')
            thrust = Variable('thrust', 40000, 'lbf', 'Engine Thrust')
            
            constraints.extend([            
                #set the velocity limits
                V[iclimb1] <= climbspeed,
                V[iclimb1] >= Vstall,
                
                #climb rate constraints
                RC[iclimb1] + 0.5 * (V[iclimb1]**3) * rho[iclimb1] * S / W_start[iclimb1] * Cd0 +
                W_start[iclimb1] / S * 2 * K / rho[iclimb1] / V[iclimb1] >= V[iclimb1] * thrust / W_start[iclimb1],
                
                #make the small angle approximation and compute theta
                theta[iclimb1]*V[iclimb1]  == RC[iclimb1],
               
                dhft[iclimb1]  == tmin[iclimb1] * RC[iclimb1],
                
                #compute the distance traveled for each segment
                #takes into account two terms of a cosine expansion
                RngClimb[iclimb1] + .5*thours[iclimb1]*V[iclimb1]*theta[iclimb1]**2 <= thours[iclimb1]*V[iclimb1],
                
                #compute fuel burn from TSFC
                W_fuel[iclimb1]  == g * TSFC[iclimb1] * thours[iclimb1] * thrust
                ])
            
            for i in range(0, Nclimb):
                if i==0:
                    constraints.extend([
                        hft[i] <= 1500*units('ft')+dhft[i]
                        ])
                else:
                     constraints.extend([
                        hft[i] <= hft[i-1]+dhft[i]
                        ])
            
            #--------------------------------------
            #transition between climb 1 and 2 if desired

            #--------------------------------------
            #Climb #2 (over 10K, no speed limit)
            constraints.extend([            
                #set the velocity limits

                #needs to be replaced by an actual Vne and a mach number
                V[iclimb2] <= 500*units('kts'),
                V[iclimb2] >= Vstall,
                
                #climb rate constraints
                RC[iclimb2] + 0.5 * (V[iclimb2]**3) * rho[iclimb2] * S / W_start[iclimb2] * Cd0 +
                W_start[iclimb2] / S * 2 * K / rho[iclimb2] / V[iclimb2] >= V[iclimb2] * thrust / W_start[iclimb2],
                
                #make the small angle approximation and compute theta
                theta[iclimb2]*V[iclimb2]  == RC[iclimb2],
               
                dhft[iclimb2]  == tmin[iclimb2] * RC[iclimb2],
                #compute the distance traveled for each segment

                #takes into account two terms of a cosine expansion
                RngClimb[iclimb2] + .5*thours[iclimb2]*V[iclimb2]*theta[iclimb2]**2 <= thours[iclimb2]*V[iclimb2],

                #compute fuel burn from TSFC
                W_fuel[iclimb2]  == g * TSFC[iclimb2] * thours[iclimb2] * thrust
                ])
                     
            #--------------------------------------
            #cruise #1
            #Breguet Range discretized to model the cruise
            

            #---------------------------------------
            #cruise climb/decent...might switch altitude in cruise

            #---------------------------------------
            #cruise #2
            #Breguet Range discretized to model the cruise

            constraints.extend([
                #constrain the climb rate by holding altitude constant
                hft[icruise2]  == 11000*units('m'),
                
                #taylor series expansion to get the weight term
                W_fuel[icruise2]/W_end[icruise2] >= te_exp_minus1(z_bre[izbre], nterm=3),
                
                #breguet range eqn
                RngCruise[izbre] <= z_bre[izbre]*LD[icruise2]*V[icruise2]/(TSFC[icruise2]*g),
                
                #time
                thours[icruise2]*V[icruise2]  == RngCruise[izbre]
                ])
            
            #constraint on the aircraft meeting the required range
            for i in range(min(izbre), max(izbre)+1):
                constraints.extend([
                     RngCruise[i]   >= (i+1) * ReqRngCruise/Nseg
                    ])


#crude replacement for substitions
            constraints.extend([
                #substitue these values later
                TSFC[icruise2]  == 1.4*units('lb/lbf/hr'),
                TSFC[iclimb1] == c1,
                LD[icruise2]  == 10,
                V[icruise2]  == 420*units('kts'),
                TSFC[iclimb2]  == c1,
                W_e  == 90000*units('lbf'),
                W_payload  == 50000*units('lbf'),
            #    rho == 1.225*units('kg/m^3')
                ])
                   
            #---------------------------------------
            #decent

            #----------------------------------------
            #landing


        #objective is to minimize the total required fuel weight, accomplished my minimizing the total weight
        objective =  W_total

        Model.__init__(self, objective, constraints, **kwargs)

class CommercialFlight(Model):
    """
    class to link all the models necessary to simualte an actual commercial flight
    """

    def __init__(self, **kwargs):
        #create all the submodels
        
        fltprofile = FlightProfile()
        atm = Troposphere()

        self.submodels = [fltprofile, atm]

        substitutions = {
            'S': 124.58,
            '\\frac{L}{D}_{max}': 15,
            'K': 0.9,
            'C_{d_0}': .025,
            'ReqRng': 1500,
            'V_{stall}': 120,
            'c1': 30000,
            }

        lc = LinkedConstraintSet([self.submodels])

        
        Model.__init__(self, fltprofile.cost, lc, substitutions)
        
        

if __name__ == "__main__":
        comflt = CommercialFlight()

        sol = comflt.localsolve(verbosity=4)

        plt.plot(np.cumsum(sol('tmin')), sol('hft'))
        plt.title('Altitude vs Time')
        plt.ylabel('Altitude [ft]')
        plt.xlabel('Time [min]')
        plt.show()

