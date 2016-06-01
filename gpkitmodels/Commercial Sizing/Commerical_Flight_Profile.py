"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
from getatm import get_atmosphere, get_T, get_rho, get_mu
from get_sum import get_hft, get_cruiserng

#TODO

"""
the error is almost certianly due to the fact I haven't set an initial altitude....I guess it should be 1500'
"""

class CommericalAircraft(Model):
    """
    minimizes the aircraft total weight, must specify all weights except fuel weight, so in effect
    we are minimizing the fuel weight
    """
    def __init__(self):
        #define the number of discretizations for each flight segment
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

        #partial flight profile for use in development
        iclimb1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))
        icruise2 = map(int, np.linspace(iclimb1[len(iclimb1)-1] + 1, Ncruise2 + iclimb1[len(iclimb1)-1], Ncruise2))

        izbre = map(int, np.linspace(0, Ncruise2 - 1, Ncruise2))
        #---------------------------------------
        #Variable definitions
        
        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #altitude
        hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude [feet]')
        dhft = VectorVariable(Nclimb1, 'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        h = VectorVariable(Nseg, 'h', 'm', 'Altitude [meters]')

        #air properties
        a = VectorVariable(Nseg, 'a', 'm/s', 'Speed of Sound')
        rho = VectorVariable(Nseg, 'rho', get_rho, 'kg/m^3', 'Air Density', args=[h])
        mu = VectorVariable(Nseg, 'mu', get_mu, 'kg/m/s', 'Air Kinematic Viscosity', args=[h])
        T = VectorVariable(Nseg, 'T', get_T, 'K', 'Air Temperature', args=[h])

        #time
        tmin = VectorVariable(Nseg, 'tmin', 'min', 'Flight Time in Minutes')
        thours = VectorVariable(Nseg, 'thr', 'hour', 'Flight Time in Hours')

        #Range

        #HOW TO HANDLE NM
        
        RngClimb = VectorVariable(Nclimb, 'RngClimb', 'miles', 'Segment Range During Climb')
    #    RngCruise = VectorVariable(Ncruise2 + Ncruise1, 'RngCCruise', 'miles', 'Segment Range During Cruise')
    #    ReqRngCruise = Variable('ReqRngCruise', 'miles', 'Required Cruise Range')
        ReqRng = Variable('ReqRng', 30, 'miles', 'Required Mission Range')

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

        #CHECK VALUE OF CD0, k, AND K UNITS
        
        Cd0 = Variable('C_{d_0}', .025, '-', 'Aircraft Cd0')
        K = Variable('K', .9, '-', 'K for Parametric Drag Model')

        #aircraft geometry

        #CHECK THE VALUE OF S
        
        S = Variable('S', 124.58, 'm^2', 'Wing Planform Area')

        #velocitites and mach numbers
        V = VectorVariable(Nseg, 'V', 'knots', 'Aircraft Flight Speed')
        M = VectorVariable(Nseg, 'M', '-', 'Aircraft Mach Number')
        Vstall = Variable('V_{stall}', 120, 'knots', 'Aircraft Stall Speed')

        #Climb/Decent Rate
        RC = VectorVariable(Nseg, 'RC', 'feet/min', 'Rate of Climb/Decent')
        theta = VectorVariable(Nseg, '\\theta', '-', 'Aircraft Climb Angle')

        #breguet parameter
        z_bre = VectorVariable(Nseg, 'z_{bre}', '-', 'Breguet Parameter')

        #engine
        TSFC = VectorVariable(Nseg, 'TSFC', 'lbm/hr/lbf', 'Thrust Specific Fuel Consumption')

        s = Variable('s', '-', 'Feasability Analysis')
        s2 = Variable('s2', '-', 'Feasability Analysis')

        with gpkit.SignomialsEnabled():
        
            constraints = []

            #---------------------------------------
            #basic constraint definitions, apply to the entire flight
            alt10k = Variable('alt10k', 10000, 'feet', 'Altitude where 250kt Speed Limit Stops')
            alt1k = Variable('alt1k', 1500, 'feet', 'Altitude where Climb Profile Starts')
            cruiseRng = Variable('cruiseRng', get_cruiserng, 'miles', 'Crusie Range Needed', args=[ReqRng, RngClimb])
            
            constraints.extend([
                #speed of sound
                a  == (gamma * R * T)**.5,
                T>=1*units('K'),
                #convert m to ft
                hft  == h,
                #convert min to hours
                tmin  == thours ,
                #constraints on the various weights
                W_e + W_payload + W_ftotal <= W_total,
                W_start[0]  == W_total,
                W_e + W_payload <= W_end[Nseg-1],
                W_ftotal   >= sum(W_fuel),
                #altitude at end of climb segment 1...constraint comes from 250kt speed limit below 10,000'
                hft[Ntakeoff + Nclimb1 - 1] <= alt10k,
                hft[0] == alt1k,
                
                #range constraint
         #       ReqRngCruise   == cruiseRng,

                #substitute these values later
                W_e  == 90000*units('lbf'),
                W_payload  == 50000*units('lbf'),

                W_ftotal>=W_e
                ])

            #constrain the segment weights in a loop
            for i in range(1, Nseg):
                constraints.extend([
                    W_start[i]   == W_end[i-1] 
                    ])
            for i in range(0,Nseg):
                constraints.extend([
                    W_start[i]   >= W_end[i] + W_fuel[i]
                    ])
            
            #---------------------------------------
            #takeoff

            #--------------------------------------
            #Climb #1 (sub 10K subject to 250KTS speed limit)
            climbspeed = Variable('climbspeed', 250, 'kts', 'Speed Limit Under 10,000 ft')
            thrust = Variable('thrust', 49140, 'lbf', 'Engine Thrust')
           
            #cos_te = VectorVariable(Nclimb1, 'cos_te', get_te, '-', 'Cosine Taylor Expansion', args=[theta])
            
            constraints.extend([            
                #set the velocity limits
                V[iclimb1] <= climbspeed,
                V[iclimb1]   >= Vstall,
                #climb rate constraints
                RC[iclimb1] + 0.5 * (V[iclimb1]**3) * rho[iclimb1] * S / W_start[iclimb1] * Cd0 +
                W_start[iclimb1] / S * 2 * K / rho[iclimb1] / V[iclimb1] <= V[iclimb1] * thrust / W_start[iclimb1],
                #make the small angle approximation and compute theta
                theta[iclimb1]*V[iclimb1]  == RC[iclimb1],
               
                dhft[iclimb1]  == tmin[iclimb1] * RC[iclimb1],
                #compute the distance traveled for each segment

                #ASSUMES SMALL ANGLES HERE AND IN THETA COMPUTATION, FIX THIS LATER --------------------------------------------------------

        #        hft==8000,
                
                RngClimb[iclimb1]  == thours[iclimb1]*V[iclimb1],

                RngClimb[Nseg-1] >= ReqRng,

                #compute fuel burn from TSFC
                W_fuel[iclimb1]  == g * TSFC[iclimb1] * thours[iclimb1] * thrust,


                #sub these later
                TSFC[iclimb1]  == 1.6*units('lbm/lbf/hr')
                ])
            
            for i in range(0, Nclimb1):
                if i==0:
                    constraints.extend([
                        hft[i]==1500*units('ft')
                        ])
                else:
                     constraints.extend([
                        hft[i] <= hft[i-1]+dhft[i]
                        ])
            
            #--------------------------------------
            #transition between climb 1 and 2 if desired

            #--------------------------------------
            #Climb #2 (over 10K, no speed limit)

            #--------------------------------------
            #cruise #1
            #Breguet Range discretized to model the cruise
            

            #---------------------------------------
            #cruise climb/decent...might switch altitude in cruise

            #---------------------------------------
            #cruise #2
            #Breguet Range discretized to model the cruise

            """constraints.extend([
                #constrain the climb rate by holding altitude constant
                hft[icruise2]  == 35000  ,
                #taylor series expansion to get the weight term
                W_fuel[icruise2]/W_end[icruise2]   >=  te_exp_minus1(z_bre[izbre], nterm=3),
                #breguet range eqn
                RngCruise[izbre] <= z_bre[izbre]*LD[icruise2]*V[icruise2]/(TSFC[icruise2]*g)  ,
                #time
                thours[icruise2]*V[icruise2]  == RngCruise[izbre]  ,


                #substitue these values later
                TSFC[icruise2]  == 1.4*units('lbm/lbf/hr'),
                LD[icruise2]  == 10,
                V[icruise2]  == 420*units('kts')
                ])
            
            #constraint on the aircraft meeting the required range
            for i in range(min(izbre), max(izbre)+1):
                constraints.extend([
                     RngCruise[i]   >= (i+1) * ReqRngCruise/Nseg
                    ])"""
                   
            #---------------------------------------
            #decent

            #----------------------------------------
            #landing


        #objective is to minimize the total required fuel weight, accomplished my minimizing the total weight
        objective =  W_total

        m = Model(objective, constraints)

        sol = m.localsolve()

if __name__ == '__main__':
    CommericalAircraft()
