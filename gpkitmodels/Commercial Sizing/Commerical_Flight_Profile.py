"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import gpkit
import numpy as np
from gpkit import VectorVariable, Variable, Model, units, ConstraintSet, LinkedConstraintSet, SignomialEquality
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt
from atmosphere import Atmosphere
from collections import defaultdict
from gpkit.small_scripts import mag
from engine_atm import EngineOnDesign, EngineOffDesign, EngineOffDesign2, EngineOffDesign3, EngineOffDesign4, EngineOffDesign5, EngineOffDesign6
from getatm import get_atmosphere_vec

#packages just needed for plotting since this is for sweeps
import matplotlib.pyplot as plt

"""
minimizes the aircraft total weight, must specify all weights except fuel weight, so in effect
we are minimizing the fuel weight
Rate of climb equation taken from John Anderson's Aircraft Performance and Design (eqn 5.85)
"""

#altitude precomputation
#select the cruise altitude
hcruise = 30000
#develop the list of altitudes
hvec = [3625, 7875, .25*(hcruise-10000)+10000, hcruise-.25*(hcruise-10000), hcruise, hcruise]
#convert from ft to m for atmosphere model
hvec = [x * 0.3048 for x in hvec]
#get the actual atmosphere values
atmdict = get_atmosphere_vec(hvec)
Tvec = atmdict['T']  * units('K')
rhovec = atmdict['rho'] * units('kg/m^3')
pvec = atmdict['p'] * units('kPa')

#TODO
#---------------------------------------------------------------------

class CommericalMissionConstraints(Model):
    """
    class that is general constraints that apply across the mission
    """
    def __init__(self, Nclimb1, Nclimb2, Ncruise2, signomial=0, **kwargs):
        #variable local to this model
        #10,000 foot altitude
        alt10k = Variable('alt10k', 10000, 'feet', '10,000 feet')

        #altitude variables
        dhTakeoff = Variable('dh_{takeoff}', 1500, 'feet', 'Total Altitude Change During Takeoff Profile')
        dhClimb1 = Variable('dh_{climb1}', 8500, 'feet', 'Total Altitude Change Required in Climb 1')
        dhClimb2 = Variable('dh_{climb2}', 'feet', 'Total Altitude Change Required in Climb 2')
        dhftClimb1 = VectorVariable(Nclimb1, 'dhftClimb1', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        dhftClimb2 = VectorVariable(Nclimb2, 'dhftClimb2', 'feet', 'Change in Altitude Per Climb Segment [feet]')
        hClimb1 = VectorVariable(Nclimb1, 'hClimb1', 'm', 'Altitude [meters]')
        hftClimb1 = VectorVariable(Nclimb1, 'hftClimb1', 'feet', 'Altitude [feet]')
        hClimb2 = VectorVariable(Nclimb2, 'hClimb2', 'm', 'Altitude [meters]')
        hftClimb2 = VectorVariable(Nclimb2, 'hftClimb2', 'feet', 'Altitude [feet]')
        hCruise2 = VectorVariable(Ncruise2, 'hCruise2', 'm', 'Altitude [meters]')
        hftCruise2 = VectorVariable(Ncruise2, 'hftCruise2', 'feet', 'Altitude [feet]')
        
        #alttidue at the top of climb
        htoc = Variable('h_{toc}', 'ft', 'Altitude at Top of Climb')

        #weight variables
        W_payload = Variable('W_{payload}', 'N', 'Aircraft Payload Weight')
        W_e = Variable('W_{e}', 'N', 'Empty Weight of Aircraft')
        W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
        W_ftotal = Variable('W_{f_{total}}', 'N', 'Total Fuel Weight')
        W_wing = Variable('W_{wing}', 'N', 'Wing Weight')
        W_total = Variable('W_{total}', 'N', 'Total Aircraft Weight')
        W_startClimb1 = VectorVariable(Nclimb1, 'W_{startClimb1}', 'N', 'Segment Start Weight')
        W_fuelClimb1 = VectorVariable(Nclimb1, 'W_{fuelClimb1}', 'N', 'Segment Fuel Weight')
        W_endClimb1 = VectorVariable(Nclimb1, 'W_{endClimb1}', 'N', 'Segment End Weight')
        W_avgClimb1 = VectorVariable(Nclimb1, 'W_{avgClimb1}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_startClimb2 = VectorVariable(Nclimb2, 'W_{startClimb2}', 'N', 'Segment Start Weight')
        W_fuelClimb2 = VectorVariable(Nclimb2, 'W_{fuelClimb2}', 'N', 'Segment Fuel Weight')
        W_endClimb2 = VectorVariable(Nclimb2, 'W_{endClimb2}', 'N', 'Segment End Weight')
        W_avgClimb2 = VectorVariable(Nclimb2, 'W_{avgClimb2}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_startCruise2 = VectorVariable(Ncruise2, 'W_{startCruise2}', 'N', 'Segment Start Weight')
        W_fuelCruise2 = VectorVariable(Ncruise2, 'W_{fuelCruise2}', 'N', 'Segment Fuel Weight')
        W_endCruise2 = VectorVariable(Ncruise2, 'W_{endCruise2}', 'N', 'Segment End Weight')
        W_avgCruise2 = VectorVariable(Ncruise2, 'W_{avgCruise2}', 'N', 'Geometric Average of Segment Start and End Weight')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #range variables
        ReqRngCruise = Variable('ReqRngCruise', 'miles', 'Required Cruise Range')
        ReqRng = Variable('ReqRng', 'miles', 'Required Mission Range')
        RngClimb1 = VectorVariable(Nclimb1, 'RngClimb1', 'miles', 'Segment Range During Climb1')
        RngClimb2 = VectorVariable(Nclimb2, 'RngClimb2', 'miles', 'Segment Range During Climb2')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        
        constraints = []
        constraints.extend([
            #convert m to ft
            hftClimb1  == hClimb1,
            hftClimb2  == hClimb2,
            hftCruise2  == hCruise2,
            
            
            #constraints on the various weights
            #with engine weight
            TCS([W_e + W_payload + W_ftotal + numeng * W_engine + W_wing <= W_total]),
 
            W_startClimb1[0]  == W_total,
   
            TCS([W_e + W_payload + numeng * W_engine + W_wing <= W_endCruise2[Ncruise2-1]]),
            TCS([W_ftotal >= sum(W_fuelClimb1) + sum(W_fuelClimb2) + sum(W_fuelCruise2)]),

            W_startClimb2[0] == W_endClimb1[Nclimb1-1],
            W_startCruise2[0] == W_endClimb2[Nclimb2-1],

            #wing weight constraint
            #based off of a raymer weight and 737 data from TASOPT output file
            (S/(124.58*units('m^2')))**.65 == W_wing/(105384.1524*units('N')),

            W_e == .75*W_payload,
            ])
        
        with gpkit.SignomialsEnabled():
            if signomial != 1:
                constraints.extend([
                    #range constraints
                    TCS([sum(RngClimb1) + sum(RngClimb2) + ReqRngCruise >= ReqRng]),
##                    dhClimb2==20000*units('ft'),
##                    TCS([dhClimb2 + alt10k >= htoc]),
                    ])
            if signomial ==1:
                 constraints.extend([
                     excessPtoc+Vtoc*Dtoc, numeng*Fd*Vtoc
##                    TCS([sum(RngClimb1) + sum(RngClimb2) >= ReqRng]),
##                    TCS([ReqRngCruise   >= sum(RngCruise)]),
##                    TCS([dhClimb2 + alt10k >= htoc])
                    ])
            for i in range(0, Nclimb1):
                constraints.extend([
                    SignomialEquality(hftClimb1[i], 1500*units('ft')+(i+1)*dhftClimb1[i]),
                    #constrain the geometric weight average
                    W_avgClimb1[i] == (W_startClimb1[i]*W_endClimb1[i])**.5,
                    TCS([W_startClimb1[i] >= W_endClimb1[i] + W_fuelClimb1[i]]),
                    ])

            for i in range(1, Nclimb1):
                constraints.extend([
                    TCS([W_startClimb1[i] == W_endClimb1[i-1]]),
                    ])

            for i in range(1, Nclimb2):
                constraints.extend([
                    TCS([W_startClimb2[i] == W_endClimb2[i-1]]),
                    ])

            for i in range(1, Ncruise2):
                constraints.extend([
                    TCS([W_startCruise2[i] == W_endCruise2[i-1]]),
                    ])

            for i in range(0, Nclimb2):
                constraints.extend([
                    SignomialEquality(hftClimb2[i], 10000*units('ft')+(i+1)*dhftClimb2[i]),
                    W_avgClimb2[i] == (W_startClimb2[i]*W_endClimb2[i])**.5,
                    TCS([W_startClimb2[i] >= W_endClimb2[i] + W_fuelClimb2[i]]),
                    ])

            for i in range(0, Ncruise2):
                constraints.extend([
                    W_avgCruise2[i] == (W_startCruise2[i]*W_endCruise2[i])**.5,
                    TCS([W_startCruise2[i] >= W_endCruise2[i] + W_fuelCruise2[i]]),
                    ])
        
        Model.__init__(self, W_ftotal, constraints, **kwargs)
        
#---------------------------------------
#takeoff

class Climb1(Model):
    """
    class to model the climb portion of a flight, applies to all climbs below
    10,000'
    """
    def __init__(self, Nclimb1, **kwargs):
        #set the speed limit under 10,000'
        speedlimit = Variable('speedlimit', 'kts', 'Speed Limit Under 10,000 ft')

        #atmosphere
        aClimb1 = VectorVariable(Nclimb1, 'aClimb1', 'm/s', 'Speed of Sound')
        rhoClimb1 = VectorVariable(Nclimb1, '\rhoClimb1', 'kg/m^3', 'Air Density')
        pClimb1 = VectorVariable(Nclimb1, 'pClimb1', 'kPa', 'Pressure')
        muClimb1 = VectorVariable(Nclimb1, '\muClimb1', 'kg/m/s', 'Air Kinematic Viscosity')
        TClimb1 = VectorVariable(Nclimb1, 'TClimb1', 'K', 'Air Temperature')

        #aero
        CLClimb1 = VectorVariable(Nclimb1, 'C_{L_{Climb1}}', '-', 'Lift Coefficient')
        WLoadClimb1 = VectorVariable(Nclimb1, 'W_{Load_{Climb1}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DClimb1 = VectorVariable(Nclimb1, 'DragClimb1', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cd0 = Variable('C_{d_0}', '-', 'Aircraft Cd0')
        K = Variable('K', '-', 'K for Parametric Drag Model')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #thrust
        thrustClimb1 = Variable('thrustClimb1', 'N', 'Engine Thrust')

        #excess power during climb
        excessPclimb1 = VectorVariable(Nclimb1, 'Excess Power Climb1', 'W', 'Excess Power During Climb')

        #climb rate and angle
        RCClimb1 = VectorVariable(Nclimb1, 'RCClimb1', 'feet/min', 'Rate of Climb/Decent')
        thetaClimb1 = VectorVariable(Nclimb1, '\\thetaClimb1', '-', 'Aircraft Climb Angle')

        #time
        tminClimb1 = VectorVariable(Nclimb1, 'tminClimb1', 'min', 'Flight Time in Minutes')
        thoursClimb1 = VectorVariable(Nclimb1, 'thrClimb1', 'hour', 'Flight Time in Hours')

        #range
        RngClimb1 = VectorVariable(Nclimb1, 'RngClimb1', 'miles', 'Segment Range During Climb1')

        #velocitites and mach numbers
        VClimb1 = VectorVariable(Nclimb1, 'VClimb1', 'knots', 'Aircraft Flight Speed')
        MClimb1 = VectorVariable(Nclimb1, 'MClimb1', '-', 'Aircraft Mach Number')

        #altitude
        dhftClimb1 = VectorVariable(Nclimb1, 'dhftClimb1', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        dhClimb1 = Variable('dh_{climb1}', 8500, 'feet', 'Total Altitude Change Required in Climb 1')

        W_avgClimb1 = VectorVariable(Nclimb1, 'W_{avgClimb1}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_fuelClimb1 = VectorVariable(Nclimb1, 'W_{fuelClimb1}', 'N', 'Segment Fuel Weight')

        
        #TEMPORARY
        mhold1 = Variable('mhold1', '-', 'segment 1 mach number')
        mhold2 = Variable('mhold2', '-', 'segment 2 mach number')
        Thold1 = Variable('Thold1', 'K', 'segment 1 T')
        Thold2 = Variable('Thold2', 'K', 'segment 2 T ')
        phold1 = Variable('phold1', 'kPa', 'segment 1 p')
        phold2 = Variable('phold2', 'kPa', 'segment 2 p')

##        TSFCc1 = VectorVariable(Nclimb1, 'TSFC_{c1}', '1/hr', 'Thrust Specific Fuel Consumption During Climb1')
        TSFCc11 = Variable('TSFC_{c11}', '1/hr', 'Thrust Specific Fuel Consumption During Climb1 part 1')
        TSFCc12 = Variable('TSFC_{c12}', '1/hr', 'Thrust Specific Fuel Consumption During Climb1 part 2')
##        thrustc1 = VectorVariable(Nclimb1, 'thrust_{c1}', 'N', 'Thrust During Climb Segment #1')
        thrustc11 = Variable('thrust_{c11}', 'N', 'Thrust During Climb Segment #1')
        thrustc12 = Variable('thrust_{c12}', 'N', 'Thrust During Climb Segment #1')

        #non-GPkit variables
        #cruise 2 lsit
        icl1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))
        
        constraints = []
            
        constraints.extend([            
            #set the velocity limits
            VClimb1[icl1] <= speedlimit,
            VClimb1[icl1] >= Vstall,

            MClimb1[0] == mhold1,
            MClimb1[1] == mhold2,

            MClimb1 * aClimb1 == VClimb1,

            #constraint on drag and thrust
##            thrustc11 >= D[iclimb1] + W_start[iclimb1]*theta[iclimb1],
            numeng*thrustc11 >= DClimb1[0] + W_avgClimb1[0]*thetaClimb1[0],
            numeng*thrustc12 >= DClimb1[1] + W_avgClimb1[1]*thetaClimb1[1],
            #climb rate constraints
##            TCS([excessP[iclimb1]+V[iclimb1]*D[iclimb1] <= V[iclimb1]*thrustc11]),
            TCS([excessPclimb1[0]+VClimb1[0]*DClimb1[0] <= VClimb1[0]*numeng*thrustc11]),
            TCS([excessPclimb1[1]+VClimb1[1]*DClimb1[1] <= VClimb1[1]*numeng*thrustc12]),
            
            TCS([DClimb1[icl1] >= (.5*S*rhoClimb1[icl1]*VClimb1[icl1]**2)*(Cd0 + K*CLClimb1[icl1]**2)]),
            RCClimb1[icl1] == excessPclimb1[icl1]/W_avgClimb1[icl1],
            RCClimb1[icl1] >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            thetaClimb1[icl1]*VClimb1[icl1]  == RCClimb1[icl1],
           
            dhftClimb1[icl1]  == tminClimb1[icl1] * RCClimb1[icl1],

            tminClimb1 == thoursClimb1,
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
##            TCS([RngClimb[iclimb1] + .5*thours[iclimb1]*V[iclimb1]*theta[iclimb1]**2 >= thours[iclimb1]*V[iclimb1]]),
            RngClimb1[icl1] == thoursClimb1[icl1]*VClimb1[icl1],

            W_avgClimb1[icl1] == .5*CLClimb1[icl1]*S*rhoClimb1[icl1]*VClimb1[icl1]**2,
            WLoadClimb1[icl1] == .5*CLClimb1[icl1]*S*rhoClimb1[icl1]*VClimb1[icl1]**2/S,
            
            #compute fuel burn from TSFC
##            W_fuel[iclimb1]  == TSFCc1[iclimb1] * thours[iclimb1] * thrust,
            W_fuelClimb1[0]  == numeng*TSFCc11 * thoursClimb1[0] * thrustc11,
            W_fuelClimb1[1]  == numeng*TSFCc12 * thoursClimb1[1] * thrustc12,
            #compute the dh required for each climb 1 segment
            dhftClimb1[icl1] == dhClimb1/Nclimb1,

            #constrain the max wing loading
            WLoadClimb1 <= WLoadmax,
            ])

        for i in range(0, Nclimb1):
            constraints.extend([
                rhoClimb1[i] == rhovec[i],
                TClimb1[i] == Tvec[i],
                pClimb1[i] == pvec[i],
                #speed of sound
                aClimb1[i]  == (gamma * R * TClimb1[i])**.5,

                phold1 == pClimb1[0],
                phold2 == pClimb1[1],
                Thold1 == TClimb1[0],
                Thold2 == TClimb1[1],
                ])
        
        Model.__init__(self, None, constraints, **kwargs)
        
#--------------------------------------
#transition between climb 1 and 2 if desired

            
class Climb2(Model):
    """
    class to model the climb portion above 10,000'
    """
    def __init__(self, Nclimb2, **kwargs):
        #aero
        CLClimb2 = VectorVariable(Nclimb2, 'C_{L_{Climb2}}', '-', 'Lift Coefficient')
        WLoadClimb2 = VectorVariable(Nclimb2, 'W_{Load_{Climb2}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DClimb2 = VectorVariable(Nclimb2, 'DragClimb2', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cd0 = Variable('C_{d_0}', '-', 'Aircraft Cd0')
        K = Variable('K', '-', 'K for Parametric Drag Model')

        #atmosphere
        aClimb2 = VectorVariable(Nclimb2, 'aClimb2', 'm/s', 'Speed of Sound')
        rhoClimb2 = VectorVariable(Nclimb2, '\rhoClimb2', 'kg/m^3', 'Air Density')
        pClimb2 = VectorVariable(Nclimb2, 'pClimb2', 'kPa', 'Pressure')
        muClimb2 = VectorVariable(Nclimb2, '\muClimb2', 'kg/m/s', 'Air Kinematic Viscosity')
        TClimb2 = VectorVariable(Nclimb2, 'TClimb2', 'K', 'Air Temperature')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #thrust
        thrustClimb2 = Variable('thrustClimb2', 'N', 'Engine Thrust')

        #excess power during climb
        excessPclimb2 = VectorVariable(Nclimb2, 'Excess Power Climb2', 'W', 'Excess Power During Climb')

        #climb rate and angle
        RCClimb2 = VectorVariable(Nclimb2, 'RCClimb2', 'feet/min', 'Rate of Climb/Decent')
        thetaClimb2 = VectorVariable(Nclimb2, '\\thetaClimb2', '-', 'Aircraft Climb Angle')

        #time
        tminClimb2 = VectorVariable(Nclimb2, 'tminClimb2', 'min', 'Flight Time in Minutes')
        thoursClimb2 = VectorVariable(Nclimb2, 'thrClimb2', 'hour', 'Flight Time in Hours')

        #range
        RngClimb2 = VectorVariable(Nclimb2, 'RngClimb2', 'miles', 'Segment Range During Climb2')

        #velocitites and mach numbers
        VClimb2 = VectorVariable(Nclimb2, 'VClimb2', 'knots', 'Aircraft Flight Speed')
        MClimb2 = VectorVariable(Nclimb2, 'MClimb2', '-', 'Aircraft Mach Number')

        #altitude
        dhftClimb2 = VectorVariable(Nclimb2, 'dhftClimb2', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
        dhClimb2 = Variable('dh_{climb2}', 'feet', 'Total Altitude Change Required in Climb 2')

        W_fuelClimb2 = VectorVariable(Nclimb2, 'W_{fuelClimb2}', 'N', 'Segment Fuel Weight')
        W_avgClimb2 = VectorVariable(Nclimb2, 'W_{avgClimb2}', 'N', 'Geometric Average of Segment Start and End Weight')

        #TEMPORARY
        mhold3 = Variable('mhold3', '-', 'segment 3 mach number')
        mhold4 = Variable('mhold4', '-', 'segment 4 mach number')
        Thold3 = Variable('Thold3', 'K', 'segment 3 T')
        Thold4 = Variable('Thold4', 'K', 'segment 4 T')
        phold3 = Variable('phold3', 'kPa', 'segment 3 p')
        phold4 = Variable('phold4', 'kPa', 'segment 4 p')
        
##        TSFCc2 = VectorVariable(Nclimb2, 'TSFC_{c2}', '1/hr', 'Thrust Specific Fuel Consumption During Climb2')
        TSFCc21 = Variable('TSFC_{c21}', '1/hr', 'Thrust Specific Fuel Consumption During Climb2')
        TSFCc22 = Variable('TSFC_{c22}', '1/hr', 'Thrust Specific Fuel Consumption During Climb2')
##        thrustc2 = VectorVariable(Nclimb1, 'thrust_{c2}', 'N', 'Thrust During Climb Segment #2')
        thrustc21 = Variable('thrust_{c21}', 'N', 'Thrust During Climb Segment #2')
        thrustc22 = Variable('thrust_{c22}', 'N', 'Thrust During Climb Segment #2')

        #non-GPkit variables
        #climb 2 lsit
        icl2 = map(int, np.linspace(0, Nclimb2 - 1, Nclimb2))

        constraints = []
        
        constraints.extend([            
            #set the velocity limits
            #needs to be replaced by an actual Vne and a mach number
            MClimb2[icl2] <= .75,
            VClimb2[icl2] >= Vstall,
            
            MClimb2[0] == mhold3,
            MClimb2[1] == mhold4,

            VClimb2 == MClimb2 * aClimb2,

            #constraint on drag and thrust
##            thrustc21 >= D[iclimb2] + W_start[iclimb2]*theta[iclimb2],
            numeng*thrustc21 >= DClimb2[0] + W_avgClimb2[0]*thetaClimb2[0],
            numeng*thrustc22 >= DClimb2[1] + W_avgClimb2[1]*thetaClimb2[1],
            
            #climb rate constraints
##            TCS([excessP[iclimb2]+V[iclimb2]*D[iclimb2] <= V[iclimb2]*thrustc21]),
            TCS([excessPclimb2[0]+VClimb2[0]*DClimb2[0] <= VClimb2[0]*numeng*thrustc21]),
            TCS([excessPclimb2[1]+VClimb2[1]*DClimb2[1] <= VClimb2[1]*numeng*thrustc22]),
            TCS([DClimb2[icl2] >= (.5*S*rhoClimb2[icl2]*VClimb2[icl2]**2)*(Cd0 + K*CLClimb2[icl2]**2)]),
            RCClimb2[icl2] == excessPclimb2[icl2]/W_avgClimb2[icl2],
            RCClimb2[0] >= 500*units('ft/min'),
            RCClimb2[1] >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            thetaClimb2[icl2]*VClimb2[icl2]  == RCClimb2[icl2],
           
            dhftClimb2[icl2]  == tminClimb2[icl2] * RCClimb2[icl2],

            tminClimb2 == thoursClimb2,
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
##            TCS([RngClimb[iclimb2] + .5*thours[iclimb2]*V[iclimb2]*theta[iclimb2]**2 <= thours[iclimb2]*V[iclimb2]]),
            RngClimb2[icl2] <= thoursClimb2[icl2]*VClimb2[icl2],

            W_avgClimb2[icl2] == .5*CLClimb2[icl2]*S*rhoClimb2[icl2]*VClimb2[icl2]**2,      
            WLoadClimb2[icl2] == .5*CLClimb2[icl2]*S*rhoClimb2[icl2]*VClimb2[icl2]**2/S,
            
            #compute fuel burn from TSFC
            W_fuelClimb2[0]  == numeng*TSFCc21 * thoursClimb2[0] * thrustc21,
            W_fuelClimb2[1]  == numeng*TSFCc22 * thoursClimb2[1] * thrustc22,

            #compute the dh required for each climb 1 segment
            dhftClimb2[icl2] == dhClimb2/Nclimb2,

            #constrain the max wing loading
            WLoadClimb2 <= WLoadmax,
            ])

        for i in range(0, Nclimb2):
            constraints.extend([
                rhoClimb2[i] == rhovec[i],
                TClimb2[i] == Tvec[i],
                pClimb2[i] == pvec[i],
                #speed of sound
                aClimb2[i]  == (gamma * R * TClimb2[i])**.5,

                phold3 == pClimb2[0],
                phold4 == pClimb2[1],
                Thold3 == TClimb2[0],
                Thold4 == TClimb2[1],
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
    def __init__(self, Nclimb2, Ncruise2, **kwargs):
        """
    class to model the second cruise portion of a flight (if the flight is
    long enough to mandate two cruise portions)
    Model is based off of a discretized Breguet Range equation
    """
        #aero
        CLCruise2 = VectorVariable(Ncruise2, 'C_{L_{Cruise2}}', '-', 'Lift Coefficient')
        WLoadCruise2 = VectorVariable(Ncruise2, 'W_{Load_{Cruise2}}', 'N/m^2', 'Wing Loading')
        WLoadmax = Variable('W_{Load_max}', 'N/m^2', 'Max Wing Loading')
        DCruise2 = VectorVariable(Ncruise2, 'DragCruise2', 'N', 'Drag')
        Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
        Cd0 = Variable('C_{d_0}', '-', 'Aircraft Cd0')
        K = Variable('K', '-', 'K for Parametric Drag Model')

        #atmosphere
        aCruise2 = VectorVariable(Ncruise2, 'aCruise2', 'm/s', 'Speed of Sound')
        rhoCruise2 = VectorVariable(Ncruise2, '\rhoCruise2', 'kg/m^3', 'Air Density')
        pCruise2 = VectorVariable(Ncruise2, 'pCruise2', 'kPa', 'Pressure')
        muCruise2 = VectorVariable(Ncruise2, '\muCruise2', 'kg/m/s', 'Air Kinematic Viscosity')
        TCruise2 = VectorVariable(Ncruise2, 'TCruise2', 'K', 'Air Temperature')

        #aircraft geometry
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #number of engines
        numeng = Variable('numeng', '-', 'Number of Engines')

        #time
        tminCruise2 = VectorVariable(Ncruise2, 'tminCruise2', 'min', 'Flight Time in Minutes')
        thoursCruise2 = VectorVariable(Ncruise2, 'thrCruise2', 'hour', 'Flight Time in Hours')

        #thrust
        thrustCruise2 = Variable('thrustCruise2', 'N', 'Engine Thrust')

        #velocitites and mach numbers
        VCruise2 = VectorVariable(Ncruise2, 'VCruise2', 'knots', 'Aircraft Flight Speed')
        MCruise2 = VectorVariable(Ncruise2, 'MCruise2', '-', 'Aircraft Mach Number')

        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')
        
##        TSFCcr2 = VectorVariable(Ncruise2, 'TSFC_{cr2}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        TSFCcr21 = Variable('TSFC_{cr21}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        TSFCcr22 = Variable('TSFC_{cr22}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        D1 = Variable('D1', 'N', 'Drag for cruise part 1')
        D2 = Variable('D2', 'N', 'Drag for cruise part 2')

        thrustcr21 = Variable('thrust_{cr21}', 'N', 'Thrust During Cruise Segment #2')
        thrustcr22 = Variable('thrust_{cr22}', 'N', 'Thrust During Cruise Segment #2')

        constraints = []
        
        #defined here for linking purposes
        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')
        Fd = Variable('F_D', 'N', 'Design Thrust')

        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        #parameter to make breguent range eqn gp compatible
        z_brec2 = VectorVariable(Ncruise2, 'z_{brec2}', '-', 'Breguet Parameter')

        #range variables
        ReqRngCruise = Variable('ReqRngCruise', 'miles', 'Required Cruise Range')
        RngCruise2 = VectorVariable(Ncruise2, 'RngCruise2', 'miles', 'Segment Range During Cruise2')

        #top of climb parameters for engine sizing
        RCtoc = Variable('RC @ TOC', 'feet/min', 'Rate of Climb at Top of Climb')
        Dtoc = Variable('Drag @ TOC', 'N', 'Drag at Top of Climb')
        excessPtoc = Variable('Excess Power @ TOC', 'W', 'Excess Power at Top Of Climb')
        Vtoc = Variable('V @ TOC', 'knots', 'Aircraft Flight Speed at TOC')
        htoc = Variable('h_{toc}', 'ft', 'Altitude at Top of Climb')
        CLtoc = Variable('C_{Ltoc}', '-', 'Lift Coefficient @ TOC')

        #altitude
        hCruise2 = VectorVariable(Ncruise2, 'hCruise2', 'm', 'Altitude [meters]')
        hftCruise2 = VectorVariable(Ncruise2, 'hftCruise2', 'feet', 'Altitude [feet]')

        W_fuelCruise2 = VectorVariable(Ncruise2, 'W_{fuelCruise2}', 'N', 'Segment Fuel Weight')
        W_avgCruise2 = VectorVariable(Ncruise2, 'W_{avgCruise2}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_avgClimb2 = VectorVariable(Nclimb2, 'W_{avgClimb2}', 'N', 'Geometric Average of Segment Start and End Weight')
        W_endCruise2 = VectorVariable(Ncruise2, 'W_{endCruise2}', 'N', 'Segment End Weight')

        #TEMPORARY
        mhold5 = Variable('mhold5', '-', 'segment 5 mach number')
        mhold6 = Variable('mhold6', '-', 'segment 6 mach number')
        Thold5 = Variable('Thold5', 'K', 'segment 5 T')
        Thold6 = Variable('Thold6', 'K', 'segment 6 T')
        phold5 = Variable('phold5', 'kPa', 'segment 5 p')
        phold6 = Variable('phold6', 'kPa', 'segment 6 p')

        #non-GPkit variables
        #cruise 2 lsit
        izbre = map(int, np.linspace(0, Ncruise2 - 1, Ncruise2))
            
        constraints.extend([
            MCruise2[izbre] == 0.8,

            MCruise2 * aCruise2 == VCruise2,
            
##                P0 == p[Nclimb],
            T0 == 280*units('K'),
            
             #climb rate constraints for engine sizing at TOC
            excessPtoc+Vtoc*Dtoc  <= numeng*Fd*Vtoc,
            RCtoc == excessPtoc/W_avgClimb2[Nclimb2-1],
            RCtoc == 500*units('ft/min'),
            Vtoc == VCruise2[0],

            #compute the drag
            TCS([Dtoc >= (.5*S*rhoCruise2[0]*Vtoc**2)*(Cd0 + K*(W_avgClimb2[Nclimb2-1]/(.5*S*rhoCruise2[0]*Vtoc**2))**2)]),
##                TCS([D[icruise2] >= (.5*S*rho[icruise2]*V[icruise2]**2)*(Cd0 + K*(W_start[icruise2]/(.5*S*rho[icruise2]*V[icruise2]**2))**2)]),
            TCS([DCruise2[0] >= (.5*S*rhoCruise2[0]*VCruise2[0]**2)*(Cd0 + K*(W_avgCruise2[0]/(.5*S*rhoCruise2[0]*VCruise2[0]**2))**2)]),
            TCS([DCruise2[1] >= (.5*S*rhoCruise2[1]*VCruise2[1]**2)*(Cd0 + K*(W_avgCruise2[1]/(.5*S*rhoCruise2[1]*VCruise2[1]**2))**2)]),
            DCruise2[0] == numeng * thrustcr21,
            DCruise2[1] == numeng * thrustcr22,

            W_avgCruise2[izbre] == .5*CLCruise2[izbre]*S*rhoCruise2[izbre]*VCruise2[izbre]**2,
            WLoadCruise2[izbre] == .5*CLCruise2[izbre]*S*rhoCruise2[izbre]*VCruise2[izbre]**2/S,
            
            #constrain the climb rate by holding altitude constant
            hftCruise2[izbre]  == htoc,
            
            #taylor series expansion to get the weight term
            TCS([W_fuelCruise2[izbre]/W_endCruise2[izbre] >= te_exp_minus1(z_brec2[izbre], nterm=3)]),

            #breguet range eqn
            TCS([z_brec2[0] >= (numeng*TSFCcr21*thoursCruise2[0]*DCruise2[0])/W_avgCruise2[0]]),
            TCS([z_brec2[1] >= (numeng*TSFCcr22*thoursCruise2[1]*DCruise2[1])/W_avgCruise2[1]]),

            #time
            thoursCruise2[izbre]*VCruise2[izbre]  == RngCruise2[izbre],
            tminCruise2 == thoursCruise2,

            #constrain the max wing loading
            WLoadCruise2 <= WLoadmax,
            ])
        
        #constraint on the aircraft meeting the required range
        for i in range(min(izbre), max(izbre)+1):
            constraints.extend([
                TCS([RngCruise2[i] == ReqRngCruise/(Ncruise2)])
                ])

        for i in range(0, Ncruise2):
            constraints.extend([
                rhoCruise2[i] == rhovec[i+4],
                TCruise2[i] == Tvec[i+4],
                pCruise2[i] == pvec[i+4],
                #speed of sound
                aCruise2[i]  == (gamma * R * TCruise2[i])**.5,

                phold5 == pCruise2[0],
                phold6 == pCruise2[1],
                Thold5 == TCruise2[0],
                Thold6 == TCruise2[1],
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
        #defining the number of segments
        Ntakeoff = 0
        Nclimb1 = 2
        Nclimb2 = 2
        Ncruise1 = 0
        Ncruiseclimb = 0
        Ncruise2 = 2
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

        Nseg = Nclimb + Ncruise + Nres + Ntakeoff + Ndecent + Nlanding

        #define all the submodels
        cmc = CommericalMissionConstraints(Nclimb1, Nclimb2, Ncruise2, 0)
        climb1 = Climb1(Nclimb1)
        climb2 = Climb2(Nclimb2)
        cruise2 = Cruise2(Nclimb2, Ncruise2)
##        atm = Atmosphere(Nseg)
        eonD = EngineOnDesign()
        eoffD = EngineOffDesign()
        eoffD2 = EngineOffDesign2()
        eoffD3 = EngineOffDesign3()
        eoffD4 = EngineOffDesign4()
        eoffD5 = EngineOffDesign5()
        eoffD6 = EngineOffDesign6()
        for i in range(Nseg):
            None
 
        substitutions = {      
            'W_{payload}': .6*44000*9.8*units('N'),
            'V_{stall}': 120,
            'ReqRng': 1000,
            'C_{d_0}': .02,
            'K': 0.05,
            'h_{toc}': hcruise,
            'speedlimit': 250,
            'numeng': 2,
            'dh_{climb2}': hcruise-10000,
            'W_{Load_max}': 6664,

            #substitutions for global engine variables
            'G_f': 1,
            'N_{{bar}_Df}': 1,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\pi_{d}': .99,
            '\pi_{fn}': .98,
            '\pi_{tn}': .99,
            '\pi_{b}': .94,
            '\pi_{f_D}': 1.5,
            '\pi_{lc_D}': 3,
            '\pi_{hc_D}': 10
            }
        #for engine on design must link T0, P0, F_D,TSFC w/TSFC from icruise 2
        
        self.submodels = [cmc, climb1, climb2, cruise2, eonD, eoffD, eoffD2, eoffD3, eoffD4, eoffD5, eoffD6]

        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_{c11}_Climb1': 'TSFC_E_EngineOffDesign', 'thrust_{c11}_Climb1': 'F_EngineOffDesign',
                                'TSFC_{c12}_Climb1': 'TSFC_E2_EngineOffDesign2', 'thrust_{c12}_Climb1': 'F_2_EngineOffDesign2',
                                'thrust_{c21}_Climb2': 'F_3_EngineOffDesign3','TSFC_{c21}_Climb2': 'TSFC_E3_EngineOffDesign3',
                                'thrust_{c22}_Climb2': 'F_4_EngineOffDesign4','TSFC_{c22}_Climb2': 'TSFC_E4_EngineOffDesign4',
                                'TSFC_{cr21}_Cruise2': 'TSFC_E5_EngineOffDesign5', 'thrust_{cr21}': 'F_{spec5}_EngineOffDesign5',
                                'TSFC_{cr22}_Cruise2': 'TSFC_E6_EngineOffDesign6', 'thrust_{cr22}': 'F_{spec6}_EngineOffDesign6',
                                'mhold1': 'M_0_1', 'mhold2': 'M_0_2', 'mhold3': 'M_0_3', 'mhold4': 'M_0_4', 'phold1': 'P_0_1',
                                'phold2': 'P_0_2','phold3': 'P_0_3','phold4': 'P_0_4','phold5': 'P_0_5','phold6': 'P_0_6',
                                'P_0': 'P_0_6'})#, 'Thold1': 'T_0_1'})#,
##                                'Thold2': 'T_0_2','Thold3': 'T_0_3','Thold4': 'T_0_4','Thold5': 'T_0_5','Thold6': 'T_0_6',})

        lc = LinkedConstraintSet(constraints, exclude={'T_0', 'P_0', 'M_0', 'a_0', 'u_0', 'P_{t_0}', 'T_{t_0}', 'h_{t_0}', 'P_{t_1.8}',
                                                       'T_{t_1.8}', 'h_{t_1.8}', 'P_{t_2}', 'T_{t_2}', 'h_{t_2}', 'P_{t_2.1}','T_{t_2.1}',
                                                       'h_{t_2.1}', 'P_{t_2.5}', 'T_{t_2.5}', 'h_{t_2.5}', 'P_{t_3}', 'T_{t_3}',
                                                       'h_{t_3}','P_{t_7}', 'T_{t_7}', 'h_{t_7}', '\pi_f','\pi_{lc}','\pi_{hc}', 'P_{t_4}'
                                                       ,'h_{t_4}','T_{t_4}','P_{t_4.1}','T_{t_4.1}','h_{t_4.1}','f','h_{t_4.5}',
                                                       'P_{t_4.5}','T_{t_4.5}','P_{t_4.9}','T_{t_4.9}','h_{t_4.9}','P_{t_5}','T_{t_5}',
                                                       'h_{t_5}','\pi_{HPT}','\pi_{LPT}','P_8','P_{t_8}','h_{t_8}','h_8','T_{t_8}',
                                                       'T_{8}','P_6','P_{t_6}','T_{t_6}','T_{6}','h_{t_6}','h_6','F_8','F_6','F','F_{sp}'
                                                       ,'I_{sp}','TSFC_E','u_6','u_8','m_{core}','T_2','P_2','u_2','\rho_2.5','T_{2.5}',
                                                       'P_{2.5}','P_{t_2.5}','u_{2.5}','T_{t_2.5}','h_{t_2.5}','h_{2.5}','M_8','M_6',
                                                       'M_5','M_7','M_2','M_{2.5}','F_{sp}','T_{2}','h_{2}','T_{6}','T_{8}','T_{5}',
                                                       'T_{7}','P_{5}','P_0','fp1','u_7','\rho_7','u_5','\rho_5','m_{f}','m_{fan}',
                                                       'm_{tild_f}','p_{tildf}','N_f','m_{{tild}_sf}','p_{{tild}_sf}','N_1','m_{lc}',
                                                        'm_{hc}','u_7','M_7','\rho_7','P_{5}','T_{5}','M_5','u_5','\rho_5','T_{t_{4spec}}'
                                                        ,'m_{tild_hc}','p_{tild_lc}','N_2','m_{{tild}_shc}','p_{{tild}_shc}','m_{tild_lc}'
                                                        ,'p_{tild_lc}','m_{{tild}_slc}','p_{{tild}_slc}','P_{7}', 'F_{spec}'})

        Model.__init__(self, cmc.cost, lc, substitutions, **kwargs)

    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            constraints.append([ub*units >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= lb*units])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m


    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver=None, verbosity=0,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver, verbosity, **kwargs)
        solhold = sol
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out, solhold

    
if __name__ == '__main__':
    m = CommercialAircraft()
    sol = m.localsolve(solver="mosek", verbosity = 4, iteration_limit=100, skipsweepfailures=True)
    
##    sol, solhold = m.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=100)
    
