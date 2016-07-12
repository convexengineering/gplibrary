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
from engine_atm import EngineOnDesign, EngineOffDesign

"""
minimizes the aircraft total weight, must specify all weights except fuel weight, so in effect
we are minimizing the fuel weight
Rate of climb equation taken from John Anderson's Aircraft Performance and Design (eqn 5.85)
"""

#TODO
#link in with engine
#fix indicies on the TSFC vector variables

#-----------------------------------------------------------
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

#create index lists for the vector variables

#partial flight profile for use during development
if Nclimb1 != 0:
    iclimb1 = map(int, np.linspace(0, Nclimb1 - 1, Nclimb1))
else:
    iclimb1 = 0
if Nclimb2 != 0:
    iclimb2 = map(int, np.linspace(iclimb1[len(iclimb1)-1] + 1, Nclimb2 + iclimb1[len(iclimb1)-1], Nclimb2))
else:
    iclimb2 =0
    
if Ncruise2 != 0:
    icruise2 = map(int, np.linspace(iclimb2[len(iclimb2)-1] + 1, Ncruise2 + iclimb2[len(iclimb2)-1], Ncruise2))
##    icruise2 = map(int, np.linspace(0, Ncruise2-1, Ncruise2))
else:
    icruise2 = 0

izbre = map(int, np.linspace(0, Ncruise - 1, Ncruise))

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
if Nclimb != 0:
    dhft = VectorVariable(Nclimb, 'dhft', 'feet', 'Change in Altitude Per Climb Segment [feet]') 
h = VectorVariable(Nseg, 'h', 'm', 'Altitude [meters]')
hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude [feet]')
htoc = Variable('h_{toc}', 'ft', 'Altitude at Top of Climb')
dhClimb1 = Variable('dh_{climb1}', 8500, 'feet', 'Total Altitude Change Required in Climb 1')
dhClimb2 = Variable('dh_{climb2}', 'feet', 'Total Altitude Change Required in Climb 2')
dhTakeoff = Variable('dh_{takeoff}', 1500, 'feet', 'Total Altitude Change During Takeoff Profile')

#time
tmin = VectorVariable(Nseg, 'tmin', 'min', 'Flight Time in Minutes')
thours = VectorVariable(Nseg, 'thr', 'hour', 'Flight Time in Hours')

#Range
if Nclimb != 0:
    RngClimb = VectorVariable(Nclimb, 'RngClimb', 'miles', 'Segment Range During Climb')
if Ncruise !=0:
    RngCruise = VectorVariable(Ncruise2 + Ncruise1, 'RngCruise', 'miles', 'Segment Range During Cruise')
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

if Nclimb != 0:
    excessP = VectorVariable(Nclimb, 'Excess Power', 'W', 'Excess Power During Climb')
D = VectorVariable(Nseg, 'Drag', 'N', 'Drag')
Dtoc = Variable('Drag @ TOC', 'N', 'Drag at Top of Climb')
excessPtoc = Variable('Excess Power @ TOC', 'W', 'Excess Power at Top Of Climb')

#aircraft geometry
S = Variable('S', 'm^2', 'Wing Planform Area')

#velocitites and mach numbers
V = VectorVariable(Nseg, 'V', 'knots', 'Aircraft Flight Speed')
M = VectorVariable(Nseg, 'M', '-', 'Aircraft Mach Number')
Vstall = Variable('V_{stall}', 'knots', 'Aircraft Stall Speed')
Vtoc = Variable('V @ TOC', 'knots', 'Aircraft Flight Speed at TOC')

#Climb/Decent Rate
if Nclimb != 0:
    RC = VectorVariable(Nclimb, 'RC', 'feet/min', 'Rate of Climb/Decent')
    theta = VectorVariable(Nclimb, '\\theta', '-', 'Aircraft Climb Angle')
RCtoc = Variable('RC @ TOC', 'feet/min', 'Rate of Climb at Top of Climb')

#breguet parameter
if Ncruise != 0:
    z_bre = VectorVariable(Ncruise, 'z_{bre}', '-', 'Breguet Parameter')

#engine
#TSFC = VectorVariable(Nseg, 'TSFC', '1/hr', 'Thrust Specific Fuel Consumption')

#currently sets the value of TSFC, just a place holder
c1 = Variable('c1', '-', 'Constant')

#defien thrust variable
thrust = Variable('thrust', 'N', 'Engine Thrust')

#define the amount of thrust required at the sizing point
thrustsizing = Variable('thrustsizing', 'N', 'Engine Thrust for Sizing')

#set the speed limit under 10,000'
speedlimit = Variable('speedlimit', 250, 'kts', 'Speed Limit Under 10,000 ft')

#temporary args
#pass in a 1 for testing climb segment 1

class CommericalMissionConstraints(Model):
    """
    class that is general constraints that apply across the mission
    """
    def __init__(self, test=0, **kwargs):
        #variable local to this model
        alt10k = Variable('alt10k', 10000, 'feet', '10,000 feet')
        
        constraints = []
        constraints.extend([
            #speed of sound
            a  == (gamma * R * T)**.5,

            #comptue the mach number
            M == V/a,
            
            #convert m to ft
            hft  == h,
            
            #convert min to hours
            tmin  == thours ,
            
            #constraints on the various weights
            TCS([W_e + W_payload + W_ftotal <= W_total]),
            W_start[0]  == W_total,
            TCS([W_e + W_payload <= W_end[Nseg-1]]),
            TCS([W_ftotal >= sum(W_fuel)]),

            rho[iclimb1] == 1.225*units('kg/m^3'),
            T[iclimb1] == 273*units('K'),
            rho[iclimb2] == .7*units('kg/m^3'),
            T[iclimb2] == 250*units('K'),
            rho[icruise2] == .3*units('kg/m^3'),
            T[icruise2] == 230*units('K'),
            ])
        
        with gpkit.SignomialsEnabled():
            if test != 1:
                constraints.extend([
                    #range constraints
                    TCS([sum(RngClimb) + ReqRngCruise >= ReqRng]),
                    TCS([dhClimb2 + alt10k >= htoc])
                    ])
            if test ==1:
                 constraints.extend([
                    TCS([sum(RngCruise) >= ReqRng]),
##                    TCS([ReqRngCruise   >= sum(RngCruise)]),
##                    SignomialEquality(dhClimb2, htoc - 10000*units('ft'))
                    ])
            for i in range(0, Nclimb1):
                constraints.extend([
                    SignomialEquality(hft[i], 1500*units('ft')+(i+1)*dhft[i])
                    ])

            for i in range(0, Nclimb2):
                constraints.extend([
                    SignomialEquality(hft[i+Nclimb1], 10000*units('ft')+(i+1)*dhft[i+Nclimb1])
                    ])

        #constrain the segment weights in a loop
        for i in range(1, Nseg):
            constraints.extend([
                W_start[i] == W_end[i-1] 
                ])
        for i in range(0, Nseg):
            constraints.extend([
                TCS([W_start[i] >= W_end[i] + W_fuel[i]])
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
        TSFCc1 = VectorVariable(Nclimb1, 'TSFC_{c1}', '1/hr', 'Thrust Specific Fuel Consumption During Climb1')
        #Climb #1 (sub 10K subject to 250KTS speed limit)
        constraints = []

        constraints.extend([            
            #set the velocity limits
            V[iclimb1] <= speedlimit,
            V[iclimb1] >= Vstall,

            #constraint on drag and thrust
            thrust >= D[iclimb1] + W_start[iclimb1]*theta[iclimb1],
            
            #climb rate constraints
            TCS([excessP[iclimb1]+V[iclimb1]*D[iclimb1] <= V[iclimb1]*thrust]),
            TCS([D[iclimb1] >= (.5*S*rho[iclimb1]*V[iclimb1]**2)*(Cd0 + K*(W_start[iclimb1]/(.5*S*rho[iclimb1]*V[iclimb1]**2))**2)]),
            RC[iclimb1] == excessP[iclimb1]/W_start[iclimb1],
            
            #make the small angle approximation and compute theta
            theta[iclimb1]*V[iclimb1]  == RC[iclimb1],
           
            dhft[iclimb1]  == tmin[iclimb1] * RC[iclimb1],
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
            TCS([RngClimb[iclimb1] + .5*thours[iclimb1]*V[iclimb1]*theta[iclimb1]**2 <= thours[iclimb1]*V[iclimb1]]),
            
            #compute fuel burn from TSFC
            W_fuel[iclimb1]  == TSFCc1[iclimb1] * thours[iclimb1] * thrust,

            #compute the dh required for each climb 1 segment
            dhft[iclimb1] == dhClimb1/Nclimb1,

            TSFCc1[iclimb1] == c1*units('1/hr')
            ])
            
        Model.__init__(self, None, constraints, **kwargs)
        
#--------------------------------------
#transition between climb 1 and 2 if desired

            
class Climb2(Model):
    """
    class to model the climb portion above 10,000'
    """
    def __init__(self, **kwargs):
        TSFCc2 = VectorVariable(Nclimb2, 'TSFC_{c2}', '1/hr', 'Thrust Specific Fuel Consumption During Climb2')
        
        constraints = []
        
        constraints.extend([            
            #set the velocity limits
            #needs to be replaced by an actual Vne and a mach number
            M[iclimb2]<= 1,
            V[iclimb2] >= Vstall,

            #constraint on drag and thrust
            thrust >= D[iclimb2] + W_start[iclimb2]*theta[iclimb2],
            
            #climb rate constraints
            TCS([excessP[iclimb2]+V[iclimb2]*D[iclimb2] <= V[iclimb2]*thrust]),
            TCS([D[iclimb2] >= (.5*S*rho[iclimb2]*V[iclimb2]**2)*(Cd0 + K*(W_start[iclimb2]/(.5*S*rho[iclimb2]*V[iclimb2]**2))**2)]),
            RC[iclimb2] == excessP[iclimb2]/W_start[iclimb2],
            RC[iclimb2] >= 500*units('ft/min'),
            
            #make the small angle approximation and compute theta
            theta[iclimb2]*V[iclimb2]  == RC[iclimb2],
           
            dhft[iclimb2]  == tmin[iclimb2] * RC[iclimb2],
            #compute the distance traveled for each segment
            #takes into account two terms of a cosine expansion
            TCS([RngClimb[iclimb2] + .5*thours[iclimb2]*V[iclimb2]*theta[iclimb2]**2 <= thours[iclimb2]*V[iclimb2]]),
            
            #compute fuel burn from TSFC
            W_fuel[iclimb2]  == TSFCc2[iclimb1] * thours[iclimb2] * thrust,

            #compute the dh required for each climb 1 segment
            dhft[iclimb2] == dhClimb2/Nclimb2,

            TSFCc2[iclimb1] == c1*units('1/hr')
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
        """
    class to model the second cruise portion of a flight (if the flight is
    long enough to mandate two cruise portions)
    Model is based off of a discretized Breguet Range equation
    """
    def __init__(self, **kwargs):
##        TSFCcr2 = VectorVariable(Ncruise2, 'TSFC_{cr2}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        TSFCcr21 = Variable('TSFC_{cr21}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        TSFCcr22 = Variable('TSFC_{cr22}', '1/hr', 'Thrust Specific Fuel Consumption During Cruise2')
        
        constraints = []
        
        #defined here for linking purposes
        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')
        Fd = Variable('F_D', 'N', 'Design Thrust')

        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        
        with gpkit.SignomialsEnabled():
            
            constraints.extend([
                #constrain flight speeds with drag
                TCS([thrust >= D[icruise2]]),

##                P0 == p[Nclimb],
                T0 == T[Nclimb],
                
                #climb rate constraints for engine sizing at TOC
                SignomialEquality(excessPtoc+Vtoc*Dtoc, Fd*Vtoc),
                RCtoc == excessPtoc/W_start[Nclimb],
                RCtoc == 1000*units('ft/min'),
                Vtoc == V[icruise2],

                #compute the drag
                SignomialEquality(Dtoc, (.5*S*rho[Ncruise2]*Vtoc**2)*(Cd0 + K*(W_start[Nclimb]/(.5*S*rho[Ncruise2]*Vtoc**2))**2)),
                TCS([D[icruise2] >= (.5*S*rho[icruise2]*V[icruise2]**2)*(Cd0 + K*(W_start[icruise2]/(.5*S*rho[icruise2]*V[icruise2]**2))**2)]),
                
                #constrain the climb rate by holding altitude constant
                hft[icruise2]  == htoc,
                
                #taylor series expansion to get the weight term
                TCS([W_fuel[icruise2]/W_end[icruise2] >= te_exp_minus1(z_bre[izbre], nterm=3)]),

                #breguet range eqn
##                TCS([RngCruise[izbre] <= z_bre[izbre]*LD[icruise2]*V[icruise2]/(TSFCcr2[iclimb1])]),
                TCS([RngCruise[0] <= z_bre[0]*LD[4]*V[4]/(TSFCcr21)]),
                TCS([RngCruise[1] <= z_bre[1]*LD[5]*V[5]/(TSFCcr22)]),
                #time
                thours[icruise2]*V[icruise2]  == RngCruise[izbre],

##                TSFCcr22 == c1*units('1/hr'),
##                TSFCcr21 == c1*units('1/hr'),
                ])
        
        #constraint on the aircraft meeting the required range
        for i in range(min(izbre), max(izbre)+1):
            constraints.extend([
                TCS([RngCruise[i] == ReqRngCruise/(Ncruise2)])
                ])
            
        constraints.extend([
            #substitue these values later
            LD[icruise2]  == 10,
            M[icruise2] == 0.8
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
        cmc = CommericalMissionConstraints(0)
        climb1 = Climb1()
        climb2 = Climb2()
        cruise2 = Cruise2()
        atm = Atmosphere(Nseg)
        eonD = EngineOnDesign()
        eoffD = EngineOffDesign()
        
        for i in range(Nseg):
            None
        
        substitutions = {      
            'W_{e}': 44000*9.8*units('N'),
            'W_{payload}': 20000*9.8*units('N'),
            'V_{stall}': 120,
            '\\frac{L}{D}_{max}': 15,
            'ReqRng': 1000,
            'C_{d_0}': .05,
            'K': 0.10,
            'S': 124.58,
            'c1': 2,
            'h_{toc}': 30000,
            'thrust': 60000*units('lbf'),

##            'P_0': 19.8,  
##            'M_0': 0.8,
##            'T_{t_4}': 1400,
##            '\pi_f': 1.5,
##            '\pi_{lc}': 3,
##            '\pi_{hc}': 10,
##            '\pi_{d}': .99,
##            '\pi_{fn}': .98,
##            '\pi_{tn}': .99,
##            '\pi_{b}': .94,
##            'alpha': 8,
##            'alphap1': 9,
##            'M_{4a}': 1,    #choked turbines
##            'M_2': .4,
##            'M_{2.5}': .5,
##            'hold_{2}': 1+.5*(1.398-1)*.4**2,
##            'hold_{2.5}': 1+.5*(1.354-1)*.5**2,
##            'T_{ref}': 288.15,
##            'P_{ref}': 101.325,

            'G_f': 1,
            'N_{{bar}_Df}': 1,
            'F_{spec}': 8.0e+04 ,
            'T_{t_{4spec}}': 1300,
            }
        #for engine on design must link T0, P0, F_D,TSFC w/TSFC from icruise 2
        
        self.submodels = [cmc, climb1, climb2, cruise2, eonD, eoffD]

        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_{cr21}_Cruise2': 'TSFC_E_EngineOnDesign', 'TSFC_{cr22}_Cruise2': 'TSFC_E_EngineOffDesign'})

        lc = LinkedConstraintSet(constraints, exclude={'T_0', 'P_0', 'M_0', 'a_0', 'u_0', 'P_{t_0}', 'T_{t_0}', 'h_{t_0}', 'P_{t_1.8}',
                                                       'T_{t_1.8}', 'h_{t_1.8}', 'P_{t_2}', 'T_{t_2}', 'h_{t_2}', 'P_{t_2.1}','T_{t_2.1}'
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
                                                        ,'p_{tild_lc}','m_{{tild}_slc}','p_{{tild}_slc}'})

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
        return out

    
if __name__ == '__main__':
    m = CommercialAircraft()
    sol = m.localsolve(kktsolver = "ldl", verbosity = 4)
    
##    print sol('Drag')
##    print sol('thrust')
##    print sol('Drag')
##    print sol('thrust').to('N')
    
##    sol = m.determine_unbounded_variables(m,verbosity=4, iteration_limit=50)
    
##    print np.cumsum(sol('tmin'))
##    plt.plot(np.cumsum(sol('tmin')), sol('hft'))
##    plt.title('Altitude vs Time')f
##    plt.ylabel('Altitude [ft]')
##    plt.xlabel('Time [min]')
##    plt.show()
    
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
