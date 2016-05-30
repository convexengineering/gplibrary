"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import numpy as np
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1

#TODO
#make discretization lists for the breguet parmeter

class CommericalAircraft(Model):
    def __init__(self):
        #define the number of discretizations for each flight segment
        Ntakeoff = 0
        Nclimb1 = 0
        Ntrans = 0
        Nclimb2 = 0
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
        icruise2 = map(int, np.linspace(0, Ncruise2 - 1, Ncruise2))
        
        #---------------------------------------
        #Variable definitions
        
        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #time
        t = VectorVariable(Nseg, 't', 'min', 'flight time in minutes')

        #altitude
        hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude')
        h = VectorVariable(Nseg, 'h', 'm', 'Altitude')

        #Range

        #HOW TO HANDLE NM
        
        Rng = VectorVariable(Nseg, 'range', 'miles', 'Segment Range')
        ReqRng = Variable('ReqRng', 'miles', 'Required Mission Range')

        #aircraft weights
        W_e = Variable('W_{e}', 'lbf', 'Empty Weight of Aircraft')
        W_start = VectorVariable(Nseg, 'W_{f}', 'lbf', 'Segment Start Weight')
        W_ftotal = Variable('W_{f_{total}}', 'lbf', 'Total Fuel Weight')
        W_end = VectorVariable(Nseg, 'W_{end}', 'lbf', 'Segment End Weight')
        W_payload = Variable('W_{payload}', 'lbf', 'Aircraft Payload Weight')
        W_total = Variable('W_{mf}', 'lbf', 'Total Aircraft Weight')

        #aero
        LD = VectorVariable(Nseg, '\\frac{L}{D}', '-', 'Lift to Drag')
        LDmax = Variable('\\frac{L}{D}_{max}', 15, '-', 'Maximum Lift to Drag')

        #CHECK VALUE OF CD0, k, AND K UNITS
        
        Cd0 = Variable('C_{d_0}', .012, '-', 'Aircraft Cd0')
        K = Variable('K', 1, '-', 'K for Parametric Drag Model')

        #aircraft geometry

        #CHECK THE VALUE OF S
        
        S = Variable('S', 120, 'm^2', 'Wing Planform Area')

        #velocitites and mach numbers
        V = VectorVariable(Nseg, 'V', 'knots', 'Aircraft Flight Speed')
        M = VectorVariable(Nseg, 'M', '-', 'Aircraft Mach Number')

        #Climb/Decent Rate
        RC = VectorVariable(Nseg, 'RC', 'feet/s', 'Rate of Climb/Decent')
        theta = VectorVariable(Nseg, '\\theta', '-', 'Aircraft Climb Angle')

        #breguet parameter
        z_bre = VectorVariable(Ncruise1+Ncruise2, 'z_{bre}', '-', 'Breguet Parameter')

        #engine
        TSFC = VectorVariable(Nseg, 'TSFC', 0.6/3600, 'lbm/hr/lbf', 'Thrust Specific Fuel Consumption')
        
        constraints = []

        #---------------------------------------
        #basic constraint definitions
        constraints.extend([
            #constraint on the aircraft meeting the required range
            Rng[Nseg-1] >= ReqRng,
            #convert m to ft
            hft == h,
            #constraints on the various weights
            W_e + W_payload + W_ftotal <= W_total,
            W_start[0] == W_total,
            W_e + W_payload <= W_end[Nseg-1] 
            ])

        #---------------------------------------
        #takeoff

        #---------------------------------------
        #Altitude Constraints

        #--------------------------------------
        #Climb #1 (sub 10K subject to 250KTS speed limit)

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

        constraints.extend([
            #constrain the climb rate by holding altitude constant
            hft[icruise2] == 35000,
            #taylor series expansion to get the weight term
            W_start[icruise2]/W_end[icruise2] >=  te_exp_minus1(z_bre[icruise2], nterm=3),
            #breguet range eqn
            Rng[icruise2] <= z_bre[icruise2]*LD*V/(TSFC*g)
            ])

        #---------------------------------------
        #decent

        #----------------------------------------
        #landing


        #objective is to minimize the total required fuel weight
        objective = 1/W_ftotal

        m = Model(objective, constraints)

        sol = m.solve()

if __name__ == '__main__':
    CommericalAircraft()
