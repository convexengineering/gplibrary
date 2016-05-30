"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
import numpy as np
from gpkit import Variable, Model, units

class CommericalAircraft(Model):
    def __init__(self):
        #define the number of discretizations for each flight segment
        Ntakeoff = 0
        Nclimb1 = 0
        Ntrans = 0
        Nclimb2 = 0
        Ncruise1 = 0
        Ncruiseclimb = 0
        Ncruise2 = 50
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
##        itakeoff = np.linspace(0, Ntakeoff - 1, Ntakeoff)
##        iclimb1 = np.linspace(Ntakeoff, itakeoff[len(itakeoff)-1]+Nclimb1, Nclimb1)
##        itrans = np.linspace(iclimb1[len(iclimb1)-1] + 1, Ntrans + iclimb1[len(iclimb1)-1], Ntrans)
##        iclimb2 = np.linspace(itrans[len(itrans)-1] + 1, Nclimb2 + itrans[len(itrans)-1], Nclimb2)
##        icruise1 = np.linspace(iclimb2[len(iclimb2)-1] + 1, Ncruise1 + iclimb2[len(iclimb2)-1], Ncruise1)
##        icruiseclimb = np.linspace(icruise1[len(icruise1)-1] + 1, Ncruiseclimb + icruise1[len(icruise1)-1], Ncruiseclimb)
##        icruise2 = np.linspace(icruiseclimb[len(icruiseclimb)-1] + 1, Ncruise2 + icruiseclimb[len(icruiseclimb)-1], Ncruise2)
##        idecent = np.linspace(icruise2[len(icruise2)-1] + 1, Ndecent + icruise2[len(icruise2)-1], Ndecent)
##        ilanding = np.linspace(idecent[len(idecent)-1] + 1, Nlanding + idecent[len(idecent)-1], Nlanding)
##        iresclimb = np.linspace(ilanding[len(ilanding)-1] + 1, Nresclimb + ilanding[len(ilanding)-1], Nresclimb)
##        iresdecent = np.linspace(iresclimb[len(iresclimb)-1] + 1, Nresdecent + iresclimb[len(iresclimb)-1], Nresdecent)
##        ireslanding = np.linspace(iresdecent[len(iresdecent)-1] + 1, Nreslanding + iresdecent[len(iresdecent)-1], Nreslanding)
##        ireshold = np.linspace(ireslanding[len(ireslanding)-1] + 1, Nreshold + ireslanding[len(ireslanding)-1], Nresdecent)

        #partial flight profile for use in development
        icruise2 = np.linspace(0, Ncruise2 -1, Ncruise2)

        #---------------------------------------
        #Variable definitions
        
        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #altitude
        hft = VectorVariable(Nseg, 'hft', 'feet', 'Altitude')
        h = VectorVariable(Nseg, 'h', 'm', 'Altitude')

        #Range
        Rng = VectorVariable(Nseg, 'range', 'natucical_miles', 'Range')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Required Mission Range')

        #aircraft weights
        W_e = Variable('W_{e}', 'lbf', 'Empty Weight of Aircraft')
        W_f = VectorVariable(Nseg, 'W_{f}', 'lbf', 'Weight of Fuel in Aircraft')
        W_p = Variable('W_{p}', 'lbf', 'Aircraft Payload Weight')
        W_mf = Variable('W_{mf}', 'lbf', 'Aricraft Weight Minus Fuel Weight')

        #aero
        LD = VectorVariable(Nseg, '\\frac{L}{D}', '-', 'Lift to Drag')
        LDmax = Variable('\\frac{L}{D}_{max}', 15, 'Maximum Lift to Drag')

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
        
        constraints = []

        #---------------------------------------
        #basic constraint definitions
        constraints.extend([
            #constraint on the aircraft weight minus fuel
            W_mf >= W_e + W_p
            #constraint on the aircraft meeting the required range
            Rng[] >= ReqRng
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
            RC[] = 0,
            theta[] = 0,
            Rng[] <= ReqRng,
            #taylor series expansion to get the weight term
            W_f[]/W_f[] >=  te_exp_minus1(z_bre, nterm=3),
            #breguet range eqn
            R <= z_bre*LD*V/(TSFC*g),
            Rng[] + R <= Rng[],
            
            )]

        #---------------------------------------
        #decent

        #----------------------------------------
        #landing


        #objective is to minimize the required fuel weight
        objective = 1/W_f

        m = Model(objective, constraints)

        sol = m.solve()

if __name__ == '__main__':
    CommericalAircraft()
