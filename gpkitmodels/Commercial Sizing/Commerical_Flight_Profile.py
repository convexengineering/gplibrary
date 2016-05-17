"""Standard tube and wing commercial aircraft sizing"""
from numpy import pi
from gpkit import Variable, Model, Units

class CommericalAircraft(Model):
    def __init__(self):
        #define the number of discretizations for each flight segment
        Ntakeoff = 5
        Nclimb1 = 10
        Ntrans = 5
        Nclimb2 = 10
        Ncruise1 = 50
        Ncruiseclimb = 10
        Ncruise2 = 50
        Ndecent = 5
        Nlanding = 5
        Nresclimb = 20
        Nresdecent = 5
        Nreslanding = 5
        Nreshold = 10

        Nclimb = Nclimb1 + Nclimb2 + Ncruiseclimb
        Ncruise = Ncruise1 + Ncruise2
        Nres = Nresclimb + Nresdecent + Nreslanding + Nreshold

        Nseg = Nclimb + Ncruise + Nres + Ntakeoff + Ntrans + Ndecent + Nlanding

        #---------------------------------------
        #Variable definitions
        
        #physical constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational Acceleration')
        gamma = Variable('\gamma', 1.4, '-', 'Air Specific Heat Ratio')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')

        #aircraft weights
        W_e = Variable('W_{e}', 'lbf', 'Empty Weight of Aircraft')
        W_f = VectorVariable(Nseg, 'W_{f}', 'lbf', 'Weight of Fuel in Aircraft')
        W_p = Variable('W_{p}', 'lbf', 'Aircraft Payload Weight')

        #velocitites and mach numbers
        V = VectorVariable(Nseg, 'V', 'knots', 'Aircraft Flight Speed')
        M = VectorVariable(Nseg, 'M', '-', 'Aircraft Mach Number')
        

        constraints = []

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

        #---------------------------------------
        #decent

        #----------------------------------------
        #landing
