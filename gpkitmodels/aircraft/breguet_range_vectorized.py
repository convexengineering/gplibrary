"""Breguet range model"""
from gpkit import Variable, VectorVariable, Model, units, SignomialEquality, SignomialsEnabled
from gpkit.constraints.tight import TightConstraintSet as TCS
from gpkit.tools import te_exp_minus1
import numpy as np

class BreguetRange(Model):
    """
    Breguet Range Model

    Assumptions
    -----------
    Fuel burn varies linearily with thrust.

    Arguments
    ----------
    n: the number of discretized segments, default value is 1
    """
    def __init__(self, n=1, **kwargs):

        #Aero parameters
        D = VectorVariable(n, 'D', 'N', 'Drag')
        Cd0 = Variable('Cd0', '-', 'Fuselage Profile Drag Coefficient')
        K = Variable('K', '-', 'Parametric Drag Model Parameter')
        LD = VectorVariable(n, 'L/D', '-', 'Lift to Drag Ratio')
        Cl = VectorVariable(n, 'Cl', '-', 'Segment Lift Coefficient')
        cdp = VectorVariable(n, 'c_{d_p}', '-', 'Airfoil Profile Drag for NC 130 (has transonic effects)')
        cdp_r = VectorVariable(n, 'c_{d_p_r}', '-', 'Transonic Drag Rise Hold Variable')

        #Aircraft parameters
        S = Variable('S', 'm^2', 'Wing Planform Area')
        
        #Atmosphere
        rho = Variable('rho', 'kg/m^3', 'Air Density')
        T = Variable('T', 'K', 'Air Temperature')
        gamma = Variable('\\gamma', 1.4, '-', 'Specific Heat Ratio for Air')
        R = Variable('R', 287, 'J/kg/K', 'Gas Constant for Air')
        a = Variable('a', 'kts', 'Speed of Sound')

        #Engine
        TSFC   = Variable('TSFC', '1/hr',
                          'Thrust specific fuel consumption')

        #velocities
        V = VectorVariable(n, 'V', 'knots', 'Segment Velocity')

        #Mach number
        M = VectorVariable(n, 'M', '-', 'Mach Number')

        #Weights
        W_start = VectorVariable(n, 'W_start', 'N', 'Segment Start Weight')
        W_end = VectorVariable(n, 'W_end', 'N', 'Segment End Weight')
        W_fuel = VectorVariable(n, 'W_fuel', 'N', 'Segment Fuel Weight')
        W_fuelTotal = Variable('W_fuelTotal', 'N', 'Segment Fuel Weight')
        W_e = Variable('W_{e}', 'N', 'Operating empty weight')
        W_avg = VectorVariable(n, 'W_avg', 'N', 'Average Segment Weight')

        #Constants
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        #Range Variables
        Rng = VectorVariable(n, 'Req', 'nautical_miles', 'Segment Range')
        ReqRng = Variable('ReqRng', 'nautical_miles', 'Total Required Range')
        
        #Breguet Paramter
        z_bre = VectorVariable(n, 'z_{bre}', '-', 'Breguet parameter')

        #Time
        t = VectorVariable(n, 't', 'hr', 'Segment Flight Time')

        #set up index list
        i = map(int, np.linspace(0, n - 1, n))

        with SignomialsEnabled():
            constraints = [
                #Segment weight constraints
                W_avg[i] == (W_start[i] * W_end[i])**.5,
                TCS([W_start[i] >= W_end[i] + W_fuel[i]]),

                #total weight constraint
                W_start[0] <= W_fuelTotal + W_e,

                #end weight constraint
                W_end[n-1] == W_e,

                #compute the speed of sound
                a == (gamma * R * T)**.5,

                #compute the transonic drag rise
                cdp_r >= (1.02458748e10 * Cl**15.587947404823325 * M**156.86410659495155 +
                2.85612227e-13 * Cl**1.2774976672501526 * M**6.2534328002723703 +
                2.08095341e-14 * Cl**0.8825277088649582 * M**0.0273667615730107 +
                1.94411925e+06 * Cl**5.6547413360261691 * M**146.51920742858428),

                cdp**6.5 >= cdp_r,

                #drag constraint
                TCS([D[i] >= (.5*S*rho*V[i]**2)*(Cd0 + K * Cl[i]**2 + cdp)]),

                #constraint on the lift coefficient, assumes steady level flight
                W_avg[i] == .5 * Cl[i] * S * rho * V[i]**2,

                #Breguet Range parameter constraints
                TCS([W_fuel[i]/W_end[i] >= te_exp_minus1(z_bre[i],3)]),
                TCS([z_bre[i] >= TSFC * t[i] * D[i]/ W_avg[i]]),

                #compute the Mach number
                M == V/a,

                #constraint on the segment time
                t[i] * V[i] == Rng[i],

                #total range is equal to the sum of the segment ranges
                TCS([ReqRng <= sum(Rng)]),

                #total fuel burn is equal to the sum of segment fuel burns
                TCS([W_fuelTotal <= sum(W_fuel)]),
                ]

        for j in range(1,n):
            constraints.extend([
                W_start[j] == W_end[j-1],
                ])

        #The objective is to minimze total fuel burn
        objective = W_fuelTotal  # Minimze total fuel burn

        #build the model
        Model.__init__(self, objective, constraints)

if __name__ == '__main__':
       m = BreguetRange(3)

       m.substitutions.update({
            'S': 125,                 #approx a B737 wing area
            'ReqRng': 1000,            #1,000 mile required range
            'W_{e}': 40000*9.81,        #approx empty weight of B737 in N
            'rho': .31,               #air density at 12,000m (~40,000')
            'T': 216,
            'Cd0': .008,               #setting profile drag coefficient
            'K': .015,                #setting parametric drag model coefficient
            'TSFC': 0.5,              #setting segment TSFC
        })

       sol = m.localsolve(solver='mosek',verbosity=4)
