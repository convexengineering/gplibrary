"""Breguet range model"""
from gpkit import Variable, VectorVariable, Model, units, SignomialEquality, SignomialsEnabled
from gpkit.constraints.tight import TightConstraintSet as TCS
from gpkit.tools import te_exp_minus1
import numpy as np


class BreguetRange(Model):
    """Breguet Range Model

    Assumptions
    -----------
    Fuel burn varies linearily with thrust and is independent of velocity.

    arguemnt n to init is the number of discretized segments, default value is 1
    """
    def __init__(self, n=1, **kwargs):

        #Aero parameters
        D = VectorVariable(n, 'D', 'N', 'Drag')
        Cd0 = Variable('Cd0', '-', 'Profile Drag Coefficient')
        K = Variable('K', '-', 'Parametric Drag Model Parameter')
        LD = VectorVariable(n, 'L/D', '-', 'Lift to Drag Ratio')
        Cl = VectorVariable(n, 'Cl', '-', 'Segment Lift Coefficient')

        #Aircraft parameters
        S = Variable('S', 'm^2', 'Wing Planform Area')
        
        #Atmosphere
        rho = Variable('rho', 'kg/m^3', 'Air Density')

        #Engine
        TSFC   = Variable('TSFC', '1/hr',
                          'Thrust specific fuel consumption')

        #velocities
        V = VectorVariable(n, 'V', 'knots', 'Segment Velocity')
        Vmax  = Variable('V_{max}', 'knots', 'Maximum velocity')

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
        print i

        with SignomialsEnabled():
            constraints = [
                        #Segment weight constraints
                        W_avg[i] == (W_start[i] * W_end[i])**.5,
                        TCS([W_start[i] >= W_end[i] + W_fuel[i]]),

                        #total weight constraint
                        SignomialEquality(W_start[0] , W_fuelTotal + W_e),

                        #end weight constraint
                        W_end[n-1] == W_e,

                        #drag constraint
                        TCS([D[i] >= (.5*S*rho*V[i]**2)*(Cd0 + K * Cl[i]**2)]),

                        #constraint on the lift coefficient, assumes steady level flight
                        W_avg[i] == .5 * Cl[i] * S * rho * V[i]**2,

                        #Breguet Range parameter constraints
                        TCS([W_fuel[i]/W_end[i] >= te_exp_minus1(z_bre[i],3)]),
                        TCS([z_bre[i] >= TSFC * t[i] * D[i]/ W_avg[i]]),
                        W_e==40000*9.8*units('N'),
                        D[i]==10000*units('N'),
                        TSFC == .5*units('1/hr'),
                        t[i]==.5*units('hr'),
                        W_start[i] <= 1000000*units('N'),
                        W_avg[i] == 80000*units('N'),
                        #velocity constraint
                        V[i] <= Vmax,

                        #constraint on the segment time
                        t[i] * V[i] == Rng[i],

                        #total range is equal to the sum of the segment ranges
                        TCS([ReqRng == sum(Rng)]),

                        #total fuel burn is equal to the sum of segment fuel burns
                        TCS([W_fuelTotal >= sum(W_fuel)]),
                          ]

        #The objective is to minimze total fuel burn
        objective = W_fuelTotal  # Minimze total fuel burn

        #substitutions
        substitutions = []

        #build the model
        Model.__init__(self, objective, constraints, substitutions)
        
    def test(self, BR):
        BR.substitutions.update({
            'S': 125,                 #approx a B737 wing area
            'ReqRng': 500,            #1,000 mile required range
            'W_e': 40000*9.81,        #approx empty weight of B737 in N
            'rho': .31,               #air density at 12,000m (~40,000')
            'Cd0': .02,               #setting profile drag coefficient
            'K': .015,                #setting parametric drag model coefficient
            'TSFC': 0.5,              #setting segment TSFC
            'V_{max}': 420,           #set the max velocity limit
        })

        BR.solve(solver = "mosek", verbosity = 4)

if __name__ == '__main__':
       m = BreguetRange()
##    m.test(m)
       m.substitutions.update({
            'S': 125,                 #approx a B737 wing area
            'ReqRng': 500,            #1,000 mile required range
            'W_e': 40000*9.81,        #approx empty weight of B737 in N
            'rho': .31,               #air density at 12,000m (~40,000')
            'Cd0': .02,               #setting profile drag coefficient
            'K': .015,                #setting parametric drag model coefficient
            'TSFC': 0.5,              #setting segment TSFC
            'V_{max}': 420,           #set the max velocity limit
        })

       sol = m.solve(solver='mosek',verbosity=4)
       print sol.table()
