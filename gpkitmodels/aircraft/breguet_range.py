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
    Fuel burn varies linearily with thrust and is independent of velocity.
    """
    def __init__(self, **kwargs):
        #Variable definitions
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumtion')

        #breguet range parameter
        z_bre = Variable('z_bre', '-', 'Breguet Range Parameter')

        #aero
        Cl = Variable('Cl', '-', 'Segment Lift Coefficient')
        D = Variable('D', 'N', 'Drag')
        Cd0 = Variable('Cd0', '-', 'Profile Drag Coefficient')
        K = Variable('K', '-', 'Parametric Drag Model Parameter')
        S = Variable('S', 'm^2', 'Wing Planform Area')

        #weights
        W_fuel = Variable('W_fuel', 'N', 'Segment Fuel Weight')
        W_end = Variable('W_end', 'N', 'Segment End Weight')
        W_avg = Variable('W_avg', 'N', 'Average Segment Weight')
        W_start = Variable('W_start', 'N', 'Segment Start Weight')

        #atmosphere
        rho = Variable('rho', 'kg/m^3', 'Air Density')

        #speed
        V = Variable('V', 'kts', 'Cruise Airspeed')
        
        #Range
        Range = Variable('Range', 'mi', 'Required Segment Range')

        #flight time
        t = Variable('t', 'hr', 'Segment Flight Time')
        
        constraints = []

        #write out all required constraints
        constraints.extend([
            #Breguet Range parameter constraints
            TCS([W_fuel/W_end >= te_exp_minus1(z_bre,3)]),
            TCS([z_bre >= TSFC * t * D/ W_end]),

            #constraint on the lift coefficient, assumes steady level flight
            W_avg == .5 * Cl * S * rho * V**2,

            #drag constraint
            TCS([D >= (.5*S*rho*V**2)*(Cd0 + K * Cl**2)]),

            #constrain the starting weight such that the aircraft is only losing weight due to fuel burn
            TCS([W_start >= W_end + W_fuel]),

            #average weight is the geometric mean of the start and end weights
            W_avg == (W_start * W_end)**.5,

            #constraint on the segment time
            t * V == Range,
            ])

        substitutions = []

        #build the model
        Model.__init__(self, W_fuel, constraints, substitutions)

if __name__ == '__main__':
       m = BreguetRange()

       m.substitutions.update({
            'S': 125,                 #approx a B737 wing area
            'Range': 1000,            #1,000 mile required range
            'W_end': 40000*9.81,      #approx empty weight of B737 in N
            'V': 420,                 #set the cruise speed
            'rho': .31,               #air density at 12,000m (~40,000')
            'Cd0': .02,               #setting profile drag coefficient
            'K': .015,                #setting parametric drag model coefficient
            'TSFC': 0.5,              #setting segment TSFC
            'W_limit': 100000*9.81,
        })

       sol = m.solve()
