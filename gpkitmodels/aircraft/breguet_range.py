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
    def __init__(self, **kwargs):
        #Variable definitions
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumtion')
        z_bre = Variable('z_bre', '-', 'Breguet Range Parameter')
        D = Variable('D', 'N', 'Drag')
        rho = Variable('rho', 'kg/m^3', 'Air Density')
        V = Variable('V', 'kts', 'Cruise Airspeed')
        Cd0 = Variable('Cd0', '-', 'Profile Drag Coefficient')
        K = Variable('K', '-', 'Parametric Drag Model Parameter')
        W_fuel = Variable('W_fuel', 'N', 'Segment Fuel Weight')
        W_end = Variable('W_end', 'N', 'Segment End Weight')
        W_avg = Variable('W_avg', 'N', 'Average Segment Weight')
        W_start = Variable('W_start', 'N', 'Segment Start Weight')
        Range = Variable('Range', 'mi', 'Required Segment Range')
        t = Variable('t', 'hr', 'Segment Flight Time')
        Cl = Variable('Cl', '-', 'Segment Lift Coefficient')
        S = Variable('S', 'm^2', 'Wing Planform Area')
        LD = Variable('L/D', '-', 'Lift to Drag Ratio')
        W_limit = Variable('W_limit', 'N', 'Non-Physical Weight Limit Needed for Convergence')
        
        
        constraints = []
        
        #write out all required constraints
        constraints.extend([
                #Breguet Range parameter constraints
                TCS([W_fuel/W_end >= te_exp_minus1(z_bre,3)]),
                TCS([z_bre >= TSFC * t * D/ W_avg]),
                
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
                
                #non-physical weight limit for demonstration purposes, relevant only in the second Breguet formulation
                W_start <= W_limit,
            ])
        
        substitutions = []
        
        #build the model
        Model.__init__(self, W_fuel, constraints, substitutions)

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
            'Range': 1000,            #1,000 mile required range
            'W_end': 40000*9.81,      #approx empty weight of B737 in N
            'V': 420,                 #set the cruise speed
            'rho': .31,               #air density at 12,000m (~40,000')
            'Cd0': .02,               #setting profile drag coefficient
            'K': .015,                #setting parametric drag model coefficient
            'TSFC': 0.5,              #setting segment TSFC
            'W_limit': 100000*9.81,
        })

       sol = m.solve(solver='mosek',verbosity=4)
       print sol.table()
