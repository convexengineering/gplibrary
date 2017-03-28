"Models for atmospheric quantities"
import numpy as np
from gpkit import units
from gpkit import VectorVariable, Variable, Model, units, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import TightConstraintSet as TCS
##from gpkit.nomials.nomial_math import SignomialEqualityTrust
# pylint: disable=bad-whitespace

GAS_CONSTANT = 287  # [J/(kg*K)]
Lval = 0.0098   # [K/m]
th = 9.81/(GAS_CONSTANT*Lval) # [-]

class Atmosphere(Model):
    """
    Density, pressure, temperature, and dynamic viscosity as a function of altitude based on
    standard atmosphere model for the Troposphere and Tropopause. Constraints, excluding dynamic
    viscosity and temperature, are developed from a gpfit and atmosphereic data from a standard atmosphere model.

    Temperature is derived from the ideal gas law and dynamic viscosity comes from Sutherland's equation
    in the Sutherland class below.

    Assumptions: Atmosphere transition at ~11km

    References:
    Anderson, Introduction to Flight
    https://en.wikipedia.org/wiki/Density_of_air#Altitude
    http://www.digitaldutch.com/atmoscalc/index.htm

    """
    def __init__(self, N, **kwargs):
        #indicies of the vector variables
        index = map(int, np.linspace(1, N-1, N))
        #define variables
        hu = VectorVariable(N, 'hu', '-', 'Altitude w/out units')
        h = VectorVariable(N, 'h', 'm', 'Altitude [meters]')
        mu = VectorVariable(N, '\\mu', 'kg/(m*s)', 'Dynamic viscosity')
        p = VectorVariable(N, 'p', 'Pa', 'Pressure')
        pu = VectorVariable(N, 'pu', '-', 'Pressure w/out Units')
        T = VectorVariable(N, 'T', 'K', 'Temperature')
        rho = VectorVariable(N, '\rho', 'kg/m^3', 'Air Density w/out Units')
        rhou = VectorVariable(N, '\rhou', '-', 'Air Density')
        R = Variable('R', GAS_CONSTANT, 'J/kg/K', 'Gas Constant for Air')
        m = Variable('m', 1, 'm', '1 Meter')
        kgm3 = Variable('kgm3', 1, 'kg/m^3', '1 kg/m^3')
        pa = Variable('pa', 1, 'Pa', '1 Pa')


        L    = Variable('L', Lval, 'K/m', 'Temperature lapse rate')
        p_0  = Variable('p_0', 101325, 'Pa', 'Pressure at sea level')
        T_0  = Variable('T_0', 288.15, 'K', 'Temperature at sea level')

        
        constraints = []
    
        constraints += [
            #constrain temperature with the ideal gas law
##            T == p/(rho*R),

            #bad way of handling units problem
##            hu*m == h,
##            pu*pa ==p,
##            rhou*kgm3 == rho,

            h[0] == 1104.9*units('m'),
            h[1] == 2400.3*units('m'),
            h[2] == 4572.0*units('m'),
            h[3] == 7620.0*units('m'),
            h[4] == 9144.0*units('m'),
            h[5] == 9144.0*units('m'),

            rho == p/(R*T),
            p == 100000*units('Pa'),
##            T == 270*units('K')
            ]
 
        
        with SignomialsEnabled():
                for i in range(N):
                    constraints.extend([
                        #pressure constraint
##                        SignomialEquality(pu[i]**(-.58), (0.000162 * hu[i]**-6.78e-06+ 1.54e-10 * hu[i]**-2+ 2.87e-05 * (hu[i])**0.00355+ 0.000812 *
##                                  (hu[i])**0.00022+ 0.000248 * (hu[i])**8.7e-05+ 1.41e-15 * (hu[i])**2.88 + 5.66e-08 * (hu[i])**1.06)),

##                        SignomialEquality(pu[i]**(-.0401), 5.59e-07 * (hu[i])**1.19 + 0.629 * (hu[i])**0.000386),

##                        pu[i] **-.58 >= (0.000162 * hu[i]**-6.78e-06+ 1.54e-10 * hu[i]**-2+ 2.87e-05 * (hu[i])**0.00355+ 0.000812 *
##                                  (hu[i])**0.00022+ 0.000248 * (hu[i])**8.7e-05+ 1.41e-15 * (hu[i])**2.88 + 5.66e-08 * (hu[i])**1.06),

                        #density constraint
##                        SignomialEquality(rhou[i]**(-2.15), (0.149 * (hu[i])**-0.00179+ 4.74e-10 * (hu[i])**2.47+ 0.132 * (hu[i])**-0.00142 + 3.55e-28 *
##                                   (hu[i])**6.89+ 0.0882 * (hu[i])**-0.00141+ 0.151 * (hu[i])**-0.000969+ 0.126 * (hu[i])**-0.00141+ 0.000232 *
##                                   (hu[i])**0.933)),

##                        rhou[i]**(-2.15) >= (0.149 * (hu[i])**-0.00179+ 4.74e-10 * (hu[i])**2.47+ 0.132 * (hu[i])**-0.00142 + 3.55e-28 *
##                                   (hu[i])**6.89+ 0.0882 * (hu[i])**-0.00141+ 0.151 * (hu[i])**-0.000969+ 0.126 * (hu[i])**-0.00141+ 0.000232 *
##                                   (hu[i])**0.933)
                       SignomialEquality(T[i] + L*h[i], T_0),
##                       T_0 >= T[i] + L*h[i],
##                       (p[i]/p_0) == (T[i]/T_0)**th,
##                        rhou == 1,
##                        pu == 10000,

                        ])

        Model.__init__(self, sum(rho), constraints, **kwargs)

if __name__ == '__main__':
    m = Atmosphere(6)
    sol = m.localsolve(solver="mosek", verbosity = 4, iteration_limit=100, skipsweepfailures=True)
