"Models for atmospheric quantities"
import numpy as np
from gpkit import units
from gpkit import VectorVariable, Variable, Model, units, SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import TightConstraintSet as TCS
# pylint: disable=bad-whitespace

GAS_CONSTANT = 287  # [J/(kg*K)]
Lval = 0.0065   # [K/m]
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
        hu   = VectorVariable(N, 'hu', '-', 'Altitude w/out units')
        h = VectorVariable(N, 'hClimb1', 'm', 'Altitude [meters]')
        mu  = VectorVariable(N, '\\mu', 'kg/(m*s)', 'Dynamic viscosity')
        p   = VectorVariable(N, 'pClimb1', 'Pa', 'Pressure')
        pu   = VectorVariable(N, 'pu', '-', 'Pressure w/out Units')
        T   = VectorVariable(N, 'TClimb1', 'K', 'Temperature')
        rho = VectorVariable(N, '\rhoClimb1', 'kg/m^3', 'Air Density w/out Units')
        rhou = VectorVariable(N, '\rhou', '-', 'Air Density')
        R = Variable('R', GAS_CONSTANT, 'J/kg/K', 'Gas Constant for Air')
        m = Variable('m', 1, 'm', '1 Meter')
        kgm3 = Variable('kgm3', 1, 'kg/m^3', '1 kg/m^3')
        pa =Variable('pa', 1, 'Pa', '1 Pa')


        L    = Variable('L', Lval, 'K/m', 'Temperature lapse rate')
        p_0  = Variable('p_0', 101325, 'Pa', 'Pressure at sea level')
        T_0  = Variable('T_0', 288.15, 'K', 'Temperature at sea level')

        
        constraints = []
    
        constraints += [
            #constrain temperature with the ideal gas law
            T == p/(rho*R),

            #bad way of handling units problem
            hu*m == h,
            pu*pa ==p,
            rhou*kgm3 == rho,
            ]
        
        with SignomialsEnabled():
                for i in range(N):
                    constraints.extend([
                        #pressure constraint
##                        SignomialEquality(pu[i]**-.58, (0.000162 * hu[i]**-6.78e-06+ 1.54e-10 * hu[i]**-2+ 2.87e-05 * (hu[i])**0.00355+ 0.000812 *
##                                  (hu[i])**0.00022+ 0.000248 * (hu[i])**8.7e-05+ 1.41e-15 * (hu[i])**2.88 + 5.66e-08 * (hu[i])**1.06)),
##
##                        #density constraint
##                        SignomialEquality(rhou[i]**(-2.15), (0.149 * (hu[i])**-0.00179+ 4.74e-10 * (hu[i])**2.47+ 0.132 * (hu[i])**-0.00142 + 3.55e-28 *
##                                   (hu[i])**6.89+ 0.0882 * (hu[i])**-0.00141+ 0.151 * (hu[i])**-0.000969+ 0.126 * (hu[i])**-0.00141+ 0.000232 *
##                                   (hu[i])**0.933)),

                        rhou ==1,
                        pu==100000,

##                        pu[i]**-.58 >= (0.000162 * hu[i]**-6.78e-06+ 1.54e-10 * hu[i]**-2+ 2.87e-05 * (hu[i])**0.00355+ 0.000812 *
##                                   (hu[i])**0.00022+ 0.000248 * (hu[i])**8.7e-05+ 1.41e-15 * (hu[i])**2.88 + 5.66e-08 * (hu[i])**1.06),
##
##                        #density constraint
##                        rhou[i]**(-2.15) >= (0.149 * (hu[i])**-0.00179+ 4.74e-10 * (hu[i])**2.47+ 0.132 * (hu[i])**-0.00142 + 3.55e-28 *
##                                   (hu[i])**6.89+ 0.0882 * (hu[i])**-0.00141+ 0.151 * (hu[i])**-0.000969+ 0.126 * (hu[i])**-0.00141+ 0.000232 *
##                                   (hu[i])**0.933),
                        ])


##                p[index] == 0.195 * (h[index]*units('1/Pa'))**0.702,
                    

##                               TCS([pu**-.58 <= (0.000162 * hu**-6.78e-06+ 1.54e-10 * hu**-2+ 2.87e-05 * (hu)**0.00355+ 0.000812 *
##                          (hu)**0.00022+ 0.000248 * (hu)**8.7e-05+ 1.41e-15 * (hu)**2.88 + 5.66e-08 * (hu)**1.06)]),
##                p[index] == 0.195 * (h[index]*units('1/Pa'))**0.702,
##                #density constraint
##                TCS([rhou**(-2.15) >= (0.149 * (hu)**-0.00179+ 4.74e-10 * (hu)**2.47+ 0.132 * (hu)**-0.00142 + 3.55e-28 *
##                           (hu)**6.89+ 0.0882 * (hu)**-0.00141+ 0.151 * (hu)**-0.000969+ 0.126 * (hu)**-0.00141+ 0.000232 *
##                           (hu)**0.933)]),


                # Temperature lapse rate constraint
##                TCS([T_0 <= T + L*h]),
##                T >= 216.65*units('K'),
                
                # Pressure-altitude relation
##                (p/p_0)**(1/th) == T/T_0,

                # Ideal gas law
##                rho == p/(R*T),
                
            #tomorrow try a one term fit for rho
##            array_sec = [SignomialEquality(rhou[i]**(-2.15), (0.149 * (hu[i])**-0.00179+ 4.74e-10 * (hu[i])**2.47+ 0.132 * (hu[i])**-0.00142 + 3.55e-28 *
##                           (hu[i])**6.89+ 0.0882 * (hu[i])**-0.00141+ 0.151 * (hu[i])**-0.000969+ 0.126 * (hu[i])**-0.00141+ 0.000232 *
##                           (hu[i])**0.933)) for i in range(len(hu))]

##            constraints.extend(array_sec)
        
##        su = Sutherland(N)
##        lc = su.link(constraints)

        Model.__init__(self, None, constraints, **kwargs)

class Sutherland(Model):
    """
    Dynamic viscosity (mu) as a function of temperature

    References:
    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
        atmos/atmos.html
    http://www.cfd-online.com/Wiki/Sutherland's_law
    """
    def __init__(self, N, **kwargs):
        index = map(int, np.linspace(1, N-1, N))
        T_s = VectorVariable(N, 'T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = VectorVariable(N, 'C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')
        T   = VectorVariable(N, 'T', 'K', 'Temperature')
        t_plus_ts_approx[index] = (T[index] + T_s[index]).mono_approximation({T: 288.15,
                                                         T_s: T_s.value})

        objective = mu
        constraints = [t_plus_ts_approx * mu == C_1 * T**1.5]

        Model.__init__(self, objective, constraints, **kwargs)

    @classmethod
    def test(cls):
        m = cls(2)
        m.substitutions.update({"T": 288})
        m.solve()

if __name__ == "__main__":
    Sutherland.test()
