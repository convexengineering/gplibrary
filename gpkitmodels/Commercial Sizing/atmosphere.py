"Models for atmospheric quantities"
import numpy as np
from gpkit import units
from gpkit import Variable, Model, units, SignomialsEnabled
from gpkit.constraints.tight import TightConstraintSet as TCS
# pylint: disable=bad-whitespace

GAS_CONSTANT = 287  # [J/(kg*K)]

hu   = Variable('hu', '-', 'Altitude w/out units')
h = Variable('h', 'm', 'Altitude [meters]')
mu  = Variable('\\mu', 'kg/(m*s)', 'Dynamic viscosity')
p   = Variable('p', 'Pa', 'Pressure')
T   = Variable('T', 'K', 'Temperature')
rho = Variable('\rho', 'kg/m^3', 'Air Density')
R = Variable('R', GAS_CONSTANT, 'J/kg/K', 'Gas Constant for Air')


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
    def __init__(self, **kwargs):
        constraints = []

        constraints += [
            #constrain temperature with the ideal gas law
            T == p/(rho*R),

            #bad way of handling units problem
            hu*units('m') == h,
            
            #pressure constraint
            TCS([p**-.58*units('Pa^.58') >= (0.000162 * hu**-6.78e-06+ 1.54e-10 * hu**-2+ 2.87e-05 * (hu)**0.00355+ 0.000812 *
                      (hu)**0.00022+ 0.000248 * (hu)**8.7e-05+ 1.41e-15 * (hu)**2.88 + 5.66e-08 * (hu)**1.06)]),
            
            #density constraint
            TCS([rho**-2.15*units('(kg/m^3)^2.15') >= (0.149 * (hu)**-0.00179+ 4.74e-10 * (hu)**2.47+ 0.132 * (hu)**-0.00142 + 3.55e-28 *
                       (hu)**6.89+ 0.0882 * (hu)**-0.00141+ 0.151 * (hu)**-0.000969+ 0.126 * (hu)**-0.00141+ 0.000232 *
                       (hu)**0.933)]),
                        ]
        su = Sutherland()
        lc = su.link(constraints)

        Model.__init__(self, None, lc, **kwargs)

class Sutherland(Model):
    """
    Dynamic viscosity (mu) as a function of temperature

    References:
    http://www-mdp.eng.cam.ac.uk/web/library/enginfo/aerothermal_dvd_only/aero/
        atmos/atmos.html
    http://www.cfd-online.com/Wiki/Sutherland's_law
    """
    def __init__(self, **kwargs):

        T_s = Variable('T_s', 110.4, "K", "Sutherland Temperature")
        C_1 = Variable('C_1', 1.458E-6, "kg/(m*s*K^0.5)",
                       'Sutherland coefficient')

        t_plus_ts_approx = (T + T_s).mono_approximation({T: 288.15,
                                                         T_s: T_s.value})

        objective = mu
        constraints = [t_plus_ts_approx * mu == C_1 * T**1.5]

        Model.__init__(self, objective, constraints, **kwargs)

    @classmethod
    def test(cls):
        m = cls()
        m.substitutions.update({"T": 288})
        m.solve()

if __name__ == "__main__":
    Sutherland.test()
