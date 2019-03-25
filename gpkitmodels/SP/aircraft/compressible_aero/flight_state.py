"""Flight state model."""

from gpkit import Model, Variable, units
from gpkit import SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight
from gpkit.constraints.loose import Loose

EPSILON = 1e-9    # a small number

class FlightState(Model):
    def setup(self):
        ### Declare Geometric Programming variables ###
        mach = Variable('mach', '', 'Freestream Mach number')
        p_dyn = Variable('p_dyn', 'pascal', 'Freestream dynamic pressure')
        p_static = Variable('p_static', 30e3, 'pascal', 'Freestream static pressure')
        gamma_air = Variable('gamma_air', 1.40, '', 'Ratio of specific heats of air')
        beta_pg = Variable('beta_pg', '', 'Prandtl-Glauert beta parameter')

        ### Declare Geometric Programming constraints. ###
        # Dynamic pressure
        # Source: https://en.wikipedia.org/wiki/Dynamic_pressure#Compressible_flow
        constraints = [
            p_dyn == 0.5 * gamma_air * p_static * mach**2,

            # Limit on Mach number range of model
            Loose([mach <= 1]),
            Loose([beta_pg >= EPSILON]),
            Loose([beta_pg <= 1]),
        ]

        with SignomialsEnabled():
            constraints += [
                # Prandtl-Glauert parameter for subsonic flow.
                # Instead of using a SignomialEquality, use two 'opposing'
                # inequalities, so that the model can simplify to a GP
                # model when Mach is specified (which it often is).
                Tight([beta_pg**2 >= 1 - mach**2]),
                Tight([-beta_pg**2 >= -1 + mach**2]),
            ]
        return constraints
