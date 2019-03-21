"""Airfoil lift vs. drag vs. Mach model."""

from math import pi
from gpkit import Model, Variable, units
from gpkit.constraints.tight import Tight
from gpkit.constraints.loose import Loose

class DummyAirfoil(Model):
    """Model for aerodynamics of a airfoil.

       This is a dummy/prototype model. I made it by guessing values to
       replicate the c_d vs Mach behavior shown in figure 8.34 of Drela's book [1].

       A better version of this model would be to run GPfit on MSES data.

        References:
            [1] Drela, Mark. "Flight Vehicle Aerodynamics," MIT Press,
                2014.
        """
    def setup(self, flight_state):
        ### Declare Geometric Programming variables ###
        mach_perp = Variable('mach_perp', '', 'Mach number of velocity component perpendicular to leading edge')
        c_l = Variable('c_l', '', '2D Lift coefficient')
        c_dp = Variable('c_dp', '', '2D Drag coefficient due to pressure')
        c_df = Variable('c_df', 50e-4, '', '2D Drag coefficient due to friction')
        c_la_2d = Variable('c_la_2d', 'radian**-1', 'Airfoil 2d d C_L/ d alpha')

        # Variables internal to the model
        _c_dp0 = Variable('_c_dp0', '', '2D Drag coefficient due to pressure at low Mach number')
        _mach_d = Variable('_mach_d', '', 'Drag divergence Mach number')
        _mach_model_max = Variable('_mach_model_max', 0.8, '', 'Maximum Mach number at which this model is valid')
        _c_l_model_max = Variable('_c_l_model_max', 1.1, '', 'Maximum lift coefficient at which this model is valid')

        # Drag divergence exponent
        div_exp = 32

        ### Declare Geometric Programming constraints. ###
        constraints = [
            # Piecewise c_dp vs. mach model
            c_dp >= _c_dp0,
            c_dp >= _c_dp0 * (mach_perp / _mach_d)**div_exp,

            # Pressure drag vs. lift fit
            # This fit was manually tuned to qualitatively match the
            # left-hand plot in Drela's [1] figure 8.34.
            # There first two terms are a quadratic fit. The third term
            # has c_l raised to a high exponent to cause a rapid
            # increase in drag above c_l = 0.9.
            Tight([_c_dp0 >= 20.5e-4 * c_l**2 + 14.9e-4 + 1e-4 * (c_l/0.9)**25], reltol=1e-3),

            # Drag divergence Mach number vs. lift linear fit
            # Note: If the 2nd piece of the c_dp constraint is not tight,
            # i.e. at low mach, this constraint will not become tight and
            # _mach_d will take on some arbitrary value. 
            _mach_d + 0.0667 * c_l <= 0.783,

            # Mach number limit of model
            Loose([mach_perp <= _mach_model_max]),

            # c_l limit of model
            Loose([c_l <= _c_l_model_max]),

            # Lift slope vs. Mach number from Prandtl-Glauert
            # TODO should this beta be based on mach or mach_perp?
            c_la_2d == 0.95 * 2 * pi / flight_state['beta_pg'],
        ]
        constraints += [flight_state]
        return constraints
