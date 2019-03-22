"""Wing aerodynamcis models for small, fast aircraft."""

from math import pi
from gpkit import Model, Variable, units
from gpkit import SignomialsEnabled, SignomialEquality
from gpkit.constraints.tight import Tight
from gpkit.constraints.loose import Loose

class WingAeroSubsonic(Model):
    def setup(self, flight_state, airfoil, limit_sweep_to_valid_range=True,
              pressure_drag_cos_exponent=1):
        """Model for aerodynamics of a wing.

        Valid at speeds up to the drag diverence Mach number,
        or maybe to Mach 1 for a swept wing.

        The lift model is based on 4.1.3.2 of the DATCOM [1]
        and Eq. 12.6 of Raymer [2].

        Beware: It seems GPkit's autosweep does not work with this model
         - Matt, 2019-03-06 14:40.

        Arguments:
            flight_state (FlightState): A gpkit Model describing the
                static pressure, Mach number, etc. at which the wing is flying.
            limit_sweep_to_valid_range (boolean): If true, include an upper limit
                constraint on wing sweep, to prevent sweep from going to high
                values at which the model is not valid. The model uses Taylor series
                approximations of trigonometric functions of sweep, and these
                are not accurate for very high sweep angles.
            pressure_drag_cos_exponent (scalar): The exponent on the `cos(sweep)`
                term multiplying the pressure drag coefficient. Drela argues that
                this exponent should be 3 (See [3], eqn. 8.183). However, Hoerner [4]
                argues that it should be 1.

        References:
            [1] USAF Stability and Control DATCOM.
            [2] Raymer, Daniel. "Aircraft Design: a Conceptual Approach,"
                5th edition, AIAA, 2012.
            [3] Drela, Mark. "Flight Vehicle Aerodynamics," MIT Press,
                2014.
            [4] Hoerner, S. F. "Fluid Dynamic Drag,"
                1965.
        """
        """
        TODO: add the effects of taper & wing load distribution.
        Perhaps this will be a useful source:
        "Estimating the Oswald Factor From Basic Aircraft Geometrical Parameters"
        http://www.fzt.haw-hamburg.de/pers/Scholz/OPerA/OPerA_PRE_DLRK_12-09-10_MethodOnly.pdf
        """

        #pylint: disable=invalid-name
        self.pressure_drag_cos_exponent = pressure_drag_cos_exponent

        ### Declare Geometric Programming variables ###
        lift = Variable('lift', 'newton', 'Wing lift force')
        drag = Variable('drag', 'newton', 'Wing drag force')

        # Wing dimensions
        S_ref = Variable('S_ref', 'm^2', 'Wing reference area')
        span = Variable('span', 'm', 'Wing span')
        aspect = Variable('aspect', '', 'Wing aspect ratio')
        sweep_c2 = Variable('sweep_c2', '', 'Wing sweep angle at half chord')

        # Flight conditions
        mach_perp = airfoil['mach_perp']
        alpha = Variable('alpha', 'radian', 'Wing angle of attack')

        # Lift coefficients
        C_L = Variable('C_L', '', 'Wing lift coefficient')
        C_La = Variable('C_La', 'radian**-1', 'Wing d C_L/ d alpha')
        alpha_max = Variable('alpha_max', 12, 'degree', 'Max angle of attack for linear regime.')

        # Drag coefficients
        C_D = Variable('C_D', '', 'Wing drag coefficient')
        C_Di = Variable('C_Di', '', 'Wing induced drag coefficient')
        span_eff = Variable('span_eff', 0.9, '', 'Span efficiency factor, aka e')

        # Helper variables
        _t1 = Variable('_t1', '', 'Helper variable for posynomial in sqrt')
        _t2 = Variable('_t2', '', 'Helper variable for tan^2 Taylor series')
        _helper_cos = Variable('_helper_cos', '', 'Helper variable for cos Taylor series')


        ### Declare Geometric Programming constraints. ###
        # Geometric constraints
        constraints_geom = [
            aspect == span**2 / S_ref,
            # Bound sweep from below to avoid dual infeasibility as
            # sweep --> 0.
            sweep_c2 >= 1e-9,
            # Sweep physically can't be more than 90 degrees.
            sweep_c2 <= pi / 2,
        ]
        if limit_sweep_to_valid_range:
            # Our Taylor series approximation for tan^2(sweep)
            # is inaccurate above about 70 degrees.
            constraints_geom += Loose([sweep_c2 <= 70 * pi / 180])

        # Lift model constraints
        constraints_lift = [
            # Definition of C_L
            lift == C_L * S_ref * flight_state['p_dyn'],

            # Definition of C_la
            C_L == C_La * alpha,
            airfoil['c_l'] == airfoil['c_la_2d'] * alpha,

            # DATCOM lift model
            # Re-arranged for GP, using helper variables
            Tight([C_La**-1 >= (1 + _t1**0.5) / (pi * aspect)]),

            # DATCOM lift model
            # stuff in the square root
            Tight([_t1 >= 1 + (aspect * pi / (airfoil['c_la_2d'] / units('radian')))**2
                   * (1 + _t2 / flight_state['beta_pg']**2)],
                   reltol=1e-3),

            # DATCOM lift model
            # Taylor series approximation for tan**2(sweep)
            Tight([_t2 >= sweep_c2**2
                   + 2/3. * sweep_c2**4
                   + 17/45. * sweep_c2**6
                   + 62/315. * sweep_c2**8
                   + 1382/14175. * sweep_c2**10
                   + 21844/467775. * sweep_c2**12],
                   reltol=1e-3),

            # Limit of angle of attack range of linear model
            alpha <= alpha_max,
        ]

        constriants_drag = [
            # M_\perp = M_\infty \cos(\Lambda)
            mach_perp == flight_state['mach'] * _helper_cos,

            # Induced drag
            C_Di == C_L**2 / (pi * aspect * span_eff),

            # Drag buildup for the wing
            # TODO Drela suggests that the cosine term should be cos**3,
            # but Hoerner suggests that it should be cos.
            # We need to validate this model against AVL or wind tunnel data
            C_D >= (airfoil['c_df'] + (airfoil['c_dp'] * _helper_cos**self.pressure_drag_cos_exponent)
                + C_Di),

            # Definition of C_D
            drag == flight_state['p_dyn'] * S_ref * C_D,
        ]

        with SignomialsEnabled():
            constriants_drag += [
                # Taylor series for cos(\Lambda)
                Tight([_helper_cos >= 1. - sweep_c2**2 / 2.
                       + sweep_c2**4 / 24.
                       - sweep_c2**6 / 720.
                       + sweep_c2**8 / 40320.
                       - sweep_c2**10 / 3628800.
                       + sweep_c2**12 / 479001600.
                      ], reltol=1e-3),
            ]

        constraints = constraints_geom + constraints_lift + constriants_drag
        constraints += [airfoil]
        return constraints
