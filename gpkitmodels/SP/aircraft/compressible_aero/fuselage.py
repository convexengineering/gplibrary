"""Fuselage aerodynamics model for small, fast aircraft."""

from math import pi
from gpkit import Model, Variable, parse_variables
from gpkit.constraints.tight import Tight

class FuselageGeom(Model):
    def setup(self):
        ### Declare Geometric Programming variables ###
        # Fuselage dimensions
        length = Variable('length', 'meter', 'Fuselage length')
        diam = Variable('diam', 'meter', 'Fuselage diameter')
        fineness = Variable('fineness', '', 'Fuselage fineness ratio (length/diameter)')
        S_ref = Variable('S_ref', 'm^2', 'Fuselage cross section reference area')

        ### Declare Geometric Programming constraints. ###
       # Geometric constraints
        constraints_geom = {
            'definition of fineness ratio':
                fineness == length / diam,
            'area of a circle':
                S_ref == (pi / 4) * diam**2,
        }
        constraints = constraints_geom
        return constraints

    def dynamic(self, mach_regime, *args):
        """Create a dynamic/performance model for the fuselage.

        See https://gpkit.readthedocs.io/en/latest/modelbuilding.html#multipoint-analysis-modeling
        for a description of the dynamic model paradigm.

        Arguments:
            mach_regime (string): The Mach number regime in which the fuselage is operating
                Different dynamic models are needed for different Mach regimes.

        All other arguments are passed to the constructor of the dynamic model.
        """
        if mach_regime == 'subsonic':
            return FuselageAeroSubsonic(*args)
        else:
            raise NotImplementedError()


class FuselageAeroSubsonic(Model):
    def setup(self, fuse_geom, flight_state):
        """Model for aerodynamics of a fuselage at speeds up to Mach 1.

        Assumes a "missile-like" nose + cylinder fuselage.

        Based on Chapter 2.4 of Fleeman [1]

        References:
            [1] Fleeman, Eugene L, "Missile Design and System Engineering",
                American Institute of Aeronautics and Astronautics, 2012.
            [2] Jerger, J. J., "Systems Preliminary Design," Principles
                of Guided Missile Design, 1960.

        TODO "Unbounded" tags, ask Ned how to use them.
        """
        self.fuse_geom = fuse_geom
        self.flight_state = flight_state

        ### Declare Geometric Programming variables ###
        drag_fuse = Variable('drag_fuse', 'newton', 'Fuselage drag')

        # Drag coefficients
        C_D0_fuse = Variable('C_D0_fuse', '', 'Fuselage drag coefficient at zero lift')
        C_D0_fric = Variable('C_D0_fric', '', 'Fuselage skin friction drag coeff')
        C_D0_base = Variable('C_D0_base', '', 'Fuselage base drag coefficient')

        # Parameters "internal" to the model
        # Parameter in the skin friction equation provided in [1], which is originally
        # from Jerger [2].
        _jerger_param = Variable('_jerger_param', 0.053, 'lbf^0.2 ft^-0.2', 'Skin friction parameter')


        ### Declare Geometric Programming constraints. ###
        # Drag model constraints
        constraints_drag = [
            # Definition of C_D for compressible flow [1].
            drag_fuse == C_D0_fuse * fuse_geom['S_ref'] * flight_state['p_dyn'],

            # Fuselage drag is the sum of friction and base drag.
            # (no wave drag in the subsonic model).
            Tight([C_D0_fuse >= C_D0_fric + C_D0_base]),

            # Base drag model from Fleeman Ch. 2.4 [1]
            Tight([C_D0_base >= 0.12 + 0.13 * flight_state['mach']**2]),

            # Friction drag model that Fleeman cites from Jerger [2].
            C_D0_fric == (_jerger_param * fuse_geom['fineness']
                * (flight_state['mach'] / (flight_state['p_dyn'] * fuse_geom['length']))**0.2),
        ]

        constraints = constraints_drag
        constraints += [fuse_geom]
        constraints += [flight_state]
        return constraints
