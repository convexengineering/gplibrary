"""Fuselage aerodynamics model for small, fast aircraft."""

from math import pi
from gpkit import Model, Variable, parse_variables
from gpkit.constraints.tight import Tight

class FuselageAeroSubsonic(Model):
    def setup(self):
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
        ### Declare Geometric Programming variables ###
        drag_fuse = Variable('drag_fuse', 'newton', 'Fuselage drag')

        # Fuselage dimensions
        len_fuse = Variable('len_fuse', 'meter', 'Fuselage length')
        diam_fuse = Variable('diam_fuse', 'meter', 'Fuselage diameter')
        fineness_fuse = Variable('fineness_fuse', '', 'Fuselage fineness ratio (length/diameter)')
        S_ref = Variable('S_ref', 'm^2', 'Fuselage cross section reference area')

        # Flight conditions
        mach = Variable('mach', '', 'Freestream Mach number')
        p_dyn = Variable('p_dyn', 'pascal', 'Freestream dynamic pressure')
        p_static = Variable('p_static', 'pascal', 'Freestream static pressure')
        _gamma_air = Variable('_gamma_air', 1.40, '', 'Ratio of specific heats of air')

        # Drag coefficients
        C_D0_fuse = Variable('C_D0_fuse', '', 'Fuselage drag coefficient at zero lift')
        C_D0_fric = Variable('C_D0_fric', '', 'Fuselage skin friction drag coeff')
        C_D0_base = Variable('C_D0_base', '', 'Fuselage base drag coefficient')

        # Parameters "internal" to the model
        # Parameter in the skin friction equation provided in [1], which is originally
        # from Jerger [2].
        _jerger_param = Variable('_jerger_param', 0.053, 'lbf^0.2 ft^-0.2', 'Skin friction parameter')


        ### Declare Geometric Programming constraints. ###
        # Geometric constraints
        constraints_geom = [
            fineness_fuse == len_fuse / diam_fuse,
            S_ref == (pi / 4) * diam_fuse**2,
        ]

        # Drag model constraints
        constraints_drag = [
            # Dynamic pressure
            # Source: https://en.wikipedia.org/wiki/Dynamic_pressure#Compressible_flow
            p_dyn == 0.5 * _gamma_air * p_static * mach**2,

            # Definition of C_D for compressible flow [1].
            drag_fuse == C_D0_fuse * S_ref * p_dyn,

            # Fuselage drag is the sum of friction and base drag.
            # (no wave drag in the subsonic model).
            Tight([C_D0_fuse >= C_D0_fric + C_D0_base]),

            # Base drag model from Fleeman Ch. 2.4 [1]
            Tight([C_D0_base >= 0.12 + 0.13 * mach**2]),

            # Friction drag model that Fleeman cites from Jerger [2].
            C_D0_fric == _jerger_param * fineness_fuse * (mach / (p_dyn * len_fuse))**0.2,
        ]

        constraints = constraints_geom + constraints_drag
        return constraints
