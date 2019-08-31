" box spar "
from gpkit import parse_variables, SignomialsEnabled
from gpkitmodels.GP.aircraft.wing.boxspar import BoxSpar as BoxSparGP

#pylint: disable=exec-used, undefined-variable, unused-argument, invalid-name

class BoxSpar(BoxSparGP):
    """ Box Spar Model

    Variables of length N-1
    -----------------------
    J                       [m^4]       spar x polar moment of inertia

    """
    @parse_variables(__doc__, globals())
    def setup(self, N, surface):
        self.boxspar = BoxSparGP.setup(self, N=N, surface=surface)

        cave = self.cave
        tau = self.tau
        w = self.w
        tshear = self.tshear

        with SignomialsEnabled():
            constraints = [J <= cave*tau*w*tshear/3*(cave*tau + w)]

        return self.boxspar, constraints
