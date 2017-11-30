from gpkit import Model, parse_variables

class Foam(Model):
    """ Foam material properties

    Variables
    ---------
    rho         0.036       [g/cm^3]        foam density

    LaTex Strings
    -------------
    rho             \\rho_{\\mathrm{foam}}

    """
    def setup(self):
        exec parse_variables(Foam.__doc__)
