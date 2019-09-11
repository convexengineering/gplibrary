from gpkit import Model, parse_variables

class FoamHD(Model):
    """ Foam high density material properties

    Constants
    ---------
    rho         0.036       [g/cm^3]        foam density

    LaTex Strings
    -------------
    rho             \\rho_{\\mathrm{foam}}

    """
    @parse_variables(__doc__, globals())
    def setup(self):
        pass

class FoamLD(Model):
    """ Foam low density material properties

    Constants
    ---------
    rho         0.024       [g/cm^3]        foam density

    LaTex Strings
    -------------
    rho             \\rho_{\\mathrm{foam}}

    """
    @parse_variables(__doc__, globals())
    def setup(self):
        pass
