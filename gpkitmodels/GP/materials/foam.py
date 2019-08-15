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
    def setup(self):
        exec(parse_variables(FoamHD.__doc__))

class FoamLD(Model):
    """ Foam low density material properties

    Constants
    ---------
    rho         0.024       [g/cm^3]        foam density

    LaTex Strings
    -------------
    rho             \\rho_{\\mathrm{foam}}

    """
    def setup(self):
        exec(parse_variables(FoamLD.__doc__))
