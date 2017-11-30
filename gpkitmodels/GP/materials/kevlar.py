from gpkit import Model, parse_variables

class Kevlar(Model):
    """ Kevlar Material Properties

    Variables
    ---------
    rho         0.049       [g/m^3]         density of Kevlar
    tmin        0.012       [in]            minimum gauge thickness
    tau         200         [MPa]           torsional stress limit

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{Kevlar}}
    tmin        t_{\\mathrm{min-Kevlar}}
    tau         \\tau_{\\mathrm{Kevlar}}

    """
    def setup(self):
        exec parse_variables(Kevlar.__doc__)

