from gpkit import Model, parse_variables

class Ti6Al4V(Model):
    """ Solid Ti64 Material Properties

    Constants
    ---------
    rho         4.512       [g/cm^3]        density of CFRP
    tmin        0.3048      [mm]            minimum gauge thickness
    tau         900         [MPa]           torsional stress limit
    E           104         [GPa]           Youngs modulus
    sigma       900         [MPa]           max stress
    G           40          [GPa]           shear modulus

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{CFRP}}
    tmin        t_{\\mathrm{min-CFRP}}
    tau         \\tau_{\\mathrm{CFRP}}

    """
    def setup(self):
        exec parse_variables(Ti6Al4V.__doc__)