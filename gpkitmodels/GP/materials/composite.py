from gpkit import Model, parse_variables

class CFRPFabric(Model):
    """ Carbon Fiber Reinforced Plastic Fabric Material Properties

    Constants
    ---------
    rho         1.6         [g/cm^3]        density of CFRP
    tmin        0.3048      [mm]            minimum gauge thickness
    tau         570         [MPa]           torsional stress limit
    E           150         [GPa]           Youngs modulus
    sigma       400         [MPa]           max stress
    G           2           [GPa]           shear modulus

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{CFRP}}
    tmin        t_{\\mathrm{min-CFRP}}
    tau         \\tau_{\\mathrm{CFRP}}

    """
    @parse_variables(__doc__, globals())
    def setup(self):
        pass

class CFRPUD(Model):
    """ Carbon Fiber Reinforced Plastic Unidirectional Material Properties

    Constants
    ---------
    rho         1.6         [g/cm^3]        density of CFRP
    E           137         [GPa]           Youngs Modulus of CFRP
    sigma       1700        [MPa]           maximum stress limit of CFRP
    tmin        0.1         [mm]            minimum gague thickness

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{CFRP}}
    tmin        t_{\\mathrm{min-CFRP}}
    sigma       \\sigma_{\\mathrm{CFRP}}

    """
    @parse_variables(__doc__, globals())
    def setup(self):
        pass

class Kevlar(Model):
    """ Kevlar Material Properties

    Constants
    ---------
    rho         0.049       [g/cm^3]        density of Kevlar
    tmin        0.012       [in]            minimum gauge thickness
    tau         200         [MPa]           torsional stress limit

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{Kevlar}}
    tmin        t_{\\mathrm{min-Kevlar}}
    tau         \\tau_{\\mathrm{Kevlar}}

    """
    @parse_variables(__doc__, globals())
    def setup(self):
        pass
