from gpkit import Model, parse_variables

class CFRPFabric(Model):
    """ Carbon Fiber Reinforced Plastic Fabric Material Properties

    Variables
    ---------
    rho         1.6         [g/m^3]         density of CFRP
    tmin        0.012       [in]            minimum gauge thickness
    tau         570         [MPa]           torsional stress limit

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{CFRP}}
    tmin        t_{\\mathrm{min-CFRP}}
    tau         \\tau_{\\mathrm{CFRP}}

    """
    def setup(self):
        exec parse_variables(CFRPFabric.__doc__)

class CFRPUD(Model):
    """ Carbon Fiber Reinforced Plastic Unidirectional Material Properties

    Variables
    ---------
    rho         1.6         [g/m^3]         density of CFRP
    E           2e7         [psi]           Youngs Modulus of CFRP
    sigma       1700        [MPa]           maximum stress limit of CFRP

    LaTex Strings
    -------------
    rho         \\rho_{\\mathrm{CFRP}}
    tmin        t_{\\mathrm{min-CFRP}}
    sigma       \\sigma_{\\mathrm{CFRP}}

    """
    def setup(self):
        exec parse_variables(CFRPUD.__doc__)

