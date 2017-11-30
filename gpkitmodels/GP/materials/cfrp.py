from gpkit import Model, parse_variables

class CFRP(Model):
    """ Carbon Fiber Reinforced Plastic Material Properties

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
        exec parse_variables(CFRP.__doc__)

