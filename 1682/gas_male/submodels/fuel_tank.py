" fuel tank "
from gpkit import Model, Variable

class FuelTank(Model):
    """
    Returns the weight of the fuel tank.  Assumes a cylinder shape with some
    fineness ratio
    """
    def __init__(self, Wfueltot, **kwargs):

        W = Variable("W", "lbf", "fuel tank weight")
        f = Variable("f", 0.03, "-", "fraction fuel tank weight to fuel weight")
        mfac = Variable("m_{fac}", 1.1, "-", "fuel volume margin factor")
        rhofuel = Variable("\\rho_{fuel}", 6.01, "lbf/gallon",
                           "density of 100LL")
        Vol = Variable("\\mathcal{V}", "ft^3", "fuel tank volume")

        constraints = [W >= f*Wfueltot,
                       Vol/mfac >= Wfueltot/rhofuel,
                      ]

        Model.__init__(self, None, constraints, **kwargs)
