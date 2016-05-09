"Minimize x + 1/x while keeping x greater than x_min."
from gpkit import Variable, Model


class SimpleGP(Model):

    def __init__(self, **kwargs):
        xmin = Variable('x_{min}', 1)
        x = Variable('x')

        constraints = [x >= xmin]
        Model.__init__(self, x + 1/x, constraints, **kwargs)
