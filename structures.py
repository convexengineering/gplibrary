from numpy import pi
from gpkit import Model, Variable


class ColumnBuckling(Model):
    """Euler column buckling, from https://en.wikipedia.org/wiki/Buckling

    For both ends pinned (hinged, free to rotate), K = 1.0.
    For both ends fixed, K = 0.50.
    For one end fixed and the other end pinned, K ~= 0.699.
    For one end fixed and the other end free to move laterally, K = 2.0.
    """
    def setup(self):
        F = Variable("F", "lbf", "vertical load on column")
        E = Variable("E", "GPa", "modulus of elasticity")
        I = Variable("I", 1e-9, "m^4", "area moment of inertia")
        L = Variable("L", 1, "m", "unsupported length of column")
        K = Variable("K", 2.0, '-', "column effective length factor")
        return 1/F, [F <= pi**2*E*I / (K*L)**2]


if __name__ == "__main__":
    from materials import Aluminum6061
    M = ColumnBuckling().merge(Aluminum6061())
    M.solve()
