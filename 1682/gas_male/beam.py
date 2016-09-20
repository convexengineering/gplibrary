"""
A simple beam example with fixed geometry. Solves the discretized
Euler-Bernoulli beam equations for a constant distributed load
"""
import numpy as np
from gpkit import Variable, Model, VectorVariable

def return_loading(N):
    y = np.linspace(0, 1, N)
    q = np.sqrt(1 - y**2)
    return q

class Beam(Model):
    """Discretization of the Euler beam equations for a distributed load.

    Arguments
    ---------
    N : int
        Number of finite elements that compose the beam.
    L : float
        [m] Length of beam.
    EI : float
        [N m^2] Elastic modulus times cross-section's area moment of inertia.
    q : float or N-vector of floats
        [N/m] Loading density: can be specified as constants or as an array.
    """
    def __init__(self, N=4, **kwargs):

        wingload = return_loading(N)
        MTOW = Variable("MTOW", 150, "lbf", "Max take off weight")
        b = Variable("b", 25, "ft", "wing span")
        EI = Variable("EI", "N*m^2")
        dx = Variable("dx", "m", "Length of an element")
        q = VectorVariable(N, "q", wingload*4*MTOW.value/np.pi/b.value, "N/m",
                           "Distributed load at each point")
        V = VectorVariable(N, "V", "N", "Internal shear")
        V_tip = Variable("V_{tip}", 0.00001, "N", "Tip loading")
        M = VectorVariable(N, "M", "N*m", "Internal moment")
        M_tip = Variable("M_{tip}", 0.00001, "N*m", "Tip moment")
        th = VectorVariable(N, "\\theta", "-", "Slope")
        th_base = Variable("\\theta_{base}", 0.000001, "-", "Base angle")
        w = VectorVariable(N, "w", "m", "Displacement")
        w_base = Variable("w_{base}", 0.000001, "m", "Base deflection")
        delta_tipmax = Variable("\\delta_{tip-max}", 0.2, "-",
                                "Max tip deflection ratio")
        # below: trapezoidal integration to form a piecewise-linear
        #        approximation of loading, shear, and so on
        # shear and moment increase from tip to base (left > right)
        constraints = [V[:-1] >= V[1:] + 0.5*dx*(q[:-1] + q[1:]),
                       V[-1] >= V_tip,
                       M[:-1] >= M[1:] + 0.5*dx*(V[:-1] + V[1:]),
                       M[-1] >= M_tip,
                       th[0] >= th_base,
                       th[1:] >= th[:-1] + 0.5*dx*(M[1:] + M[:-1])/EI,
                       w[0] >= w_base,
                       w[1:] >= w[:-1] + 0.5*dx*(th[1:] + th[:-1]),
                       b/2 == (N-1)*dx,
                       w[-1]/(b/2) <= delta_tipmax]
        # slope and displacement increase from base to tip (right > left)
        # minimize tip displacement (the last w)
        Model.__init__(self, EI, constraints, **kwargs)


    def test(self):
        sol = self.solve()
        L, EI, q = sol("L"), sol("EI"), sol("q")
        N = len(q)
        q = q[0]  # assume uniform loading for the check below
        x = np.linspace(0, L, N)  # position along beam
        w_gp = sol("w")  # deflection along beam
        w_exact = q/(24.*EI) * x**2 * (x**2 - 4*L*x + 6*L**2)  # analytic soln

if __name__ == "__main__":
    Beam().test()
