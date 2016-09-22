"""
A simple beam example with fixed geometry. Solves the discretized
Euler-Bernoulli beam equations for a constant distributed load
"""
import numpy as np
from gpkit import Variable, Model, VectorVariable

def return_loading(lam, N):
    eta = np.linspace(0, 1, N)
    c = 2/(1+lam)*(1+(lam-1)*eta)
    return c

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

        cb = return_loading(0.5, N)
        cbavg = (cb[:-1] + cb[1:])/2
        W_cent = Variable("W_{cent}", 124, "lbf", "center weight")
        N_fac = Variable("N_{fac}", 5, "count", "safety load factor")
        b = Variable("b", 23, "ft", "wing span")
        I = VectorVariable(N-1, "I", "m^4", "spar x moment of inertia")
        E = Variable("E", 181, "GPa", "carbon fiber young's modulus")
        dx = Variable("dx", "m", "Length of an element")
        cbar = VectorVariable(N, "\\bar{c}", cb, "-",
                              "normalized distributed load at each point")
        h = VectorVariable(N-1, "h", "ft", "spar height")
        t = VectorVariable(N-1, "t", "ft", "thickness")
        tau = Variable("\\tau", 0.115, "-", "wing t/c ratio")
        S = Variable("S", 21.38, "ft**2", "wing area")
        w = VectorVariable(N-1, "w", "in", "spar width")
        Q = VectorVariable(N, "Q", "N/m", "")
        V = VectorVariable(N, "V", "N", "Internal shear")
        V_tip = Variable("V_{tip}", 1e-10, "N", "Tip loading")
        M = VectorVariable(N, "M", "N*m", "Internal moment")
        M_tip = Variable("M_{tip}", 1e-10, "N*m", "Tip moment")
        th = VectorVariable(N, "\\theta", "-", "Slope")
        th_base = Variable("\\theta_{base}", 1e-10, "-", "Base angle")
        delta = VectorVariable(N, "\\delta", "m", "Displacement")
        delta_base = Variable("\\delta_{base}", 1e-10, "m", "Base deflection")
        delta_tipmax = Variable("\\delta_{tip-max}", 0.2, "-",
                                "Max tip deflection ratio")
        m = VectorVariable(N-1, "m", "kg", "spar mass")
        m_tot = Variable("m_{tot}", "kg", "spar mass")
        rho = Variable("\\rho", 1.76, "g/cm**3", "carbon fibre density")
        sigma = Variable("\\sigma", 475e6, "Pa", "stress limit of CF")
        # below: trapezoidal integration to form a piecewise-linear
        #        approximation of loading, shear, and so on
        # shear and moment increase from tip to base (left > right)
        constraints = [Q >= cbar*W_cent*N_fac/b,
                       I <= w*h**3/12,
                       m >= rho*w*h*b/2/(N-1),
                       m_tot >= m.sum(),
                       h <= t,
                       t <= S/b*cbavg*tau,
                       V[:-1] >= V[1:] + 0.5*dx*(Q[:-1] + Q[1:]),
                       V[-1] >= V_tip,
                       M[:-1] >= M[1:] + 0.5*dx*(V[:-1] + V[1:]),
                       M[-1] >= M_tip,
                       M[:-1]/w/h**2 <= sigma,
                       th[0] >= th_base,
                       th[1:] >= th[:-1] + 0.5*dx*(M[1:] + M[:-1])/E/I,
                       delta[0] >= delta_base,
                       delta[1:] >= delta[:-1] + 0.5*dx*(th[1:] + th[:-1]),
                       b/2 == (N-1)*dx,
                       delta[-1]/(b/2) <= delta_tipmax]
        # slope and displacement increase from base to tip (right > left)
        # minimize tip displacement (the last w)
        Model.__init__(self, m_tot, constraints, **kwargs)

if __name__ == "__main__":
    beam = Beam(5)
    sol = beam.solve("mosek")
