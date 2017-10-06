" constant taper chord "
import numpy as np

def c_bar(lam, N):
    "returns wing chord lengths for constant taper wing"
    eta = np.linspace(0, 1, N)
    c = (lam-1)*eta + 1
    cmac = 2.0/3.0*(1+lam+lam**2.0)/(1+lam)
    Sbar = 1 + lam
    return c, eta, cmac, Sbar
