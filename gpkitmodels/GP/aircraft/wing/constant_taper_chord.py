" constant taper chord "
import numpy as np

def c_bar(lam, N):
    "returns wing chord lengths for constant taper wing"
    eta = np.linspace(0, 1, N)
    c = 2/(1+lam)*(1+(lam-1)*eta)
    cbarmac = 2./3*(1+lam+lam**2)/(1+lam)
    deta = np.diff(eta)
    return c, eta, deta, cbarmac

