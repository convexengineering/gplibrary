"Functions that return atmospheric quantities given an altitude or range of altitudes"
import numpy as np

def get_atmosphere(h):
    """Returns atmospheric properties for a given altitude up to 20km"""

    g = 9.81        # [m/s^2]
    L = 0.0065      # [K/m]
    R = 287         # [J/(kg*K)]
    T_0 = 288.15    # [K]
    p_0 = 101325    # [Pa]
    C_1  = 1.458E-6 # [kg/(m*s*K**0.5)]
    S = 110.4       # [K]

    if h <= 11000:
        th = g/(R*L)          # [-]
        T = T_0 - L*h          # [K]
        p = p_0*(T/T_0)**th  # [Pa]
    elif h <= 20000:
        T = 216.7   # [K]
        p_1 = 22630 # [Pa]
        p = p_1*np.exp(-g/(R*T)*(h-11000)) # [Pa]
    else:
        raise Exception('Model only valid up to 20km')

    rho = p/(R*T)         # [kg/m^3]
    mu = C_1*T**1.5/(T+S) # [kg/(m*s)]

    return T, p, rho, mu

def get_atmosphere_vec(hvec):
    """Returns atmospheric properties for a vector of altitudes"""

    Tvec, pvec, rhovec, muvec = [], [], [], []
    for h in hvec:
        T, p, rho, mu = get_atmosphere(h)
        Tvec += [T]
        pvec += [p]
        rhovec += [rho]
        muvec += [mu]

    A = {
         'T': Tvec,
         'p': pvec,
         'rho': rhovec,
         'mu': muvec
        }
    return A

def get_T(h):
    T, _, _, _ = get_atmosphere(h)
    return T

def get_rho(h):
    _, _, rho, _ = get_atmosphere(h)
    return rho
        
def get_mu(h):
    _, _, _, mu = get_atmosphere(h)
    return mu
