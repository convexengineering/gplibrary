import numpy as np
from numpy import sin, arctan, tan, cos, arcsin, arccos, deg2rad, rad2deg
import matplotlib.pyplot as plt

def get_Eirr(latitude, day):
    """
    day is juilian day, measured from Jan 1st
    latitude is in degrees
    Returns:
    -------
    ESirr: Solar energy per unit area of from the sun Whr/m^2
    tday: time of daylight
    tnight: time of daylight
    """
    assert isinstance(day, int)

    beta = 2*np.pi*(day-1)/365
    lat = deg2rad(latitude)
    delta = (0.006918 - 0.399912*cos(beta) + 0.070257*sin(beta) -
             0.006758*cos(2*beta) + 0.000907*sin(2*beta) -
             0.002697*cos(3*beta) + 0.00148*sin(3*beta))
    tstart = 12/np.pi*arccos(-tan(delta)*tan(lat))
    tend = -tstart
    t = np.linspace(tstart, tend, 50)
    costhsun = sin(delta)*sin(lat) + cos(delta)*cos(lat)*cos(2*np.pi*t/24)

    r0 = 149.597e6 # avg distance from earth to sun km
    Reo = r0*(1 + 0.017*sin(2*np.pi*(day-93)/365))
    Esun = 63372630 # energy from sun surface W/m^2
    Rsun = 695842 # radius of the sun, km
    E0 = Esun*4*np.pi*Rsun**2/4/np.pi/Reo**2
    tau = np.exp(-0.175/costhsun)
    E = E0*costhsun# *tau
    fig, ax = plt.subplots()
    ax.plot(t, E)
    ax.set_xlabel("Time [hr]")
    ax.set_ylabel("Available Solar Power [W/m$^2$]")
    ax.grid()
    fig.savefig("lat45.pdf")
    return np.trapz(E)*(abs(tend-tstart))/50.0, tstart*2, 24-tstart*2
