" air properties function of altitude "
from numpy import sqrt

def get_airvars(altitude):
    """
    inputs
    ------
    altitude: atmosphereic altitude ft

    output
    ------
    rho: air density
    mu: dynamic viscosity
    """

    if not hasattr(altitude, "__len__"):
        altitude = [altitude]

    rho = []
    mu = []
    for a in altitude:
        h = a*0.3048
        Tsl = 288.15 # sea level temperature [K]
        Tatm = Tsl - 0.0065*h # temperature at atmosphere [K]

        Rspec = 287.058 # specific gas constant of air [J/kg/K]
        psl = 101325 # pressure at sea level [Pa]
        musl = 1.458e-6

        rho.append(psl*Tatm**(5.257-1)/Rspec/(Tsl**5.257))
        mu.append(musl*sqrt(Tatm)/(1+110.4/Tatm))

    if len(rho) == 1:
        rho = rho[0]
        mu = mu[0]

    return rho, mu
