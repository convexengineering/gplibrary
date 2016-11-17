" air properties function of altitude "

def get_density(altitude):
    """
    inputs
    ------
    altitude: atmosphereic altitude ft

    output
    ------
    rho: air density
    """

    if not hasattr(altitude, "__len__"):
        altitude = [altitude]

    rho = []
    for a in altitude:
        h = a*0.3048
        Tsl = 288.15 # sea level temperature [K]
        Tatm = Tsl - 0.0065*h # temperature at atmosphere [K]

        Rspec = 287.058 # specific gas constant of air [J/kg/K]
        psl = 101325 # pressure at sea level [Pa]

        rho.append(psl*Tatm*(5.257-1)/Rspec/(Tsl**5.257))

    if len(rho) == 1:
        rho = rho[0]

    return rho
