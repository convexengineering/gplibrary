import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

def get_windspeed(latitude, perc, altitude):
    filename = None
    h = altitude*0.3048
    p = 101325/100*(1 - 2.25577e-5*h)**5.25588
    pressures = [550, 650, 100, 30]
    for mb in pressures:
        if abs(mb-p) < 15:
            filename = "wind%d.csv" % mb
    if filename:
        df = pd.read_csv(filename)
        wind = df[df["Latitude"] == latitude]["perc%d" % perc].item()
    else:
        print "Altitude not close enough to known values"
        wind = 20
    return wind

if __name__ == "__main__":
    fig, ax = plt.subplots()
    for p in [80, 90, 95]:
        wind = []
        for l in np.arange(0, 70, 1):
            wind.append(get_windspeed(l, p, 50000))
        ax.plot(wind, np.arange(0, 70, 1))

    ax.set_xlabel("Latitude [deg]")
    ax.set_ylabel("Wind speed [m/s]")
    ax.set_ylim([0, 70])
    ax.grid()
    ax.legend(["%d Percentile Winds" % a for a in [80, 90, 95]])
    fig.savefig("latvswindh50.pdf")
