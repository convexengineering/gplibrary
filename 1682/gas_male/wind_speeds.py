import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
plt.rcParams.update({'font.size':19})

def get_windspeed(latitude, perc, altitude):
    filename = None
    h = altitude*0.3048
    p = 101325.0/100*(1 - 2.25577e-5*h)**5.25588
    pressures = [5, 10, 30] + range(50, 1050, 50)
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
            wind.append(get_windspeed(l, p, 16000))
        ax.plot(wind, np.arange(0, 70, 1))

    ax.set_ylabel("Latitude [deg]")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylim([0, 70])
    ax.grid()
    ax.legend(["%d Percentile Winds" % a for a in [80, 90, 95]], loc=2, fontsize=15)
    fig.savefig("latvswindh16.pdf", bbox_inches="tight")

    fig, ax = plt.subplots()
    for p in [80, 90, 95]:
        wind = []
        for l in np.arange(0, 70, 1):
            wind.append(get_windspeed(l, p, 50000))
        ax.plot(wind, np.arange(0, 70, 1))

    ax.set_ylabel("Latitude [deg]")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylim([0, 70])
    ax.grid()
    ax.legend(["%d Percentile Winds" % a for a in [80, 90, 95]], loc=2, fontsize=15)
    fig.savefig("latvswindh50.pdf", bbox_inches="tight")

    fig, ax = plt.subplots()
    alt = [350, 1800, 3200, 4800, 6400, 8000, 9900, 11700, 13800, 15961, 18288, 20800, 23500, 26600, 30000, 34000, 38600, 44300, 51800, 63400, 71000, 85000]
    for p in [80, 90, 95]:
        wind = []
        for h in alt:
            wind.append(get_windspeed(45, p, h))
        ax.plot(wind, alt)
    ax.set_ylabel("Altitude [ft]")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylim([0, 100000])
    ax.grid()
    ax.legend(["%d Percentile Winds" % a for a in [80, 90, 95]], loc=1, fontsize=15)
    fig.savefig("altvswindl45.pdf", bbox_inches="tight")

    fig, ax = plt.subplots()
    for p in [80, 90, 95]:
        wind = []
        for h in alt:
            wind.append(get_windspeed(30, p, h))
        ax.plot(wind, alt)
    ax.set_ylabel("Altitude [ft]")
    ax.set_xlabel("Wind speed [m/s]")
    ax.set_ylim([0, 100000])
    ax.grid()
    ax.legend(["%d Percentile Winds" % a for a in [80, 90, 95]], loc=1, fontsize=15)
    fig.savefig("altvswindl30.pdf", bbox_inches="tight")
