import pandas as pd
import numpy as np
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
from gpfit.fit import fit
plt.rcParams.update({'font.size':19})

np.random.seed(0)

# Fitting BSFC vs. RPM
df = pd.read_csv('Dataset_BSFC_kgKwh.csv')
RPM = df['RPM']
RPMmax = np.amax(RPM)
BSFC = df['BSFC']
BSFCmin = np.amin(BSFC)
x = np.array(log(RPM/RPMmax))
y = np.array(log(BSFC/BSFCmin))
Type = 'SMA'
K = 2

cstrt, rmserror = fit(x,y,K,Type)
print "RMS error = %.4f" % rmserror
yfit = cstrt.evaluate(x)

fig, ax = plt.subplots()
ax.plot(RPM, BSFC, "o", markerfacecolor="None")
ax.plot(RPM, np.exp(yfit)*BSFCmin)
ax.set_xlabel("$RPM$")
ax.set_ylabel("$BSFC$ [lb/hp/hr]")
ax.legend(["Manufacture Data", "GP approximation"])
ax.grid()
fig.savefig("rpmtobsfcfit.pdf", bbox_inches="tight")

# Fitting Power vs. RPM
df = pd.read_csv('Dataset_Power_Kw.csv')
RPM = df['RPM']
RPM_max = np.amax(RPM)
P = df['P']
P_max = np.amax(P)
y = np.array(log(P/P_max))
x = np.array(log(RPM/RPM_max))
Type = 'SMA'
K = 1
cstrt, rmserror = fit(x,y,K,Type)
print "RMS error = %.4f" % rmserror
yfit = cstrt.evaluate(x)

fig, ax = plt.subplots()
ax.plot(RPM, P, "o", markerfacecolor="None")
ax.plot(RPM, np.exp(yfit)*P_max)
ax.set_xlabel("$RPM$")
ax.set_ylabel("Shaft Power [kW]")
ax.legend(["Manufacture Data", "GP approximation"])
ax.grid()
fig.savefig("rpmtopowerfit.pdf", bbox_inches="tight")

# Fitting Torque vs. RPM
# df = pd.read_csv('Dataset_Torque_Nm.csv')
# RPM = df['RPM']
# RPM_max = np.amax(RPM)
# Q = df['Q']
# Q_max = np.amax(Q)
# logRPM = log(RPM/RPM_max)
# logQ = log(Q/Q_max)
# Type = 'SMA'
# K = 1
# cstrt, rms_error = fit(logRPM,logQ,K,Type)

# fit lapse rate
df = pd.read_csv("DF35_maxPvh.csv")
u = df["ft"]
w = df["kW"]/max(df["kW"])
x = np.array(u)
y = np.array(w)

A = np.vstack([x, np.ones(len(x))]).T
m, c = np.linalg.lstsq(A, y)[0]
print "Equation: y = %.4gx + %.4f" % (m, c)
fig, ax = plt.subplots()
ax.plot(x, y, 'o', label='RCV Engine Data', markerfacecolor="None")
ax.plot(x, m*x + c, label='Fitted Line')
ax.set_ylabel("Engine Lapse Rate")
ax.set_xlabel("Altitude [ft]")
ax.legend()
ax.grid()
fig.savefig("lapseline.pdf", bbox_inches="tight")