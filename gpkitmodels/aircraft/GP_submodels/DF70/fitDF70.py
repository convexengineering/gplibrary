import pandas as pd
import numpy as np
from numpy import logspace, log, log10
import matplotlib.pyplot as plt
from gpfit.fit import fit

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
yfit = cstrt.evaluate(x)

fig, ax = plt.subplots()
ax.plot(RPM, BSFC, "o", markerfacecolor="None")
ax.plot(RPM, np.exp(yfit)*BSFCmin)
ax.set_xlabel("$RPM$")
ax.set_ylabel("$BSFC$ [lb/hp/hr]")
ax.legend(["Manufacture Data", "GP approximation"])
ax.grid()
fig.savefig("rpmtobsfcfit.pdf")

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
cstrt, rms_error = fit(x,y,K,Type)
yfit = cstrt.evaluate(x)

fig, ax = plt.subplots()
ax.plot(RPM, P, "o", markerfacecolor="None")
ax.plot(RPM, np.exp(yfit)*P_max)
ax.set_xlabel("$RPM$")
ax.set_ylabel("Shaft Power [kW]")
ax.legend(["Manufacture Data", "GP approximation"])
ax.grid()
fig.savefig("rpmtopowerfit.pdf")

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

#h = np.delete(h,[0])
#P = np.delete(P,[0])

# Fitting Power vs. RPM
# P = df['kW']
# P_MSL = np.amax(df['kW'])
# Pnorm = P/P_MSL
# L = 1-Pnorm

# the variables you input into fit() can't have negatives here or Inf or -Inf
# logL = log(L)
# logh = log(h)



# plt.xlabel('log(h)')
# plt.ylabel('log(1-Pnorm)')
# #plt.show()

#cstrt, rms_error = fit(logh[1:-1],logL[1:-1],K,Type)


