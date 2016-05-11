import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import pdb
plt.rc('font', family='serif')
plt.rc('font', serif='Times New Roman')
plt.rcParams.update({'font.size':19})
plt.rcParams.update({'lines.linewidth':2})


def W_pay(M):
    M.substitutions.update({'MTOW': ('sweep', np.linspace(8+140, 18+140, 11))})
    #M.substitutions.update({'W_{pay}': ('sweep', np.linspace(8, 18, 11))})
    sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
    MTOW = sol('MTOW').magnitude
    W_pay = MTOW - np.multiply(np.ones([11]),140)
    #W_pay = sol('W_{pay}')
    t_station = sol('t_{station}')
    df = pd.read_csv('trimDrag.csv')
    df = df.dropna()
    CD10 = df['CD_interp'][df['W_pay_interp']==10]
    CD_norm = np.array(df['CD_interp'])/np.array(CD10)
    t_station = t_station/CD_norm


    plt.close()
    plt.plot(W_pay, t_station)
    plt.plot([10,10],[0,6], '--g')
    plt.plot([0,10],[6,6], '--g')
    plt.xlabel('Payload Weight [lbs]')
    plt.ylabel('Time on Station [days]')
    plt.grid()
    plt.axis([8, 18, 0, 10])
    plt.savefig('tvsW_pay.pdf')
