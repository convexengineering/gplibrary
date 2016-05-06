import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman')
plt.rcParams.update({'font.size':19})

def wind(M):
    M.substitutions.update({'V_{wind}': ('sweep', np.linspace(5, 25, 20))})
    sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
    V_wind = sol('V_{wind}')
    t_station = sol('t_{station}')
    

    plt.close()
    plt.plot(V_wind, t_station)
    plt.plot([20,20], [0,10], '--', color='r')
    plt.title('Endurance vs Wind Speed')
    plt.xlabel('Wind Speed on Station [m/s]')
    plt.ylabel('Time on Station [days]')
    plt.grid()
    plt.axis([5, 40, 0, 10])
    plt.savefig('tvsV_wind.png')