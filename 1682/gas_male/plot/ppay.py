import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman')
plt.rcParams.update({'font.size':19})

def P_pay(M): 
    #M.substitutions.update({'P_{pay}': ('sweep', np.linspace(5, 200, 20))})
    #M.substitutions.update({'MTOW' : np.array(148)})
    sol = M.solve(solver='mosek', verbosity=3, skipsweepfailures=True)
    P_pay = sol('P_{pay}')
    t_station = sol('t_{station}')


    plt.close()
    plt.plot(P_pay, t_station)
    plt.title('Endurance vs Payload Power')
    plt.xlabel('Payload Power [Watts]')
    plt.ylabel('Time on Station [days]')
    plt.grid()
    plt.axis([5, 200, 0, 10])
    plt.savefig('tvsP_pay.png')