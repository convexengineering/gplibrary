import matplotlib.pyplot as plt
import numpy as np
from gasmale import GasPoweredMALE

plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman')
plt.rcParams.update({'font.size':19})

def altitude(M): 
    # plot for h vs h_station
    h_station = np.linspace(14000, 22000, 20)
    t_station = []
    rho = []
    mu = []
    a = []
    MTOW = []
    P_shafttot = []
    m_dotfuel = []
    P_shaftmax = []
    RPM = []
    eta_prop = []
    BSFC = []
    CD = []
    V = []
    W_fueltot = []

    for h in h_station:
        M = GasPoweredMALE(h)
        M.substitutions.update({'S': S_fix})
        #M.substitutions.update({'S_{fuse}': Sfuse_fix})
        M.substitutions.update({'b': b_fix})
        #M.substitutions.update({'l_{fuse}': lfuse_fix})
        M.substitutions.update({'Vol_{fuse}' : Volfuse_fix+0.0001})
        M.substitutions.update({'Vol_{fuel}' : Volfuel_fix+0.0001})
        #M.substitutions.update({'h_{spar}': hspar_fix})
        #del M.substitutions["t_{station}"]
        M.cost = 1/M["t_{station}"]
        sol = M.solve('mosek', verbosity=0)
        t_station.append(sol('t_{station}').magnitude)
        rho.append(sol('\\rho')[3].magnitude)
        mu.append(sol('\\mu')[3].magnitude)
        MTOW.append(sol('MTOW').magnitude)
        P_shafttot.append(sol('P_{shaft-tot}').magnitude)
        a.append(np.sqrt(sol('T_{atm}')[3].magnitude*1.4*286))
        m_dotfuel.append(sol('m_{dot-fuel}').magnitude)
        P_shaftmax.append(sol('P_{shaft-max}').magnitude)
        RPM.append(sol('RPM').magnitude)
        eta_prop.append(sol('\\eta_{prop}').magnitude)
        BSFC.append(sol('BSFC').magnitude)
        CD.append(sol('C_D').magnitude)
        V.append(sol('V').magnitude)
        W_fueltot.append(sol('W_{fuel-tot}').magnitude)


    plt.close()
    plt.plot(h_station, t_station)
    plt.title('Endurance vs Altitude')
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('Time on Station [days]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 10])
    plt.savefig('tvsh_station.png')
    
    plt.close()
    plt.plot(h_station, W_fueltot)
    plt.title('Fuel Weight vs Altitude')
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('Fuel Weight [lbs]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 100])
    plt.savefig('W_fuelvsh_station.png')
    
    plt.close()
    plt.plot(h_station, MTOW)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('MTOW [lbs]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 200])
    plt.savefig('MTOWvsh_station.png')
    
    plt.close()
    plt.plot(h_station, P_shafttot)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('P_shafttot [hp]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 5])
    plt.savefig('P_shafttotvsh_station.png')
    
    plt.close()
    plt.plot(h_station, m_dotfuel)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('m_dotfuel [lb/sec]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 0.0002])
    plt.savefig('m_dotfuelvsh_station.png')
    
    plt.close()
    plt.plot(h_station, P_shaftmax)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('P_shaftmax [hp]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 5])
    plt.savefig('P_shaftmaxvsh_station.png')
    
    plt.close()
    plt.plot(h_station, BSFC)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('BSFC [kg/kW/hr]')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 1.5])
    plt.savefig('BSFCvsh_station.png')
    
    plt.close()
    plt.plot(h_station, eta_prop)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('eta_prop')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0.7, 0.8])
    plt.savefig('eta_propvsh_station.png')
    
    plt.close()
    plt.plot(h_station, RPM)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('RPM')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 9000])
    plt.savefig('RPMvsh_station.png')
    
    plt.close()
    plt.plot(h_station, CD)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('CD')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 0.2])
    plt.savefig('CDvsh_station.png')
    
    plt.close()
    plt.plot(h_station, V)
    plt.xlabel('Altitude at Station [ft]')
    plt.ylabel('V')
    plt.grid()
    plt.axis([min(h_station), max(h_station), 0, 30])
    plt.savefig('Vvsh_station.png')