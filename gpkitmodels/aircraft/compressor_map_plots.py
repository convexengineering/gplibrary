from engine_components import FanMap, LPCMap
from gpkit import Model, units
import numpy as np
import matplotlib.pyplot as plt
from engine import EngineOnDesign
import warnings
from gpkit.small_scripts import mag

warnings.filterwarnings("error")

engineOnD = EngineOnDesign()
    
sol = engineOnD.localsolve(verbosity = 1, kktsolver="ldl")

mhtD, mltD, mFanBarD, mlcD, NlpcD, NhpcD, A5, A7 = engineOnD.sizing(sol)

#-----------------------------------------------------------------------------
#Fan Case


##mvec = np.linspace(100,1114,15)
##
##nvec = np.linspace(1,1.1,15)
##
##ntot = []
##pitot = []
##mtot = []
##
##for i in range(len(nvec)):
##    nans=[]
##    pians=[]
##    mans=[]
##    for j in range(len(mvec)):
##        fanmap = FanMap()
##        fanmap.substitutions.update({
##            'T_{ref}': 288.15,
##            'T_{t_2}': sol('T_{t_2}'),
##            'm_{fan}': mvec[j],
##            'm_{fan_D}': 10*sol('m_{core}'),
##            'P_{ref}': 101.325,
##            'P_{t_2}': sol('P_{t_2}'),
##            '\pi_{f_D}': sol('\pi_f'),
##            'N_{{bar}_Df}': 1,
##            'N_f': nvec[i],
##            'm_{fan_bar_D}': mFanBarD,
##            })
##        try:
##            solf = fanmap.localsolve(kktsolver="ldl", verbosity = 0)
##            nans.append(mag(solf('N_{barf}')))
##            mans.append(solf('m_{f}')/mFanBarD)
##            pians.append(mag(solf('\pi_f')))
##            #print solf('m_{f}')/mFanBarD
##        except (ValueError, RuntimeWarning):
##            print ["failed soln", nvec[i], mvec[j]]
##        
##    for i in range(len(mans)):
##        mans[i]= float(mans[i])
##    for i in range(len(pians)):
##        pians[i]= float(pians[i])
##    for i in range(len(nans)):
##        nans[i]= float(nans[i])
##    ntot.append(nans)
##    mtot.append(mans)
##    pitot.append(pians)
##
##for i in range(len(ntot)):
##       plt.plot(mtot[i],pitot[i], '-k')
##plt.title('Fan Compressor Map')
##plt.xlabel('Normalized Corrected Mass Flow')
##plt.ylabel('Pressure Ratio')
##plt.show()
##plt.show()

#------------------------------------------------------------
#LPC Case

mvec = np.linspace(10,252,15)

nvec = np.linspace(.9,1.4,15)

ntot = []
pitot = []
mtot = []

for i in range(len(nvec)):
    nans=[]
    pians=[]
    mans=[]
    for j in range(len(mvec)):
        lpcmap = LPCMap()
        lpcmap.substitutions.update({
            'T_{ref}': 288.15,
            'T_{t_2}': sol('T_{t_2}'),
            'm_{core}': mvec[j],
            'm_{core_D}': sol('m_{core}'),
            'P_{ref}': 101.325,
            'P_{t_2}': sol('P_{t_2}'),
            '\pi_{lc_D}': sol('\pi_{lc}'),
            'N_{{bar}_D_lc}': 1,
            'N_1': nvec[i],
            'm_{lc_D}': mFanBarD,
            })
        try:
            sollc = lpcmap.localsolve(kktsolver="ldl", verbosity = 0)
            nans.append(mag(sollc('N_{bar_lc}')))
            mans.append(sollc('m_{lc}')/mlcD)
            pians.append(mag(sollc('\pi_{lc}')))
            #print solf('m_{f}')/mFanBarD
        except (ValueError, RuntimeWarning):
            print ["failed soln", nvec[i], mvec[j]]
        
    for i in range(len(mans)):
        mans[i]= float(mans[i])
    for i in range(len(pians)):
        pians[i]= float(pians[i])
    for i in range(len(nans)):
        nans[i]= float(nans[i])
    ntot.append(nans)
    mtot.append(mans)
    pitot.append(pians)

for i in range(len(ntot)):
    plt.plot(mtot[i],pitot[i], '-k')
plt.title('Low Pressure Compressor Map')
plt.xlabel('Normalized Corrected Mass Flow')
plt.ylabel('Pressure Ratio')
plt.show()



#-------------------------------------------
#Run this once substitutions work!

##lpcmap = LPCMap()
##lpcmap.substitutions.update({
##    'T_{ref}': 1000,
##    'T_{t_2}': sol('T_{t_2}'),
##    'm_{core}': ('sweep', mvec),
##    'm_{core_D}': sol('m_{core}'),
##    'P_{ref}': 22,
##    'P_{t_2}': sol('P_{t_2}'),
##    '\pi_{lc_D}': sol('\pi_{lc}'),
##    'N_{{bar}_D_lc}': 1,
##    'N_1': ('sweep',nvec)
##    })
##sol = lpcmap.localsolve(kktsolver="ldl", verbosity = 4, skipsweepfailures=True)
##print "Done"
##1/0
##
##mvec = np.linspace(20,100,2)
##
##nvec = np.linspace(1,20,1)
##
##fanmap = FanMap()
##fanmap.substitutions.update({
##    'T_{ref}': 1000,
##    'T_{t_2}': sol('T_{t_2}'),
##    'm_{fan}': ('sweep', mvec),
##    'm_{fan_D}': sol('alpha')*sol('m_{core}'),
##    'P_{ref}': 22,
##    'P_{t_2}': sol('P_{t_2}'),
##    '\pi_{f_D}': sol('\pi_f'),
##    'N_{{bar}_Df}': 1,
##    'N_f': ('sweep',nvec)
##    })
##sol = fanmpa.localsolve(kktsolver="ldl", verbosity = 4, skipsweepfailures=True)

