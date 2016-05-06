import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman')
plt.rcParams.update({'font.size':19})

def fixed(t, V, NSeg, NClimb1, NCruise1, NClimb2, NLoiter, NCruise2,
          mStart, mEndClimb, mEndCruise, mEndClimb2, mEndLoiter):
    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, V)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 40 
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('Flight Velocity [m/s]')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('Vvsmissionnowind.png')