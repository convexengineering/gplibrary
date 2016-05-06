import numpy as np
import matplotlib.pyplot as plt
plt.rc('font', family='serif') 
plt.rc('font', serif='Times New Roman')
plt.rcParams.update({'font.size':19})

def mission(sol, t, V, NSeg, NClimb1, NCruise1, NClimb2, NLoiter, NCruise2,
            mStart, mEndClimb, mEndCruise, mEndClimb2, mEndLoiter):
    eta_prop = sol('\\eta_{prop}')
    BSFC = sol('BSFC')
    P_shaftmax = sol('P_{shaft-max}')
    P_shafttot = sol('P_{shaft-tot}')
    RPM = sol('RPM')
    V = sol('V')
    W_end = sol('W_{end}')

    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, eta_prop)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 1
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('prop efficiency')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('eta_propvsmission.png')

    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, BSFC)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 2
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('BSFC [lb/hp/hr]')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('BSFCvsmission.png')

    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, P_shaftmax)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 5
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('Max Available Shaft Power [hp]')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('P_shaftmaxvsmission.png')

    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, P_shafttot)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 5
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('Total Shaft Power [hp]')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('P_shafttotvsmission.png')

    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, RPM)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 9000
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('RPM')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('RPMvsmission.png')

    plt.close()
    fig, ax = plt.subplots()
    line, = ax.plot(t, W_end)
    ax.xaxis.set_ticks(np.arange(0,NSeg-1,1))
    labels = [item.get_text() for item in ax.get_xticklabels()]
    labels[mStart + np.round(NClimb1/2)] = 'Climb'
    labels[mEndClimb + np.round(NCruise1/2)] = 'Cruise'
    labels[mEndCruise + np.round(NClimb2/2)] = 'Climb'
    labels[mEndClimb2 + np.round(NLoiter/2)] = 'Loiter'
    labels[mEndLoiter + np.round(NCruise2/2)] = 'Cruise'
    ax.set_xticklabels(labels)
    ymin, ymax = 0, 150 
    ax.set_ylim([ymin, ymax])
    ax.set_ylabel('Aircraft Weight')
    ax.grid()
    ax.plot([mEndClimb-1,mEndClimb-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndCruise-1,mEndCruise-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndClimb2-1,mEndClimb2-1], [ymin, ymax], '--', color='r')
    ax.plot([mEndLoiter-1,mEndLoiter-1], [ymin, ymax], '--', color='r')
    plt.savefig('Wvsmission.png')
    
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
    plt.savefig('Vvsmission.png')