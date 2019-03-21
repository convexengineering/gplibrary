"""Demonstrate the wing model."""

import os.path
from subprocess import CalledProcessError
import numpy as np
from matplotlib import pyplot as plt
import gpkit
from gpkit import units, Variable, Model
from gpkit.tools.autosweep import autosweep_1d
from gpkit.small_scripts import mag

import wing
from wing_model_compare import datcom_model
import airfoil
import flight_state

"""
Issues with the model 2019-03-06 T 1223
Issue 1 maximuzing lift pushes sweep to 0, but this makes the
problem dual infeasible. Need to keep sweep at a small value.

Also there is some issue with beta_pg
 - solved with signomial constraint
"""

def tan2_taylor(x, order=8):
    """Taylor series of tan^2(x)."""
    coeffs = [0, 0, 1, 0, 2/3., 0, 17/45., 0, 62/315., 0, 1382/14175.,
              0, 21844/467775.]
    if order + 1 > len(coeffs):
        raise ValueError('I do not have coeffs up to that order.')
    val = sum([coeffs[i] * x**i for i in range(order + 1)])
    return val


def lift_plots():
    fig_dir = 'docs/figures/wing_lift_model_compare'

    aspect = 10

    flight_state_model = flight_state.FlightState()
    airfoil_model = airfoil.DummyAirfoil(flight_state_model)
    wing_model = wing.WingAeroSubsonic(flight_state_model, airfoil_model, limit_sweep_to_valid_range=False)

    wing_model.substitutions['S_ref'] = 1    # [units: meter^2]
    wing_model.substitutions['aspect'] = aspect    # [units: dimensionless]
    wing_model.substitutions['mach'] = 0.5    # [units: dimensionless]
    wing_model.substitutions['p_static'] = 30e3    # [units: pascal]
    # wing_model.substitutions['alpha'] = np.deg2rad(2)    # [units: radian]
    wing_model.substitutions['sweep_c2'] = 1e-6    # [units: dimensionless]
    wing_model.substitutions['C_D'] = 1e3    # [units: dimensionless]

    # maximize lift
    wing_model.cost = wing_model['lift']**-1

    sol = wing_model.solve()
    print sol.table()

    ### Compare the GPkit implementation to my non-GP implementation
    # of the DATCOM wing lift model###

    # Do an autosweep of the GP model vs. Mach number.
    # Confusing! GP uses "sweep" to refer to sweeping the model over
    # different values for a variable. But "sweep_c2" is also a
    # variable in the wing model.
    # varsweep = wing_model.autosweep({wing_model['mach_WingAeroSubsonic']: (0.01, 0.99)}, tol=1e-3)

    # Solve the GP at a few mach numbers
    mach_short = np.linspace(0.01, 0.80, 10)
    C_La_gp = np.zeros(len(mach_short))
    beta_gp = np.zeros(len(mach_short))
    for i in xrange(len(mach_short)):
        print('mach : {}'.format(mach_short[i]))
        wing_model.substitutions['mach'] = mach_short[i]
        sol = wing_model.solve()
        C_La_gp[i] = sol['variables']['C_La']
        beta_gp[i] = sol['variables']['beta_pg']
        # assert sol['variables']['sweep_c2'] < 1e-6

    # Evaluate the non-GP implementation of the DATCOM model
    # vs. Mach number
    # Set sweep to 0, because the GP model should bring sweep to 0 if
    # maximizing lift.
    mach = np.linspace(0.01, 0.99)
    C_La_datcom = np.array([datcom_model(0, aspect, m) for m in mach])

    # Plot lift coefficient slope vs. Mach number
    # plt.plot(varsweep['mach'], varsweep['C_La'], color='C0', linestyle='-', label='GPkit impl. of DATCOM model (autosweep)')
    plt.plot(mach_short, C_La_gp, color='C1', linestyle='', marker='+', markersize=12,
             label='GPkit impl. of DATCOM model')
    plt.plot(mach, C_La_datcom, color='black', linestyle='-.', label='DATCOM model')
    plt.legend()
    plt.xlabel('Freestream Mach number [-]')
    plt.ylabel(r'Lift slope $C_{L\alpha}$ [rad^-1]')
    plt.title('Verification of GP impl. of DATCOM wing lift model'
        + '\n at wing sweep $\\Lambda = 0$')
    plt.savefig(os.path.join(fig_dir, 'GP_mach.png'))

    # Plot the Prandtl-Glauert beta vs. Mach number
    plt.figure()
    # plt.plot(varsweep['mach'], varsweep['beta_pg'], label='GPkit impl. (autosweep)')
    plt.plot(mach_short, beta_gp, color='C1', linestyle='', marker='+', markersize=12,
             label='GPkit impl.')
    plt.plot(mach, (1 - mach**2)**0.5,
             color='black', linestyle='-.', label=r'$\beta = \sqrt{1 - M^2}$')
    plt.xlabel('Freestream Mach number [-]')
    plt.ylabel(r'Prandtl-Glauert $\beta$ [-]')
    plt.legend()
    plt.ylim([0, 1.1])

    ### Compare implementations vs. sweep ###
    mach = 0.6
    wing_model.substitutions['mach'] = mach
    # Solve the GP at a few sweeps
    sweep_short = np.linspace(0.01, np.pi/2 - 1e-3, 10)
    C_La_gp = np.zeros(len(sweep_short))
    beta_gp = np.zeros(len(sweep_short))
    for i in xrange(len(sweep_short)):
        wing_model.substitutions['sweep_c2'] = sweep_short[i]
        sol = wing_model.solve()
        C_La_gp[i] = sol['variables']['C_La']
        beta_gp[i] = sol['variables']['beta_pg']

    # Evaluate the non-GP implementation of the DATCOM model
    # vs. Mach number
    # Set sweep to 0, becuase the GP model should bring sweep to 0 if
    # maximizing lift.
    sweep = np.linspace(0.0, np.pi / 2)
    C_La_datcom = np.array([datcom_model(s, aspect, mach) for s in sweep])

    # Plot lift coefficient slope vs. Mach number
    plt.figure(figsize=(6, 8))
    plt.subplot(2, 1, 1)
    plt.plot(np.rad2deg(sweep_short), C_La_gp, color='C1', linestyle='', marker='+', markersize=12,
             label='GPkit impl. of DATCOM model')
    plt.plot(np.rad2deg(sweep), C_La_datcom, color='black', linestyle='-.', label='DATCOM model')
    plt.legend()
    plt.xlabel(r'Wing sweep angle $\Lambda$ [deg]')
    plt.ylabel(r'Lift slope $C_{L\alpha}$ [rad^-1]')
    plt.title('Verification of GP impl. of DATCOM wing lift model'
              + '\n at freestream Mach $M_\\infty={:.1f}$'.format(mach))
    plt.axvline(x=70, color='grey')

    plt.subplot(2, 1, 2)
    plt.plot(np.rad2deg(sweep[:-1]), np.tan(sweep[:-1])**2,
             color='black', linestyle='-.', label=r'$\tan^2(\Lambda)$')
    order = 12
    approx = np.array([tan2_taylor(s, order) for s in sweep])
    plt.plot(np.rad2deg(sweep), approx,
             color='C1', linestyle='-', label='Taylor series, order {:d}'.format(order))
    plt.xlabel(r'Wing sweep angle $\Lambda$ [deg]')
    plt.ylim([0, 10])
    plt.text(70, 3, r'Model limit max $\Lambda$',
             rotation=90, verticalalignment='center', horizontalalignment='left')

    plt.legend()
    plt.axvline(x=70, color='grey')
    plt.tight_layout()
    plt.subplots_adjust(top=0.9)
    plt.savefig(os.path.join(fig_dir, 'GP_sweep.png'))


def drag_plots():
    aspect = 10
    C_L = 0.3

    flight_state_model = flight_state.FlightState()
    airfoil_model = airfoil.DummyAirfoil(flight_state_model)
    wing_model = wing.WingAeroSubsonic(flight_state_model, airfoil_model, limit_sweep_to_valid_range=False)

    wing_model.substitutions['S_ref'] = 1    # [units: meter^2]
    wing_model.substitutions['aspect'] = aspect    # [units: dimensionless]
    wing_model.substitutions['p_static'] = 30e3    # [units: pascal]
    wing_model.substitutions['C_L'] = C_L    # [units: dimensionless]

    # minimize drag
    wing_model.cost = wing_model['drag']

    # Plot drag vs Mach at various sweeps
    ax1 = plt.subplot(3, 1, 1)
    plt.subplot(3, 1, 2, sharex=ax1)
    plt.subplot(3, 1, 3, sharex=ax1)

    mach = np.linspace(0.5, 0.99, 20)
    C_D = np.zeros(len(mach))
    c_l = np.zeros(len(mach))
    alpha = np.zeros(len(mach))
    sweeps = np.deg2rad([1e-6, 30, 60])
    # sweeps = np.deg2rad([1e-6, 20, 40])
    for sweep in sweeps:
        for i in range(len(mach)):
            wing_model.substitutions['mach'] = mach[i]
            wing_model.substitutions['sweep_c2'] = sweep
            try:
                sol = wing_model.solve()
                C_D[i] = sol['variables']['C_D']
                c_l[i] = sol['variables']['c_l']
                alpha[i] = sol['variables']['alpha']
            # except CalledProcessError, RuntimeWarning:
            except:
                C_D[i] = np.nan
                c_l[i] = np.nan
                alpha[i] = np.nan
        plt.subplot(3, 1, 1)
        plt.plot(mach, C_D, label='$\\Lambda = ${:.0f} deg'.format(
            np.rad2deg(sweep)))
        plt.subplot(3, 1, 2)
        plt.plot(mach, c_l, label='$\\Lambda = ${:.0f} deg'.format(
            np.rad2deg(sweep)))
        plt.subplot(3, 1, 3)
        plt.plot(mach, np.rad2deg(alpha), label='$\\Lambda = ${:.0f} deg'.format(
            np.rad2deg(sweep)))
    plt.subplot(3, 1, 1)
    plt.ylim([0, 0.025])
    plt.ylabel('$C_D$ [-]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')
    plt.legend()
    plt.title('Drag vs Mach at constant $C_L =${:.2f}'.format(C_L))

    plt.subplot(3, 1, 2)
    plt.ylabel('$c_l$ [-]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')
    plt.legend()

    plt.subplot(3, 1, 3)
    plt.ylabel('Angle of attack $\\alpha$ [deg]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')
    plt.legend()

    plt.tight_layout()
    plt.subplots_adjust(hspace=0)
    """
    Why does C_D (at const C_L) decrease with sweep at low mach?
    I think because of cos^3 term on pressure drag in Drela 8.183?
    Also, although Drela says that drag should increase monotonically with sweep
    below the critical Mach number, the plots in Hoerner indicate that
    C_D may actually decrease slightly with sweep even at lower (~0.6)
    Mach numbers.
    """


def ld_plots():
    aspect = 10

    flight_state_model = flight_state.FlightState()
    airfoil_model = airfoil.DummyAirfoil(flight_state_model)
    wing_model = wing.WingAeroSubsonic(flight_state_model, airfoil_model)

    wing_model.substitutions['S_ref'] = 1    # [units: meter^2]
    wing_model.substitutions['aspect'] = aspect    # [units: dimensionless]
    wing_model.substitutions['mach'] = 0.5    # [units: dimensionless]
    wing_model.substitutions['p_static'] = 30e3    # [units: pascal]
    wing_model.substitutions['sweep_c2'] = 1e-6    # [units: dimensionless]
    # wing_model.substitutions['lift'] = 1e3    # [units: newton]

    # maximize lift/drag
    # Also, add a light cost to sweep - there are structural and control
    # penalties on sweep which this model does not capture.
    # Without this cost on sweep, the model would optimize to high sweep
    # angles for low-Mach wings, because the drag model predicts a very slight
    # decrease in drag from sweep even at low Mach.
    wing_model.cost = wing_model['drag']/wing_model['lift'] + 0.005 * wing_model['sweep_c2']

    sol = wing_model.solve()
    print sol.table()

    # Plot sweep for best L/D vs mach
    del wing_model.substitutions['sweep_c2']
    mach = np.linspace(0.01, 0.99, 11)
    ld = np.zeros(len(mach))
    c_l = np.zeros(len(mach))
    sweep = np.zeros(len(mach))
    gp_plot_style = {'color':'C1', 'linestyle':'', 'marker':'+', 'markersize':12}

    for i in range(len(mach)):
        wing_model.substitutions['mach'] = mach[i]
        sol = wing_model.localsolve()
        ld[i] = sol['variables']['lift'] / sol['variables']['drag']
        c_l[i] = sol['variables']['c_l']
        sweep[i] = sol['variables']['sweep_c2']

    plt.subplot(3, 1, 1)
    plt.plot(mach, ld, label='GP model', **gp_plot_style)
    plt.ylabel('$L/D$')

    plt.subplot(3, 1, 2)
    plt.plot(mach, np.rad2deg(sweep), label='GP model', **gp_plot_style)
    plt.ylabel('Sweep $\\Lambda$ [deg]')

    plt.subplot(3, 1, 3)
    plt.plot(mach, c_l, label='GP model', **gp_plot_style)
    plt.ylabel('2D lift coeff. $c_l$ [-]')
    plt.xlabel('Freestream Mach number $M_\\infty$ [-]')

    plt.suptitle('Optimizing for best $L/D$')
    plt.tight_layout()


def main():
    # lift_plots()
    # drag_plots()
    ld_plots()
    plt.show()

if __name__ == '__main__':
    main()
