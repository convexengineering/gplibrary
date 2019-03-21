"""Comparison of wing aerodynamics models."""

import os.path
from math import pi
import numpy as np
from matplotlib import pyplot as plt


def datcom_model(sweep, aspect, mach, c_la_M_2d=None):
    """DATCOM model for swept wing lift in compressible flow.

    See Section 4.1.3.2 of the DATCOM [1] and Eq. 12.6 of Raymer [2].

    Arguments:
        c_la_M_2d (scalar): Lift coefficient slope of the 2D airfoil used
            in the wing, at the given Mach number [units: radian**-1].
        sweep (scalar): The sweep angle of the wing, measured at mid-chord
            [units: raadians].
        aspect (scalar): Aspect ratio of the wing [units: dimensionless].
        mach (scalar): Freestream Mach number [units: dimensionless].

    Returns:
        scalar: The 3D lift coefficient slope of the wing (d C_L / d alpha)
            [units: radian**-1].

    References:
      [1] USAF Stability and Control DATCOM.
      [2] Raymer, Daniel. "Aircraft Design: a Conceptual Approach,"
          5th edition, AIAA, 2012.
    """
    if mach >= 1:
        raise ValueError('The model is not valid for supersonic flow')
    if sweep < 0 or sweep > pi/2:
        raise ValueError('Sweep angle must be in [0, pi / 2]')
    assert aspect > 0
    assert mach >= 0

    beta = (1 - mach**2)**0.5

    if c_la_M_2d is None:
        c_la_M_2d = 0.95 * 2 * pi / beta
    C_La_3d = 2 * pi * aspect / (
        2 + (4 + 4 * aspect**2 * pi**2 / c_la_M_2d**2
             * (1 + (np.tan(sweep) / beta)**2 ))**0.5)
    return C_La_3d


def drela_model_nosweep(aspect, mach, c_la_inc_2d):
    """Drela's model for straight wing lift in conpressible flow.

    See Drela [1], Eq. 8.87.

    Arguments:
        aspect (scalar): Aspect ratio of the wing [units: dimensionless].
        mach (scalar): Freestream Mach number [units: dimensionless].
        c_la_inc_2d (scalar): The incompressible, 2D lift coefficient slope
            of the airfoil used in the wing [units: radian**-1].

    Returns:
        scalar: The 3D lift coefficient slope of the wing (d C_L / d alpha)
            [units: radian**-1].

    References:
      [1] Drela, Mark. "Flight Vehicle Aerodynamics," MIT Press,
          2014.
    """
    if mach >= 1:
        raise ValueError('The model is not valid for supersonic flow')
    assert aspect > 0
    assert mach >= 0

    beta = (1 - mach**2)**0.5

    C_La_3d = c_la_inc_2d / (
        beta + c_la_inc_2d / (pi * aspect))

    return C_La_3d


def drela_model_infspan(sweep, mach, c_la_inc_2d):
    """Drela's model for lift of a swept, infinite-span wing in conpressible flow.

    See Drela [1], Eq. 8.103.

    Arguments:
        sweep (scalar): The sweep angle of the wing, measured at mid-chord
            [units: raadians].
        mach (scalar): Freestream Mach number [units: dimensionless].
        c_la_inc_2d (scalar): The incompressible, 2d lift coefficient slope
            of the airfoil used in the wing [units: radian**-1].

    Returns:
        scalar: The 3D lift coefficient slope of the wing (d C_L / d alpha)
            [units: radian**-1].

    References:
      [1] Drela, Mark. "Flight Vehicle Aerodynamics," MIT Press,
          2014.
    """
    if mach >= 1:
        raise ValueError('The model is not valid for supersonic flow')
    assert mach >= 0
    if sweep < 0 or sweep > pi/2:
        raise ValueError('Sweep angle must be in [0, pi / 2]')

    beta = (1 - mach**2)**0.5

    C_La_3d = c_la_inc_2d * np.cos(sweep) / (
        (beta * np.cos(sweep))**2 + np.sin(sweep)**2)**0.5

    return C_La_3d


def plot_mach():
    mach_datcom = np.linspace(0, 0.99)
    mach_drela = np.linspace(0, 0.99)
    CL_datcom = np.zeros(len(mach_datcom))
    CL_drela = np.zeros(len(mach_drela))
    c_la_inc_2d = 0.95 * 2 * pi

    aspects = [4, 8, 12]
    colors = ['C0', 'C1', 'C2']

    for aspect, color in zip(aspects, colors):
        for i in xrange(len(mach_datcom)):
            CL_datcom[i] = datcom_model(sweep=0, aspect=aspect, mach=mach_datcom[i])
            CL_drela[i] = drela_model_nosweep(
                aspect=aspect, mach=mach_drela[i], c_la_inc_2d=c_la_inc_2d)

        plt.plot(mach_datcom, CL_datcom,
                 color=color, linestyle='-',
                 label='DATCOM $\\Lambda=0$, AR={:.0f}'.format(aspect))
        plt.plot(mach_drela, CL_drela,
                 color=color, linestyle='--',
                 label='Drela $\\Lambda=0$, AR={:.0f}'.format(aspect))
    plt.xlabel(r'Freestream Mach number $M_\infty$ [-]')
    plt.ylabel(r'Lift slope $C_{L\alpha}$ [rad^-1]')
    plt.ylim([0, 14])
    plt.grid()
    plt.title('Lift slope vs. Mach Number for Unswept Wings')
    plt.legend()


def plot_aspect():
    aspect = np.linspace(2, 16)
    CL_datcom = np.zeros(len(aspect))
    CL_drela = np.zeros(len(aspect))
    c_la_inc_2d = 0.95 * 2 * pi

    machs = [0, 0.8]
    colors = ['C0', 'C1']

    for mach, color in zip(machs, colors):
        for i in xrange(len(aspect)):
            CL_datcom[i] = datcom_model(sweep=0, aspect=aspect[i], mach=mach)
            CL_drela[i] = drela_model_nosweep(
                aspect=aspect[i], mach=mach, c_la_inc_2d=c_la_inc_2d)

        plt.plot(aspect, CL_datcom,
                 color=color, linestyle='-',
                 label='DATCOM $\\Lambda=0$, $M_\\infty={:.1f}$'.format(mach))
        plt.plot(aspect, CL_drela,
                 color=color, linestyle='--',
                 label='Drela $\\Lambda=0$, $M_\\infty={:.1f}$'.format(mach))
    plt.xlabel(r'Wing aspect ratio $AR$ [-]')
    plt.ylabel(r'Lift slope $C_{L\alpha}$ [rad^-1]')
    plt.ylim([0, 9])
    plt.title('Lift slope vs. Aspect Ratio for Unswept Wings')
    plt.grid()
    plt.legend()


def plot_sweep():
    sweep = np.linspace(0, pi/2)
    CL_datcom = np.zeros(len(sweep))
    CL_drela = np.zeros(len(sweep))
    c_la_inc_2d = 0.95 * 2 * pi

    machs = [0, 0.8]
    colors = ['C0', 'C1']

    for mach, color in zip(machs, colors):
        for i in xrange(len(sweep)):
            CL_datcom[i] = datcom_model(sweep=sweep[i], aspect=1e3, mach=mach)
            CL_drela[i] = drela_model_infspan(
                sweep=sweep[i], mach=mach, c_la_inc_2d=c_la_inc_2d)

        plt.plot(np.rad2deg(sweep), CL_datcom,
                 color=color, linestyle='-',
                 label='DATCOM $M_\\infty={:.1f}$'.format(mach))
        plt.plot(np.rad2deg(sweep), CL_drela,
                 color=color, linestyle='--', marker='+',
                 label='Drela $M_\\infty={:.1f}$'.format(mach))
    plt.xlabel(r'Wing sweep angle $\Lambda$ [deg]')
    plt.ylabel(r'Lift slope $C_{L\alpha}$ [rad^-1]')
    plt.ylim([0, 10.5])
    plt.title('Lift slope vs. Sweep for Infinite-Span Wings')
    plt.grid()
    plt.legend()


def main():
    fig_dir = 'docs/figures/wing_lift_model_compare'
    plt.figure(figsize=(8, 6))
    plot_mach()
    plt.savefig(os.path.join(fig_dir, 'CLa_vs_mach.png'))

    plt.figure(figsize=(8, 6))
    plot_aspect()
    plt.savefig(os.path.join(fig_dir, 'CLa_vs_aspect.png'))


    plt.figure(figsize=(8, 6))
    plot_sweep()
    plt.savefig(os.path.join(fig_dir, 'CLa_vs_sweep.png'))

    plt.show()


if __name__ == '__main__':
    main()
