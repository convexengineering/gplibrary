"""Demonstrate the fuselage models."""

import numpy as np
from matplotlib import pyplot as plt
import gpkit
from gpkit import units, Variable, Model
from gpkit.tools.autosweep import autosweep_1d
from gpkit.small_scripts import mag

import fuselage


def main():
    fuse_model = fuselage.FuselageAeroSubsonic()

    fuse_model.substitutions['diam_fuse'] = 0.1    # [units: meter]
    fuse_model.substitutions['len_fuse'] = 0.5    # [units: meter]
    fuse_model.substitutions['mach'] = 0.8    # [units: dimensionless]
    fuse_model.substitutions['p_static'] = 30e3    # [units: pascal]

    fuse_model.cost = fuse_model['drag_fuse']

    sol = fuse_model.solve()
    print sol.table()

    sweep = fuse_model.autosweep({fuse_model['mach']: (0.1, 0.99)}, tol=1e-3)
    sweep.plot(fuse_model['drag_fuse'])
    sweep.plot(fuse_model['C_D0_fuse'])
    plt.show()

if __name__ == '__main__':
    main()
