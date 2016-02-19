from gas_hale import GasPoweredHALE
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

# to modify the GasPoweredHALE model, edit gas_hale.py

m = GasPoweredHALE()
m.substitutions.update({'t': 6})
sol = m.solve()

x = sol('W_{fuel}')/sol('W')
y = sol('W_{pay}')/sol('W')
z = sol('P_{shaft}')/745.699872
print x
print y
print z

import numpy as np
m.substitutions.update({"t": ('sweep', np.linspace(4, 7, 5)),
                        "d_{footprint}": ('sweep', np.linspace(0.1,0.4,5))})
sol = m.solve(solver="cvxopt", verbosity=0, skipsweepfailures=True)

import numpy as np
m.substitutions.update({"t": ('sweep', np.linspace(4, 7, 5)),
                        "d_{footprint}": ('sweep', np.linspace(0.1,0.4,5))})
sol = m.solve(solver="cvxopt", verbosity=0, skipsweepfailures=True)

#%matplotlib inline
from gpkit.interactive.plotting import contour_array
_ = contour_array(m, "f_{airframe}", "t", ["W", "S"])