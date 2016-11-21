from numpy import pi, cos, tan
import numpy as np
from gpkit import Variable, Model, units
from gpkit.tools import te_exp_minus1
from gpkit.constraints.tight import TightConstraintSet as TCS
import matplotlib.pyplot as plt

sweep = False
n = 100
sweep_I = True
sweep_ro = False

plot = False

# Creating the model for the moment of inertia of a cylindrical cross-section
# Note: sweeps were created to check constraint tightness

class Beam(Model):
    def __init__(self):
        I     = Variable("I","m^4","Moment of inertia")
        A     = Variable("A", "m^2","Cross-sectional area")
        ro    = Variable("r_o" ,'m', 'Outer radius')
        romax = Variable('r_o_{max}','m', 'Maximum outer radius')

        constraints = [TCS([4*I**2/(A**2*ro**2) + A/pi <= ro**2]),
                        ro <= romax]

        # Minimizing cross-sectional area => min(weight)
        objective = A

        Model.__init__(self,objective,constraints)

if __name__ == '__main__':
    M = Beam()
    M.substitutions.update({'I':10**-4})
    M.substitutions.update({'r_o_{max}':0.2})
    if sweep:
        if sweep_I:
            M.substitutions.update({'I':('sweep',np.logspace(-10,-3,n))})
        if sweep_ro:
            M.substitutions.update({'r_o':('sweep',np.linspace(0.05,0.3,n))})
    sol = M.solve('mosek', skipsweepfailures = True)


if plot:
    if sweep_I:
        I = sol('I')
        A = sol('A')
        ro = sol('r_o')
        plt.semilogx(I,(4*I**2/(A**2*ro**2) + A/pi - ro**2)/(ro**2))
        plt.xlabel('Moment of inertia [m^-4]')
        plt.ylabel('Tightness of constraint')
        plt.title('Checking constraint tightness (normalized by ro**2)')
        plt.savefig('tightness.pdf')
        plt.show()
