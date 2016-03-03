from numpy import pi
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

class GasPoweredHALE(Model):
    def setup(self):

        # define number of segments
        Nseg = 3

        constraints = []

        #----------------------------------------------------
        # Weight and fuel model

        MTOW = Variable('MTOW', 'lbf', 'max take off weight')
        W_end = VectorVariable(Nseg, 'W_{end}', 'lbf', 'segment-end weight')
        W_fuel = VectorVariable(Nseg, 'W_{fuel}', 'lbf',
                                'segment-fuel weight') 
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        W_payload = Variable('W_{fix}',10,'lbf', 'Payload weight')
        W_avionics = Variable('W_{avionics}', 2, 'lbf', 'Avionics weight')
        f_airframe = Variable('f_{airframe}', 0.3, '-', 'airframe weight fraction')
        W_airframe = Variable('W_{airframe}', 'lbf', 'airframe weight')
        W_begin = W_end.left # define beginning of segment weight
        W_begin[0] = MTOW 

        # end of first segment weight + first segment fuel weight must be greater 
        # than MTOW.  Each end of segment weight must be greater than the next end
        # of segment weight + the next segment fuel weight. The last end segment
        # weight must be greater than the zero fuel weight
        constraints.extend([MTOW >= W_end[0] + W_fuel[0], 
                            W_end[:-1] >= W_end[1:] + W_fuel[1:], 
                            W_end[-1] >= W_zfw,
                            W_airframe >= f_airframe*MTOW,
                            W_zfw >= W_payload + W_avionics + W_airframe])
        
        #----------------------------------------------------
        # Steady level flight model

        CD = Variable('C_D', 0.05, '-', 'Drag coefficient')
    	CL = Variable('C_L', 1, '-', 'Lift coefficient')
        V = VectorVariable(Nseg,'V', 'm/s','cruise speed')
        rho = VectorVariable(Nseg, r'\rho', [1.2, 1.2, 1.2], 'kg/m^3', 'air density')
        S = Variable('S', 190, 'ft^2', 'wing area')
        eta_prop = VectorVariable(Nseg, r'\eta_{prop}', [0.6, 0.8, 0.6], '-',
                                  'propulsive efficiency')
        P_shaft = VectorVariable(Nseg, 'P_{shaft}', [30,30,30], 'hp', 'Shaft power')

        constraints.extend([P_shaft >= V*(W_end+W_begin)/2*CD/CL/eta_prop,
                            0.5*rho*CL*S*V**2 == (W_end*W_begin)**0.5])

        #----------------------------------------------------
        # Breguet Range
        z_bre = VectorVariable(Nseg, 'z_{bre}', '-', 'breguet coefficient')
        BSFC = Variable('BSFC', 0.5, 'lbf/hr/hp', 'brake specific fuel consumption')
        t = VectorVariable(Nseg-1, 't', [5,0.2], 'days', 'time on station')
        R = Variable('R', 100, 'nautical_miles', 'range to station')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([z_bre[1:] >= V[1:]*t[1:]*BSFC*CD/CL/eta_prop[1:],
                            z_bre[0] >= R*BSFC*CD/CL/eta_prop[0],
                            W_fuel/W_end >= te_exp_minus1(z_bre, 3)])

        objective = MTOW
        return objective, constraints

if __name__ == "__main__":
	M = GasPoweredHALE()
	M.solve()
