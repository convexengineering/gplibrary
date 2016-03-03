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
        W_pay = Variable('W_{pay}',10,'lbf', 'Payload weight')
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
                            W_zfw >= W_pay + W_avionics + W_airframe])
        
        #----------------------------------------------------
        # Steady level flight model

        CD = VectorVariable(Nseg, 'C_D', '-', 'Drag coefficient')
    	CL = VectorVariable(Nseg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(Nseg, 'V', 'm/s','cruise speed')
        rho = VectorVariable(Nseg, r'\rho', 'kg/m^3', 'air density')
        S = Variable('S', 190, 'ft^2', 'wing area')
        eta_prop = VectorVariable(Nseg, r'\eta_{prop}', [0.6, 0.8, 0.6], '-',
                                  'propulsive efficiency')
        P_shaft = VectorVariable(Nseg, 'P_{shaft}', 'hp', 'Shaft power')

        constraints.extend([P_shaft >= V*(W_end+W_begin)/2*CD/CL/eta_prop,
                            0.5*rho*CL*S*V**2 >= (W_end+W_begin)/2])

        #----------------------------------------------------
        # Breguet Range
        z_bre = VectorVariable(Nseg, 'z_{bre}', '-', 'breguet coefficient')
        BSFC = Variable('BSFC', 0.5, 'lbf/hr/hp', 'brake specific fuel consumption')
        t = VectorVariable(Nseg, 't', [0.2,5,0.2], 'days', 'time on station')
        R = Variable('R', 100, 'nautical_miles', 'range to station')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        
        constraints.extend([z_bre >= V*t*BSFC*CD/CL/eta_prop,
                            R == V[0]*t[0],
                            W_fuel/W_end >= te_exp_minus1(z_bre, 3)])

        #----------------------------------------------------
        # Aerodynamics model

        Cd0 = Variable('C_{d0}', 0.02, '-', 'Non-wing drag coefficient')
        CLmax = Variable('C_{L-max}', 1.5, '-', 'Maximum lift coefficient')
        e = Variable('e', 0.9, '-', 'Spanwise efficiency')
        AR = Variable('AR', '-', 'Aspect ratio')
        b = Variable('b', 'ft', 'Span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', 'Dynamic viscosity')
        Re = VectorVariable(Nseg, 'Re', '-', 'Reynolds number')
        Cf = VectorVariable(Nseg, 'C_f', '-', 'wing skin friction coefficient')
        Kwing = Variable('K_{wing}', 1.3, '-', 'wing form factor')
        cl_16 = Variable('cl_{16}', 0.0001, '-', 'profile stall coefficient')

        constraints.extend([CD >= Cd0 + 2*Cf*Kwing + CL**2/(pi*e*AR) + cl_16*CL**16,
                            b**2 == S*AR,
                            CL <= CLmax, 
                            Re == rho*V/mu*(S/AR)**0.5,
                            Cf >= 0.074/Re**0.2])

        #----------------------------------------------------
        # Atmosphere model

        h = VectorVariable(Nseg, 'h', [8000,15000,8000], 'ft', 'Altitude')
        gamma = Variable(r'\gamma',1.4,'-', 'Heat capacity ratio of air')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        T_sl = VectorVariable(Nseg, 'T_{sl}', [288.15,288.15,288.15], 'K',
                              'Temperature at sea level')
        L_atm = Variable('L_{atm}', 0.0065, 'K/m', 'Temperature lapse rate')
        T_atm = VectorVariable(Nseg, 'T_{atm}', 'K', 'Air temperature')
        a_atm = VectorVariable(Nseg, 'a_{atm}','m/s', 'Speed of sound at altitude')
        R_spec = Variable('R_{spec}', 287.058,'J/kg/K', 'Specific gas constant of air')
        TH = (g/R_spec/L_atm).value.magnitude  # dimensionless

        constraints.extend([#h <= [20000, 20000, 20000]*units.m,  # Model valid to top of troposphere
                            T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                            rho == p_sl*T_atm**(TH-1)/R_spec/(T_sl**TH)])
            # http://en.wikipedia.org/wiki/Density_of_air#Altitude

        objective = MTOW
        return objective, constraints

if __name__ == '__main__':
	M = GasPoweredHALE()
	M.solve()
