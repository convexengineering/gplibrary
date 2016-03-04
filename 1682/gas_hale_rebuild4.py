from numpy import pi
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

class GasPoweredHALE(Model):
    def setup(self):

        # Ok I disagree with how you have this set up. I think it's a good 
        # to break this down.  But I think the climb constraint can be built
        # into the flight segments.  I propose this: let's have 1 segement 
        # cruise out.  Let's have 5 segments on loiter to improve fidelity.
        # and lets have another cruise segment coming back.  And then we 
        # build the climb constraints using the beginning weights of the 
        # first section when we get up to cruise altitude and then the 
        # beginning of the second section when we climb to get up to station.


        # define number of segments
        NSeg = 7 # number of flight segments
        NCruise = 2 # number of cruise segments
        NLoiter = NSeg - NCruise # number of loiter segments
        iCruise = [0,-1] # cuise index
        iLoiter = [1] # loiter index
        for i in range(2,NSeg): iLoiter.append(i)
        iClimb = [0,1] # climb index

        #Note: NSeg has to be an odd number
        # defining indices of different flight segments
        #NLoiter = (NSeg-1)/2
        #if NSeg == 3:
        #    NCruise = [0,2]
        #elif NSeg == 7:
        #    Nclimb = [0,2,4,6]
        #    NCruise = [1,5]
        
        constraints = []

        #----------------------------------------------------
        # Fuel weight model 

        MTOW = Variable('MTOW', 'lbf', 'max take off weight')
        W_end = VectorVariable(NSeg, 'W_{end}', 'lbf', 'segment-end weight')
        W_fuel = VectorVariable(NSeg, 'W_{fuel}', 'lbf',
                                'segment-fuel weight') 
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        W_pay = Variable('W_{pay}',10,'lbf', 'Payload weight')
        W_avionics = Variable('W_{avionics}', 2, 'lbf', 'Avionics weight')
        f_airframe = Variable('f_{airframe}', 0.25, '-', 'airframe weight fraction')
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
                            W_airframe >= f_airframe*MTOW])
        
        #----------------------------------------------------
        # Steady level flight model

        CD = VectorVariable(NSeg, 'C_D', '-', 'Drag coefficient')
    	CL = VectorVariable(NSeg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(NSeg, 'V', 'm/s','cruise speed')
        rho = VectorVariable(NSeg, r'\rho', 'kg/m^3', 'air density')
        S = Variable('S', 16, 'ft^2', 'wing area')
        eta_prop = VectorVariable(NSeg, r'\eta_{prop}', np.linspace(0.8,0.8,NSeg), '-',
                                  'propulsive efficiency')
        P_shaft = VectorVariable(NSeg, 'P_{shaft}', 'hp', 'Shaft power')

        # Climb model
        h_dot = Variable(NSeg, 'h_{dot}', 200, 'ft/min', 'Climb rate')
        
        # Berk I kinda of disagree with the way you set up the climb model.  
        # The way you have it set up, you are combining the cruise and climb 
        # conditions.  I believe those need to be separate. Also this is 
        # saying that each flight segment has to also meet a climb constaint
        # If you look at what I did I imposed the climb condition hdot <= 
        # (T-D)*V/W for the first and second leg which will probably be the 
        # hardest.  

        constraints.extend([P_shaft >= V*(W_end+W_begin)/2*CD/CL/eta_prop, # + W_begin*h_dot/eta_prop,
                            W_begin[iClimb]*CD[iClimb]/CL[iClimb] >= h_dot*W_begin[iClimb]/V[iClimb] + 
                                                      0.5*rho[iClimb]*V[iClimb]**2*S*CD[iClimb],
                            0.5*rho*CL*S*V**2 >= (W_end+W_begin)/2])

        #----------------------------------------------------
        # Engine Model
        W_eng = Variable('W_{eng}', 'lbf', 'Engine weight')
        W_engtot = Variable('W_{eng-tot}', 'lbf', 'Installed engine weight')
        W_engref = Variable('W_{eng-ref}', 4.4107, 'lbf', 'Reference engine weight')
        P_shaftref = Variable('P_{shaft-ref}', 2.295, 'hp', 'reference shaft power')

        # Engine Weight Constraints
        constraints.extend([W_eng/W_engref >= 0.5538*(P_shaft/P_shaftref)**1.075,
                            W_engtot >= 2.572*W_eng**0.922*units('lbf')**0.078])

        #----------------------------------------------------
        # Weight breakdown
        constraints.extend([W_airframe >= f_airframe*MTOW,
                            W_zfw >= W_pay + W_avionics + W_airframe + W_engtot])
        
        #----------------------------------------------------
        # Breguet Range
        z_bre = VectorVariable(NSeg, 'z_{bre}', '-', 'breguet coefficient')
        BSFC = VectorVariable(NSeg,'BSFC', np.linspace(0.5,0.5,NSeg), 'lbf/hr/hp', 'brake specific fuel consumption')
        t = VectorVariable(NSeg, 't', 'days', 'time on station')
        t_cruise = Variable('t_{cruise}', 0.5, 'days', 'time to station')
        dt = Variable('dt', 5/NLoiter, 'days', 'time interval for loiter')
        R = Variable('R', 200, 'nautical_miles', 'range to station')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([z_bre >= V*t*BSFC*CD/CL/eta_prop,
                            R == V[iCruise]*t[iCruise],
                            t[iLoiter] == dt,
                            t[iCruise] <= t_cruise,
                            W_fuel/W_end >= te_exp_minus1(z_bre, 3)])

        #----------------------------------------------------
        # Aerodynamics model

        Cd0 = Variable('C_{d0}', 0.02, '-', 'Non-wing drag coefficient')
        CLmax = Variable('C_{L-max}', 1.5, '-', 'Maximum lift coefficient')
        e = Variable('e', 0.9, '-', 'Spanwise efficiency')
        AR = Variable('AR', '-', 'Aspect ratio')
        b = Variable('b', 'ft', 'Span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', 'Dynamic viscosity')
        Re = VectorVariable(NSeg, 'Re', '-', 'Reynolds number')
        Cf = VectorVariable(NSeg, 'C_f', '-', 'wing skin friction coefficient')
        Kwing = Variable('K_{wing}', 1.3, '-', 'wing form factor')
        cl_16 = Variable('cl_{16}', 0.0001, '-', 'profile stall coefficient')

        constraints.extend([CD >= Cd0 + 2*Cf*Kwing + CL**2/(pi*e*AR) + cl_16*CL**16,
                            b**2 == S*AR,
                            AR <= 20, # temporary constraint until we input a valid structural model
                            CL <= CLmax, 
                            Re == rho*V/mu*(S/AR)**0.5,
                            Cf >= 0.074/Re**0.2])

        #----------------------------------------------------
        # Atmosphere model

        h = VectorVariable(NSeg, 'h', 'ft', 'Altitude')
        gamma = Variable(r'\gamma',1.4,'-', 'Heat capacity ratio of air')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        T_sl = Variable('T_{sl}', 288.15, 'K', 'Temperature at sea level')
        L_atm = Variable('L_{atm}', 0.0065, 'K/m', 'Temperature lapse rate')
        T_atm = VectorVariable(NSeg, 'T_{atm}', 'K', 'Air temperature')
        a_atm = VectorVariable(NSeg, 'a_{atm}','m/s', 'Speed of sound at altitude')
        R_spec = Variable('R_{spec}', 287.058,'J/kg/K', 'Specific gas constant of air')
        TH = (g/R_spec/L_atm).value.magnitude  # dimensionless

        # ok so I issue I have with this one is that we are subjecting ourselves
        # to some altitude for cruise going out to station.  It may want to fly
        # lower than this, so we should let it tell us what it wants to fly at. 

        constraints.extend([#h <= [20000, 20000, 20000]*units.m,  # Model valid to top of troposphere
                            T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                            rho == p_sl*T_atm**(TH-1)/R_spec/(T_sl**TH),
                            #h[NLoiter] >= 15000*units('ft'), #makes sure that the loiter occurs above minimum h
                            #h[NCruise] >= h_min
                            ])
            # http://en.wikipedia.org/wiki/Density_of_air#Altitude

        #----------------------------------------------------
        # altitude constraints
        h_station = Variable('h_{station}', 15000, 'ft', 'minimum altitude at station')
        h_min = Variable('h_{min}', 1000, 'ft', 'minimum cruise altitude')

        constraints.extend([h[iLoiter] >= h_station,
                            h[iCruise] >= h_min])

        #----------------------------------------------------
        # wind speed model

        # So this wind model is based on the 95%. Which from what Oren is telling us
        # is probably overkill.  I think it would be better to break up the wind 
        # speed by directly specifying the wind speeds over time.  That way we can vary
        # it and say well, if its only 30 m/s one day of the 5 we're up there this is 
        # what happens to our endurance. 

        V_wind = VectorVariable(NLoiter, 'V_{wind}', [20,25,30,10,15], 'm/s', 'wind speed')
        #wd_cnst = Variable('wd_{cnst}', 0.0015, 'm/s/ft', 
        #                   'wind speed constant predicted by model')
        #                    #0.002 is worst case, 0.0015 is mean at 45d
        #wd_ln = Variable('wd_{ln}', 8.845, 'm/s',
        #                 'linear wind speed variable')
        #                #13.009 is worst case, 8.845 is mean at 45deg
        #h_min = Variable('h_{min}', 11800, 'ft', 'minimum height')
        #h_max = Variable('h_{max}', 20866, 'ft', 'maximum height')

        constraints.extend([#V_wind >= wd_cnst*h + wd_ln, 
                            V[iLoiter] >= V_wind])

        objective = MTOW
        return objective, constraints

if __name__ == '__main__':
	M = GasPoweredHALE()
	M.solve()
