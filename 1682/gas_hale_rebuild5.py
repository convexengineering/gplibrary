from numpy import pi
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

class GasPoweredHALE(Model):
    def setup(self):

        # define number of segments
        NSeg = 9 # number of flight segments
        NCruise = 2 # number of cruise segments
        NClimb = 2 # number of climb segments
        NLoiter = NSeg - NCruise - NClimb# number of loiter segments
        iCruise = [1,-1] # cuise index
        iLoiter = [3] # loiter index
        for i in range(4,NSeg-1): iLoiter.append(i)
        iClimb = [0,2] # climb index

        constraints = []

        #----------------------------------------------------
        # Fuel weight model 

        MTOW = Variable('MTOW', 'lbf', 'max take off weight')
        W_end = VectorVariable(NSeg, 'W_{end}', 'lbf', 'segment-end weight')
        W_fuel = VectorVariable(NSeg, 'W_{fuel}', 'lbf',
                                'segment-fuel weight') 
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        f_airframe = Variable('f_{airframe}', 0.25, '-', 'airframe weight fraction')
        W_airframe = Variable('W_{airframe}', 'lbf', 'airframe weight')
        W_begin = W_end.left # define beginning of segment weight
        W_begin[0] = MTOW 

        # Payload model
        W_pay = Variable('W_{pay}',10,'lbf', 'Payload weight')
        Vol_pay = Variable('Vol_{pay}',0.5,'ft^3','Payload volume')

        # Avionics model
        W_avionics = Variable('W_{avionics}', 2, 'lbf', 'Avionics weight')
        Vol_avionics = Variable('Vol_{avionics}',0.125,'ft^3','Avionics volume')

        # end of first segment weight + first segment fuel weight must be greater 
        # than MTOW.  Each end of segment weight must be greater than the next end
        # of segment weight + the next segment fuel weight. The last end segment
        # weight must be greater than the zero fuel weight

        constraints.extend([MTOW >= W_end[0] + W_fuel[0], 
                            W_end[:-1] >= W_end[1:] + W_fuel[1:], 
                            W_end[-1] >= W_zfw])
        
        #----------------------------------------------------
        # Steady level flight model

        CD = VectorVariable(NSeg, 'C_D', '-', 'Drag coefficient')
    	CL = VectorVariable(NSeg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(NSeg, 'V', 'm/s','cruise speed')
        rho = VectorVariable(NSeg, r'\rho', 'kg/m^3', 'air density')
        S = Variable('S', 'ft^2', 'wing area')
        eta_prop = VectorVariable(NSeg, r'\eta_{prop}', '-',
                                  'propulsive efficiency')
        P_shaft = VectorVariable(NSeg, 'P_{shaft}', 'hp', 'Shaft power')
        T = VectorVariable(NSeg, 'T', 'lbf', 'Thrust')

        # Climb model
        h_dot = Variable('h_{dot}', 120, 'ft/min', 'Climb rate')
        
        constraints.extend([P_shaft >= T*V/eta_prop,
                            T >= 0.5*rho*V**2*CD*S,
                            T[iClimb] >= 0.5*rho[iClimb]*V[iClimb]**2*CD[iClimb]*S + W_begin[iClimb]*h_dot/V[iClimb],
                            0.5*rho*CL*S*V**2 >= (W_end+W_begin)/2,
                            eta_prop[iClimb] == 0.5,
                            eta_prop[iCruise] == 0.6,
                            eta_prop[iLoiter] == 0.7
                            ])
        # Propulsive efficiency variation with different flight segments,
        # will change depending on propeller characteristics

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
        # Breguet Range
        z_bre = VectorVariable(NSeg, 'z_{bre}', '-', 'breguet coefficient')
        BSFC = Variable('BSFC', 0.7, 'lb/hr/hp',
                              'brake specific fuel consumption')
        t = VectorVariable(NSeg, 't', 'days', 'time per flight segment')
        t_cruise = Variable('t_{cruise}', 0.5, 'days', 'time to station')
        t_station = Variable('t_{station}', 5, 'days', 'time on station')
        R = Variable('R', 200, 'nautical_miles', 'range to station')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([z_bre >= V*t*BSFC*g*CD/CL/eta_prop,
                            R <= V[iCruise]*t[iCruise],
                            t[iLoiter] >= t_station/NLoiter,
                            t[iCruise[0]] <= t_cruise,
                            W_fuel/W_end >= te_exp_minus1(z_bre, 3)])

        #----------------------------------------------------
        # BSFC model
        
        #h_fuel = Variable('h_{fuel}', 42e6, 'J/kg', 'heat of combustion')
        #eta_th = VectorVariable(NSeg, r'\eta_{th}', np.linspace(0.25, 0.25, NSeg),  '-', 'thermal efficiency')
        #theta = VectorVariable(NSeg, r'\theta', np.linspace(1,1,NSeg), '-', 'throttle setting')
        #Q = VectorVariable(NSeg, 'Q',np.linspace(1.9,1.9,NSeg), 'N*m', 'engine torque')
        #Q_max = Variable('Q_{max}', 2.58, 'N*m', 'maximum torque')
        #eta_otto = Variable(r'\eta_{otto}', 0.5, '-', 'maximum thermal efficiency, otto')
        #A0 = Variable('A_0', 1.3548, '-', 'thermal coefficient A0')
        #A1 = Variable('A_1', 
        #
        #constraints.extend([BSFC >= 1/h_fuel/eta_th/theta*(Q/Q_max),
        #                    ])

        #----------------------------------------------------
        # Aerodynamics model

        CLmax = Variable('C_{L-max}', 1.5, '-', 'Maximum lift coefficient')
        e = Variable('e', 0.9, '-', 'Spanwise efficiency')
        AR = Variable('AR', '-', 'Aspect ratio')
        b = Variable('b', 'ft', 'Span')
        mu = Variable(r'\mu', 1.5e-5, 'N*s/m^2', 'Dynamic viscosity')
        Re = VectorVariable(NSeg, 'Re', '-', 'Reynolds number')
        Cf = VectorVariable(NSeg, 'C_f', '-', 'wing skin friction coefficient')
        Kwing = Variable('K_{wing}', 1.3, '-', 'wing form factor')
        cl_16 = Variable('cl_{16}', 0.0001, '-', 'profile stall coefficient')

        # fuselage drag 
        Kfuse = Variable('K_{fuse}', 1.1, '-', 'Fuselage form factor')
        S_fuse = Variable('S_{fuse}', 'ft^2', 'Fuselage surface area')
        Cffuse = Variable('C_{f-fuse}', '-', 'Fuselage skin friction coefficient')
        CDfuse = Variable('C_{D-fuse}', '-', 'fueslage drag')
        l_fuse = Variable('l_{fuse}',3,'ft', 'fuselage length')
        Refuse = Variable('Re_{fuse}', '-', 'fuselage Reynolds number')

        # landing gear
        A_rearland = Variable('A_{rear-land}', 6, 'in^2', 'rear landing gear frontal area')
        A_frontland = Variable('A_{front-land}', 6, 'in^2', 'front landing gear frontal area')
        CDland = Variable('C_{D-land}', 0.2, '-', 'drag coefficient landing gear')
        CDAland = Variable('CDA_{land}', '-', 'normalized drag coefficient landing gear')
        
        constraints.extend([CD >= CDfuse + CDAland + 2*Cf*Kwing + CL**2/(pi*e*AR)
                                + cl_16*CL**16,
                            CDAland >= (2*CDland*A_rearland + CDland*A_frontland)/S,
                            b**2 == S*AR,
                            CL <= CLmax, 
                            Re == rho*V/mu*(S/AR)**0.5,
                            Cf >= 0.074/Re**0.2,
                            CDfuse >= Kfuse*S_fuse*Cffuse/S,
                            Refuse == rho*V/mu*l_fuse,
                            Cffuse >= 0.074/Refuse**0.2,
                            ])

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

        constraints.extend([#h <= [20000, 20000, 20000]*units.m,  # Model valid to top of troposphere
                            T_sl >= T_atm + L_atm*h,     # Temp decreases w/ altitude
                            rho == p_sl*T_atm**(TH-1)/R_spec/(T_sl**TH)])
            # http://en.wikipedia.org/wiki/Density_of_air#Altitude

        #----------------------------------------------------
        # altitude constraints
        h_station = Variable('h_{station}', 15000, 'ft', 'minimum altitude at station')
        h_min = Variable('h_{min}', 5000, 'ft', 'minimum cruise altitude')

        constraints.extend([h[iLoiter] >= h_station,
                            h[iCruise] >= h_min,
                            h[iClimb] >= h_min, 
                            t[iClimb[0]]*h_dot == h_min, 
                            # still need to determine min cruise altitude, 
                            #and make these variables independent of user-input numbers
                            t[iClimb[1]]*h_dot == 10000*units('ft'),
                            ])

        #----------------------------------------------------
        # Weight breakdown

        W_cent = Variable('W_{cent}', 'lbf', 'Center aircraft weight')
        W_fuse = Variable('W_{fuse}', 'lbf', 'fuselage weight') 
        W_wing = Variable('W_{wing}', 'lbf', 'Total wing structural weight')
        m_fuse = Variable('m_{fuse}', 'kg', 'fuselage mass')
        m_cap = Variable('m_{cap}', 'kg', 'Cap mass')
        m_skin = Variable('m_{skin}','kg','Skin mass')
        m_tail = Variable('m_{tail}', 0.75, 'kg', 'tail mass')
        W_fueltot = Variable('W_{fueltot}', 'lbf', 'total fuel weight')

        W_fueltot = 0
        for j in range(0,NSeg):
            W_fueltot += W_fuel[j]

        constraints.extend([W_wing >= m_skin*g + m_cap*g,
                            W_fuse == m_fuse*g,
                            W_cent >= W_fueltot + W_pay + W_engtot + W_fuse + W_avionics,
                            W_zfw >= W_pay + W_engtot + W_fuse + W_wing + m_tail*g + W_avionics]) 

        #----------------------------------------------------
        # Structural model

        # Structural parameters
        rho_skin = Variable(r'\rho_{skin}', 0.1, 'g/cm^2', 'Wing Skin Density') 
        rho_cap = Variable(r'\rho_{cap}',1.76, 'g/cm^3', 'Density of CF cap')
        E_cap = Variable('E_{cap}', 2e7, 'psi', 'Youngs modulus of CF cap')
        sigma_cap = Variable(r'\sigma_{cap}', 475e6,'Pa', 'Cap stress') 
        
        # Structural lengths
        h_spar = Variable('h_{spar}', 'm', 'Spar height') 
        t_cap = Variable('t_{cap}', 0.028, 'in', 'Spar cap thickness') 
        #arbitrarily placed based on available cf
        w_cap = Variable('w_{cap}', 'in', 'Spar cap width')
        c = Variable('c', 'ft', 'Wing chord') #assumes straight, untapered wing

        # Structural ratios
        tau = Variable(r'\tau', 0.12,'-', 'Airfoil thickness ratio') #find better number
        LoverA = Variable('LoverA', 'lbf/ft^2', 'Wing loading')
        lambda_c = Variable(r'\lambda_c', '-', 'Taper ratio')

        # Structural areas
        A_capcent = Variable('A_{capcent}', 'm**2', 'Cap area at center')
        A_cap = Variable('A_{cap}', 'm**2', 'Cap area') #currently assumes constant area

        # Structural volumes
        Vol_cap = Variable('Vol_{cap}', 'm**3', 'Cap volume')

        # Structural evaluation parameters
        M_cent = Variable('M_cent', 'N*m', 'Center bending moment')
        F = Variable('F', 'N', 'Load on wings')
        SL = Variable('SL', 'Pa', 'Shear load') #need to add constraint
        N_Max = Variable('N_{Max}', 5,'-', 'Load factor') 
        #load rating for max number of g's
        P_cap = Variable('P_{cap}', 'N', 'Cap load')
        delta_tip = Variable(r'\delta_{tip}', 'ft', 'Tip deflection') 
        #need to add constraint

        constraints.extend([m_skin >= rho_skin*S*2,
                            F >= W_cent*N_Max,
                            c == S/b,
                            M_cent >= b*F/8,
                            P_cap >= M_cent/h_spar,
                            A_capcent >= P_cap/sigma_cap,
                            Vol_cap >= A_capcent*b/3,
                            m_cap == rho_cap*Vol_cap,
                            h_spar <= tau*c,
                            w_cap == A_capcent/t_cap,
                            LoverA == MTOW/S,
                            delta_tip == b**2*sigma_cap/(4*E_cap*h_spar),
                            delta_tip <= b/5]) 

        #----------------------------------------------------
        # Fuselage model

        # Constants
        rho_fuel = Variable(r'\rho_{fuel}', 6.01, 'lbf/gallon', 'density of 100LL')

        # Non-dimensional variables
        k1fuse = Variable('k_{1-fuse}', 2.5, '-', 'fuselage form factor 1')
        k2fuse = Variable('k-{2-fuse}', 20, '-', 'fuselage form factor 2')

        # Volumes 
        Vol_fuel = Variable('Vol_{fuel}',  'm**3', 'Fuel Volume')
        Vol_fuse = Variable('Vol_{fuse}', 'm**3', 'fuselage volume')



        constraints.extend([m_fuse >= S_fuse*rho_skin,
                            (l_fuse/k1fuse)**3 == Vol_fuse,
                            (S_fuse/k2fuse)**3 == Vol_fuse**2,
                            Vol_fuel >= W_fueltot/rho_fuel,
                            Vol_fuse >= Vol_fuel+Vol_avionics+Vol_pay])

        #----------------------------------------------------
        # wind speed model

        V_wind = VectorVariable(NLoiter,'V_{wind}',[30,20,35,15,25], 'm/s', 'wind speed')

        constraints.extend([V[iLoiter] >= V_wind])

        objective = MTOW 
        return objective, constraints

if __name__ == '__main__':
    M = GasPoweredHALE()
    sol = M.solve()

    #----------------------------------------------
    # post processing
    
    #M.substitutions.update({'BSFC': ('sweep', np.linspace(0.3,1,15))})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
    #
    #BSFC = sol('BSFC')
    #MTOW = sol('MTOW')
    #b = sol('b')

    #plt.close()
    #plt.plot(BSFC,MTOW)
    #plt.xlabel('BSFC [lbf/hp/hr]')
    #plt.ylabel('MTOW [lbf]')
    #plt.savefig('BSFCvsMTOW.png')

    #plt.close()
    #plt.plot(BSFC,b)
    #plt.xlabel('BSFC [lbf/hp/hr]')
    #plt.ylabel('b [ft]')
    #plt.savefig('BSFCvsb.png')

    #M.substitutions.update({r'\eta_{prop}':('sweep', np.linspace(0.5,1,15)), 'BSFC':0.7})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

    #eta_prop = sol(r'\eta_{prop}')
    #MTOW = sol('MTOW')
    #b = sol('b')

    #plt.close()
    #plt.plot(eta_prop,MTOW)
    #plt.xlabel('eta_prop ')
    #plt.ylabel('MTOW [lbf]')
    #plt.savefig('eta_propvsMTOW.png')

    #plt.close()
    #plt.plot(eta_prop,b)
    #plt.xlabel('eta_prop ')
    #plt.ylabel('b [ft]')
    #plt.savefig('eta_propvsb.png')
    #plt.close()
        


