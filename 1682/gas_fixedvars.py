from numpy import pi
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

class GasPoweredHALE(Model):
    def __init__(self, **kwargs):

        # define number of segments
        NSeg = 9 # number of flight segments
        NCruise = 2 # number of cruise segments
        NClimb = 2 # number of climb segments
        NLoiter = NSeg - NCruise - NClimb # number of loiter segments
        iCruise = [1, -1] # cuise index
        iLoiter = [3] # loiter index
        for i in range(4, NSeg-1): iLoiter.append(i)
        iClimb = [0, 2] # climb index

        constraints = []

        #----------------------------------------------------
        # Fuel weight model 

        MTOW = Variable('MTOW', 'lbf', 'max take off weight')
        W_end = VectorVariable(NSeg, 'W_{end}', 'lbf', 'segment-end weight')
        W_fuel = VectorVariable(NSeg, 'W_{fuel}', 'lbf',
                                'segment-fuel weight')
        W_zfw = Variable('W_{zfw}', 'lbf', 'Zero fuel weight')
        W_begin = W_end.left # define beginning of segment weight
        W_begin[0] = MTOW 

        # Payload model
        W_pay = Variable('W_{pay}', 10, 'lbf', 'Payload weight')
        Vol_pay = Variable('Vol_{pay}', 1, 'ft^3', 'Payload volume')

        # Avionics model
        W_avionics = Variable('W_{avionics}', 8, 'lbf', 'Avionics weight')
        Vol_avionics = Variable('Vol_{avionics}', 0.125, 'ft^3', 'Avionics volume')

        # end of first segment weight + first segment fuel weight must be greater 
        # than MTOW.  Each end of segment weight must be greater than the next end
        # of segment weight + the next segment fuel weight. The last end segment
        # weight must be greater than the zero fuel weight

        constraints.extend([MTOW >= W_end[0] + W_fuel[0], 
                            W_end[:-1] >= W_end[1:] + W_fuel[1:], 
                            W_end[-1] >= W_zfw])
        
        #----------------------------------------------------
        # Steady level flight model
        t = VectorVariable(NSeg, 't', 'days', 'time per flight segment')
        CD = VectorVariable(NSeg, 'C_D', '-', 'Drag coefficient')
    	CL = VectorVariable(NSeg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(NSeg, 'V', 'm/s', 'cruise speed')
        rho = VectorVariable(NSeg, r'\rho', 'kg/m^3', 'air density')
        S = Variable('S', 'ft^2', 'wing area')
        eta_prop = VectorVariable(NSeg, r'\eta_{prop}', '-',
                                  'propulsive efficiency')
        eta_propCruise = Variable(r'\eta_{prop-cruise}',0.6,'-','propulsive efficiency in cruise')
        eta_propClimb = Variable(r'\eta_{prop-climb}',0.5,'-','propulsive efficiency in climb')
        eta_propLoiter = Variable(r'\eta_{prop-loiter}',0.7,'-','propulsive efficiency in loiter')
        P_shaft = VectorVariable(NSeg, 'P_{shaft}', 'hp', 'Shaft power')
        T = VectorVariable(NSeg, 'T', 'lbf', 'Thrust')

        # Climb model
        h_dot = VectorVariable(NClimb, 'h_{dot}', 'ft/min', 'Climb rate')
        # Currently climb rate is being optimized to reduce fuel consumption. 
        # In future, could implement min climb rate. 
        
        constraints.extend([P_shaft == T*V/eta_prop, 
                            T >= 0.5*rho*V**2*CD*S, 
                            T[iClimb] >= 0.5*rho[iClimb]*V[iClimb]**2*CD[iClimb]*S + 
                                         W_begin[iClimb]*h_dot/V[iClimb], 
                            0.5*rho*CL*S*V**2 == (W_end*W_begin)**0.5,
                            eta_prop[iLoiter] == eta_propLoiter,
                            eta_prop[iCruise] == eta_propCruise,
                            eta_prop[iClimb] == eta_propClimb
                            ])
        # Propulsive efficiency variation with different flight segments, 
        # will change depending on propeller characteristics

        #----------------------------------------------------
        # altitude constraints
        h_station = Variable('h_{station}', 15000, 'ft', 'minimum altitude at station')
        h_cruise = Variable('h_{cruise}', 5000, 'ft', 'minimum cruise altitude')
        h = VectorVariable(NSeg, 'h', 
                [5000,5000,15000,15000,15000,
                 15000,15000,15000,5000], 'ft', 'altitude')
        deltah = Variable(r'\delta_h', h_station.value-h_cruise.value, 'ft', 'delta height') 
        constraints.extend([t[iClimb[0]]*h_dot[0] >= h_cruise, 
                            t[iClimb[1]]*h_dot[1] >= deltah
                            ])

        #----------------------------------------------------
        # Engine Model (DF35)

        W_engtot = Variable('W_{eng-tot}', 8, 'lbf', 'Installed engine weight')
                #conservative for 4.2 engine complete with prop, generator and structures
        FuelOilFrac = Variable('FuelOilFrac',.98,'-','Fuel-oil fraction')
        BSFC_min = Variable('BSFC_{min}', 0.32, 'kg/kW/hr', 'Minimum BSFC')
        BSFC = VectorVariable(NSeg, 'BSFC', 'lb/hr/hp', 
                              'brake specific fuel consumption') 
        RPM_max = Variable('RPM_{max}', 9000, 'rpm', 'Maximum RPM')
        RPM = VectorVariable(NSeg, 'RPM', 'rpm', 'Engine operating RPM')
        P_shaftmaxMSL = Variable('P_{shaft-maxMSL}', 5.84, 'hp', 
                                 'Max shaft power at MSL')
        P_shaftmax = VectorVariable(NSeg, 'P_{shaft-max}',
                [P_shaftmaxMSL.value*(1-(0.906**(1/0.15)*(v.value/h_station.value)**0.92)) for v in h])
        P_shaftmaxMSL = Variable('P_{shaft-maxMSL}', 5.84, 'hp', 
                                 'Max shaft power at MSL')
        P_avn = VectorVariable(NSeg, 'P_{avn}', [40,40,40,50,50,50,50,50,40], 'watts', 'avionics power')
        P_shafttot = VectorVariable(NSeg, 'P_{shaft-tot}', 'hp', 'total power need including power draw from avionics')

        V_max = VectorVariable(NSeg,'V_{max}','m/s','Maximum velocity')
        
        # Engine Power Relations
        constraints.extend([P_shaftmax >= P_shafttot, 
                            P_shafttot >= P_shaft + P_avn*1.25, #need to figure out better value for alternator efficiency
                            (BSFC/BSFC_min)**0.129 >= 2*.486*(RPM/RPM_max)**-0.141 + \
                                                      0.0268*(RPM/RPM_max)**9.62, 
                            (P_shafttot/P_shaftmax)**0.1 == 0.999*(RPM/RPM_max)**0.292, 
                            RPM <= RPM_max, 
                            ])

        #----------------------------------------------------
        # Breguet Range
        z_bre = VectorVariable(NSeg, 'z_{bre}', '-', 'breguet coefficient')
        t_cruise = Variable('t_{cruise}', 1, 'days', 'time to station')
        t_station = Variable('t_{station}', 6, 'days', 'time on station')
        R = Variable('R', 200, 'nautical_miles', 'range to station')
        R_cruise = Variable('R_{cruise}',180,'nautical_miles','range to station during climb')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')

        constraints.extend([z_bre >= P_shafttot*t*BSFC*g/W_end,
                            R_cruise <= V[iCruise[0]]*t[iCruise[0]], 
                            R <= V[iCruise[1]]*t[iCruise[1]], 
                            t[iLoiter] >= t_station/NLoiter, 
                            sum(t[[0,1,2]]) <= t_cruise, 
                            FuelOilFrac*W_fuel/W_end >= te_exp_minus1(z_bre, 3)])
        
        #----------------------------------------------------
        # Atmosphere model
        gamma = Variable(r'\gamma', 1.4, '-', 'Heat capacity ratio of air')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        rho_sl = Variable(r'\rho_{sl}',1.225, 'kg/m^3', 'density at sea level')
        T_sl = Variable('T_{sl}', 288.15, 'K', 'Temperature at sea level')
        mu_sl = Variable(r'\mu_{sl}', 1.789*10**-5, 'N*s/m^2','Dynamic viscosity at sea level')
        L_atm = Variable('L_{atm}', 0.0065, 'K/m', 'Temperature lapse rate')
        T_atm = VectorVariable(NSeg, 'T_{atm}', 'K', 'Air temperature')
        a_atm = VectorVariable(NSeg, 'a_{atm}', 'm/s', 'Speed of sound at altitude')
        mu_atm = VectorVariable(NSeg,r'\mu', 'N*s/m^2', 'Dynamic viscosity')

        R_spec = Variable('R_{spec}', 287.058, 'J/kg/K', 'Specific gas constant of air')

        # Atmospheric variation with altitude (valid from 0-7km of altitude)
        constraints.extend([(rho/rho_sl)**0.1 == 0.954*(h/h_station)**(-0.0284),
                            (T_atm/T_sl)**0.1 == 0.989*(h/h_station)**(-0.00666),
                            (mu_atm/mu_sl)**0.1 == 0.991*(h/h_station)**(-0.00529)
                            ])


        #----------------------------------------------------
        # Aerodynamics model

        CLmax = Variable('C_{L-max}', 1.5, '-', 'Maximum lift coefficient')
        e = Variable('e', 0.9, '-', 'Spanwise efficiency')
        AR = Variable('AR', '-', 'Aspect ratio')
        b = Variable('b', 'ft', 'Span')
        Re = VectorVariable(NSeg, 'Re', '-', 'Reynolds number')

        # fuselage drag 
        Kfuse = Variable('K_{fuse}', 1.1, '-', 'Fuselage form factor')
        S_fuse = Variable('S_{fuse}', 'ft^2', 'Fuselage surface area')
        Cffuse = VectorVariable(NSeg, 'C_{f-fuse}', '-', 'Fuselage skin friction coefficient')
        CDfuse = Variable('C_{D-fuse}', '-', 'fueslage drag')
        l_fuse = Variable('l_{fuse}', 'ft', 'fuselage length')
        l_cent = Variable('l_{cent}', 'ft', 'center fuselage length')
        Refuse = VectorVariable(NSeg, 'Re_{fuse}', '-', 'fuselage Reynolds number')
        Re_ref = Variable("Re_{ref}", 3e5, "-", "Reference Re for cdp")
        cdp = VectorVariable(NSeg, "c_{dp}", "-", "wing profile drag coeff")

        constraints.extend([
            CD >= (CDfuse + cdp + CL**2/(pi*e*AR))*1.3,
            cdp >= ((0.006 + 0.005*CL**2 + 0.00012*CL**10)*(Re/Re_ref)**-0.3),
            b**2 == S*AR, 
            CL <= CLmax, 
            Re == rho*V/mu_atm*(S/AR)**0.5, 
            CDfuse >= Kfuse*S_fuse*Cffuse/S, 
            Refuse == rho*V/mu_atm*l_fuse, 
            Cffuse >= 0.455/Refuse**0.3, 
            ])

        #----------------------------------------------------
        # Weight breakdown

        W_cent = Variable('W_{cent}', 'lbf', 'Center aircraft weight')
        W_fuse = Variable('W_{fuse}', 'lbf', 'fuselage weight') 
        W_wing = Variable('W_{wing}', 'lbf', 'Total wing structural weight')
        W_fueltot = Variable('W_{fuel-tot}', 'lbf', 'total fuel weight')
        W_skid = Variable('W_{skid}', 3, 'lbf', 'skid weight')
        m_fuse = Variable('m_{fuse}', 'kg', 'fuselage mass')
        m_cap = Variable('m_{cap}', 'kg', 'Cap mass')
        m_skin = Variable('m_{skin}', 'kg', 'Skin mass')
        m_tail = Variable('m_{tail}', 1*1.5, 'kg', 'tail mass')
        m_rib = Variable('m_{rib}', 1.2*1.5, 'kg','rib mass')

        constraints.extend([W_wing >= m_skin*g + 1.2*m_cap*g, 
                            W_fuse >= m_fuse*g + m_rib*g, 
                            W_fueltot >= W_fuel.sum(),
                            W_cent >= W_fueltot + W_pay + W_engtot + W_fuse + W_avionics + W_skid, 
                            W_zfw >= W_pay + W_engtot + W_fuse + W_wing + m_tail*g +
                                     W_avionics + W_skid]) 

        #----------------------------------------------------
        # Structural model

        # Structural parameters
        rho_skin = Variable(r'\rho_{skin}', 0.1, 'g/cm^2', 'Wing Skin Density') 
        rho_cap = Variable(r'\rho_{cap}', 1.76, 'g/cm^3', 'Density of CF cap')
        E_cap = Variable('E_{cap}', 2e7, 'psi', 'Youngs modulus of CF cap')
        sigma_cap = Variable(r'\sigma_{cap}', 475e6, 'Pa', 'Cap stress') 
        
        # Structural lengths
        h_spar = Variable('h_{spar}', 'm', 'Spar height') 
        t_cap = Variable('t_{cap}', 0.028, 'in', 'Spar cap thickness') 
        #arbitrarily placed based on available cf
        w_cap = Variable('w_{cap}', 'in', 'Spar cap width')
        c = Variable('c', 'ft', 'Wing chord') #assumes straight, untapered wing

        # Structural ratios
        tau = Variable(r'\tau', 0.12, '-', 'Airfoil thickness ratio') #find better number
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
        N_Max = Variable('N_{Max}', 5, '-', 'Load factor') 
        #load rating for max number of g's
        P_cap = Variable('P_{cap}', 'N', 'Cap load')
        delta_tip = Variable(r'\delta_{tip}', 'ft', 'Tip deflection') 
        delta_tip_max = Variable(r'\delta_{tip-max}', 0.2, '-',
                                 'max tip deflection ratio') 

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
                            delta_tip/b <= delta_tip_max]) 

        #----------------------------------------------------
        # Fuselage model

        # Constants
        rho_fuel = Variable(r'\rho_{fuel}', 6.01, 'lbf/gallon', 'density of 100LL')

        # Non-dimensional variables
        k1fuse = Variable('k_{1-fuse}', 2.858, '-', 'fuselage form factor 1')
        k2fuse = Variable('k-{2-fuse}', 5.938, '-', 'fuselage form factor 2')
        w_cent = Variable('w_{cent}', 'ft', 'center fuselage width')
        fr = Variable('fr', 3.5, '-', 'fineness ratio fuselage')

        # Volumes
        Vol_fuel = Variable('Vol_{fuel}', 'm**3', 'Fuel Volume')
        Vol_fuse = Variable('Vol_{fuse}', 'm**3', 'fuselage volume')

        constraints.extend([m_fuse >= S_fuse*rho_skin, 
                            l_cent == fr*w_cent,
                            l_fuse >= l_cent*1.1,
                            (l_fuse/k1fuse)**3 == Vol_fuse, 
                            (S_fuse/k2fuse)**3 == Vol_fuse**2, 
                            Vol_fuse >= l_cent*w_cent**2,
                            Vol_fuel >= W_fuel.sum()/rho_fuel, 
                            l_cent*w_cent**2 >= Vol_fuel+Vol_avionics+Vol_pay])

        #----------------------------------------------------
        # wind speed model

        V_wind = Variable('V_{wind}', 25, 'm/s', 'wind speed')
        #wd_cnst = Variable('wd_{cnst}', 0.001077, 'm/s/ft', 
        #                   'wind speed constant predicted by model')
        #                    #0.002 is worst case, 0.0015 is mean at 45d
        #wd_ln = Variable('wd_{ln}', 8.845, 'm/s',
        #                 'linear wind speed variable')
                        #13.009 is worst case, 8.845 is mean at 45deg
        #h_cruise = Variable('h_{cruise}', 11800, 'ft', 'minimum height')
        #h_max = Variable('h_{max}', 20866, 'ft', 'maximum height')

        constraints.extend([#V_wind >= wd_cnst*h + wd_ln, 
                            V[iLoiter] >= V_wind,
                            #h[NCruise] >= 11800*units('ft')
                            ])

        objective = MTOW
        Model.__init__(self,objective,constraints, **kwargs)

if __name__ == '__main__':
    M = GasPoweredHALE()
    #M.substitutions.update({M["t_{station}"]:6})
    sol = M.solve('cvxopt') 
    P_shaftmaxMSL = 4355; #W
    P_shaftmax = sol('P_{shaft-max}')
    V = sol('V')
    CD = sol('C_D')
    rho = sol(r'\rho')
    S = sol('S')*.3048**2 #m^2
    eta_prop = sol(r'\eta_{prop}')
    T = sol('T')
    rhoSL = sol(r'\rho_{sl}')
    h = sol('h_{station}')*0.3048 #m
    wd_cnst = 0.001077/.3048
    wd_ln = 8.845
    eta_propLoiter = 0.7

    #Finding maximum altitude (post-processing)

    h_max = np.ones([1,5])*15000*0.3048;
    rho_maxalt= np.zeros([1,5]);
    V_maxalt = np.ones([1,5])*25;
    P_shaftmaxalt = np.zeros([1,5])

    for i in range(0,h_max.shape[1]):
        while P_shaftmaxalt[0,i] == 0 or P_shaftmaxalt[0,i]*eta_propLoiter >= .5*rho_maxalt[0,i]*V_maxalt[0,i]**3.*CD[3+i]*S:
            h_max[0,i] = h_max[0,i] + 250.
            V_maxalt[0,i] = wd_cnst*h_max[0,i] + wd_ln
            rho_maxalt[0,i] = (0.954*(h_max[0,i]/h)**(-0.0284))**10*rhoSL
            P_shaftmaxalt[0,i] = P_shaftmaxMSL*(1-(0.906**(1/0.15)*(h_max[0,i]/h)**0.92))
            
    h_max = h_max/0.3048 # Back to feet
    print h_max

    #V_maxTO = ((P_shaftmaxMSL*0.5/0.5/rhoSL/S/CD[0]).to('m^3/s^3'))**(1./3.)
    #V_max = ((P_shaftmax*0.7/0.5/rho/S/CD).to('m^3/s^3'))**(1./3.)

    #print "Max velocity at TO: {:.1f~}".format(V_maxTO)
    #print "Max velocity at loiter: {:.1f~}".format(V_max[3])

    #----------------------------------------------
    # post processing
    
    #M.substitutions.update({'MTOW': ('sweep', np.linspace(70, 150, 15))})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
    #
    #MTOW = sol('MTOW')
    #t_station = sol('t_{station}')

    #plt.close()
    #plt.plot(MTOW, t_station)
    #plt.xlabel('MTOW [lbf]')
    #plt.ylabel('t_station [days]')
    #plt.grid()
    #plt.savefig('tvsMTOW.png')

    #M.substitutions.update({'R':('sweep', np.linspace(100, 600, 15)), 'MTOW':143})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

    #R = sol('R')
    #t_station = sol('t_{station}')

    #plt.close()
    #plt.plot(R, t_station)
    #plt.xlabel('R [nm]')
    #plt.grid()
    #plt.ylabel('t_station [days]')
    #plt.axis([0,600,0,10])
    #plt.savefig('tvsR.png')

    #M.substitutions.update({'h_{cruise}': ('sweep', np.linspace(1000,15000,15))})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

    #h_cruise = sol('h_{cruise}')
    #t_station = sol('t_{station}')
    #deltah = sol(r'\delta_h')

    #plt.close()
    #plt.plot(h_cruise, t_station)
    #plt.xlabel('h_cruise [ft]')
    #plt.grid()
    #plt.ylabel('t_station [days]')
    #plt.savefig('tvsh_cruise.png')

    #M.substitutions.update({'h_{station}':('sweep', np.linspace(15000,20000, 15)), r'\rho' 'h_{cruise}':5000})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

    #h_station = sol('h_{station}')
    #t_station = sol('t_{station}')

    #plt.close()
    #plt.plot(h_station, t_station)
    #plt.xlabel('h_station [lbf]')
    #plt.grid()
    #plt.ylabel('t_station [days]')
    #plt.savefig('tvsh_station.png')

    #M.substitutions.update({'V_{wind}':('sweep', np.linspace(1, 40, 40))})
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

    #V_wind = sol('V_{wind}')
    #t_station = sol('t_{station}')

    #plt.close()
    #plt.plot(V_wind, t_station)
    #plt.xlabel('V_wind [m/s]')
    #plt.ylabel('t_station [days]')
    #plt.axis([0,40,0,10])
    #plt.grid()
    #plt.savefig('tvsV_wind.png')

    #M.substitutions.update({'V_{wind}':('sweep', np.linspace(1, 40, 40)), 'MTOW': })
    #sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

    #V_wind89 = sol('V_{wind}')
    #t_station89 = sol('t_{station}')

    #plt.plot(V_wind83, t_station83, V_wind89, t_station89)
    #plt.xlabel('V_wind [m/s]')
    #plt.ylabel('t_station [days]')
    #plt.grid()
    #plt.savefig('tvsV_wind89.png')
