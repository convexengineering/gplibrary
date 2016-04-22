from numpy import pi
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
import gpkit
import numpy as np
gpkit.settings['latex_modelname'] = False

PLOT = False

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
        # altitude constraints
        h_station = Variable('h_{station}', 15000, 'ft',
                             'minimum altitude at station')
        h_cruise = Variable('h_{cruise}', 5000, 'ft',
                            'minimum cruise altitude')
        h = np.array([h_cruise]*2 + [h_station]*6 + [h_cruise])
        deltah = Variable('\\delta_h', h_station.value-h_cruise.value, 'ft',
                          'delta height')
        t = VectorVariable(NSeg, 't', 'days', 'time per flight segment')
        h_dot = VectorVariable(NClimb, 'h_{dot}', 'ft/min', 'Climb rate')
        h_dotmin = Variable('h_{dot-min}', 100, 'ft/min',
                            'minimum necessary climb rate')

        constraints.extend([t[iClimb[0]]*h_dot[0] >= h_cruise,
                            t[iClimb[1]]*h_dot[1] >= deltah,
                            h_dot >= h_dotmin
                           ])

        #----------------------------------------------------
        # Atmosphere model
        gamma = Variable('\\gamma', 1.4, '-', 'Heat capacity ratio of air')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        L_atm = Variable('L_{atm}', 0.0065, 'K/m', 'Temperature lapse rate')
        T_sl = Variable('T_{sl}', 288.15, 'K', 'Temperature at sea level')
        T_atm = VectorVariable(NSeg, 'T_{atm}',
                lambda T_sl, L_atm, h: T_sl - L_atm*h, 'K',
                'Air temperature', args=[T_sl, L_atm, h])
        a_atm = VectorVariable(NSeg, 'a_{atm}', 'm/s',
                               'Speed of sound at altitude')
        mu_atm = VectorVariable(NSeg,'\\mu', 'N*s/m^2', 'Dynamic viscosity')
        mu_sl = Variable('\\mu_{sl}', 1.789*10**-5, 'N*s/m^2',
                         'Dynamic viscosity at sea level')
        R_spec = Variable('R_{spec}', 287.058, 'J/kg/K',
                          'Specific gas constant of air')
        h_ref = Variable('h_{ref}', 15000, 'ft', 'ref altitude')
        rho = VectorVariable(NSeg, '\\rho', 'kg/m^3', 'air density')

        # Atmospheric variation with altitude (valid from 0-7km of altitude)
        constraints.extend([
            rho == p_sl*T_atm**(5.257-1)/R_spec/(T_sl**5.257),
            (mu_atm/mu_sl)**0.1 == 0.991*(h/h_station)**(-0.00529)
            ])

        #----------------------------------------------------
        # Fuel weight model

        MTOW = Variable('MTOW', 143, 'lbf', 'max take off weight')
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
        Vol_avionics = Variable('Vol_{avionics}', 0.125, 'ft^3',
                                'Avionics volume')

        # end of first segment weight + first segment fuel weight must be
        # greater  than MTOW.  Each end of segment weight must be greater
        # than the next end of segment weight + the next segment fuel weight.
        # The last end segment weight must be greater than the zero fuel
        # weight

        constraints.extend([MTOW >= W_end[0] + W_fuel[0],
                            W_end[:-1] >= W_end[1:] + W_fuel[1:],
                            W_end[-1] >= W_zfw])

        #----------------------------------------------------
        # Steady level flight model

        CD = VectorVariable(NSeg, 'C_D', '-', 'Drag coefficient')
        CL = VectorVariable(NSeg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(NSeg, 'V', 'm/s', 'cruise speed')
        S = Variable('S', 'ft^2', 'wing area')
        eta_prop = VectorVariable(NSeg, '\\eta_{prop}', '-',
                                  'propulsive efficiency')
        eta_propCruise = Variable('\\eta_{prop-cruise}', 0.6, '-',
                                  'propulsive efficiency in cruise')
        eta_propClimb = Variable('\\eta_{prop-climb}', 0.5, '-',
                                 'propulsive efficiency in climb')
        eta_propLoiter = Variable('\\eta_{prop-loiter}', 0.7, '-',
                                  'propulsive efficiency in loiter')
        P_shaft = VectorVariable(NSeg, 'P_{shaft}', 'hp', 'Shaft power')
        T = VectorVariable(NSeg, 'T', 'lbf', 'Thrust')

        # Climb model
        # Currently climb rate is being optimized to reduce fuel consumption.
        # In future, could implement min climb rate.

        constraints.extend([
            P_shaft == T*V/eta_prop,
            T >= 0.5*rho*V**2*CD*S,
            T[iClimb] >= (0.5*rho[iClimb]*V[iClimb]**2*CD[iClimb]*S +
                          W_begin[iClimb]*h_dot/V[iClimb]),
            0.5*rho*CL*S*V**2 == (W_end*W_begin)**0.5,
            eta_prop[iLoiter] == eta_propLoiter,
            eta_prop[iCruise] == eta_propCruise,
            eta_prop[iClimb] == eta_propClimb
            ])
        # Propulsive efficiency variation with different flight segments,
        # will change depending on propeller characteristics

        #----------------------------------------------------
        # Engine Model (DF35)

        W_engtot = Variable('W_{eng-tot}', 7.1, 'lbf',
                            'Installed engine weight')
        W_engref = Variable('W_{eng-ref}', 4.4107, 'lbf',
                            'Reference engine weight')
        FuelOilFrac = Variable('FuelOilFrac', .98, '-', 'Fuel-oil fraction')
        P_shaftref = Variable('P_{shaft-ref}', 2.295, 'hp',
                              'reference shaft power')
        BSFC_min = Variable('BSFC_{min}', 0.32, 'kg/kW/hr', 'Minimum BSFC')
        BSFC = VectorVariable(NSeg, 'BSFC', 'lb/hr/hp',
                              'brake specific fuel consumption')
        RPM_max = Variable('RPM_{max}', 9000, 'rpm', 'Maximum RPM')
        RPM = VectorVariable(NSeg, 'RPM', 'rpm', 'Engine operating RPM')
        P_shaftmaxMSL = Variable('P_{shaft-maxMSL}', 5.84, 'hp',
                                 'Max shaft power at MSL')
        P_shaftmax = VectorVariable(NSeg, 'P_{shaft-max}',
                                    lambda P_shaftmaxMSL, h_ref, h:
                                    P_shaftmaxMSL*(1-0.5178*(h/h_ref)**0.92),
                                    'hp', 'Max shaft power at altitude',
                                    args=[P_shaftmaxMSL, h_ref, h])
        #Lfactor = VectorVariable(NSeg, 'L_factor', '-',
        #                         'Max shaft power loss factor')
        P_avn = Variable('P_{avn}', 40, 'watts', 'avionics power')
        P_pay = Variable('P_{pay}', 10, 'watts', 'payload power')
        P_shafttot = VectorVariable(NSeg, 'P_{shaft-tot}', 'hp',
                'total power need including power draw from avionics')
        m_dotfuel = VectorVariable(NSeg, 'm_{dot-fuel}', 'lb/sec',
                                   'fuel flow rate')
        eta_alternator = Variable('\\eta_{alternator}', 0.8, '-',
                                  'alternator efficiency')

        # Engine Weight Constraints
        constraints.extend([
            #Lfactor == 0.906**(1/0.15)*(h/h_station)**0.92,
            #P_shaftmax/P_shaftmaxMSL + Lfactor <= 1,
            P_shaftmax >= P_shafttot,
            P_shafttot[iCruise] >= P_shaft[iCruise] + P_avn/eta_alternator,
            P_shafttot[iClimb] >= P_shaft[iClimb] + P_avn/eta_alternator,
            P_shafttot[iLoiter] >= (P_shaft[iLoiter] +
                                    (P_avn+P_pay)/eta_alternator),
            (BSFC/BSFC_min)**0.129 >= (2*.486*(RPM/RPM_max)**-0.141 +
                                       0.0268*(RPM/RPM_max)**9.62),
            (P_shafttot/P_shaftmax)**0.1 == 0.999*(RPM/RPM_max)**0.292,
            RPM <= RPM_max,
            ])

        #----------------------------------------------------
        # Breguet Range
        z_bre = VectorVariable(NSeg, 'z_{bre}', '-', 'breguet coefficient')
        t_cruise = Variable('t_{cruise}', 1, 'days', 'time to station')
        t_station = Variable('t_{station}', 'days', 'time on station')
        R = Variable('R', 200, 'nautical_miles', 'range to station')
        R_cruise = Variable('R_{cruise}', 180, 'nautical_miles',
                            'range to station during climb')

        constraints.extend([
            z_bre >= P_shafttot*t*BSFC*g/W_end,
            R_cruise <= V[iCruise[0]]*t[iCruise[0]],
            R <= V[iCruise[1]]*t[iCruise[1]],
            t[iLoiter] >= t_station/NLoiter,
            sum(t[[0,1,2]]) <= t_cruise,
            FuelOilFrac*W_fuel/W_end >= te_exp_minus1(z_bre, 3)
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
        Cffuse = VectorVariable(NSeg, 'C_{f-fuse}', '-',
                                'Fuselage skin friction coefficient')
        CDfuse = Variable('C_{D-fuse}', '-', 'fueslage drag')
        l_fuse = Variable('l_{fuse}', 'ft', 'fuselage length')
        l_cent = Variable('l_{cent}', 'ft', 'center fuselage length')
        Refuse = VectorVariable(NSeg, 'Re_{fuse}', '-',
                                'fuselage Reynolds number')
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
        # landing gear
        #A_rearland = Variable('A_{rear-land}', 6, 'in^2',
        #                      'rear landing gear frontal area')
        #A_frontland = Variable('A_{front-land}', 6, 'in^2',
        #                       'front landing gear frontal area')
        #CDland = Variable('C_{D-land}', 0.2, '-',
        #                  'drag coefficient landing gear')
        #CDAland = Variable('CDA_{land}', '-',
        #                   'normalized drag coefficient landing gear')

        #constraints.extend([
        #    CD >= (CDfuse + 2*Cf*Kwing + CL**2/(pi*e*AR) + cl_16*CL**16 +
        #           CDAland),
        #    CDAland >= (2*CDland*A_rearland + CDland*A_frontland)/S
        #    ])

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
        m_tail = Variable('m_{tail}', 1.587, 'kg', 'tail mass')
        m_rib = Variable('m_{rib}', 1.36, 'kg','rib mass')
        W_fueltank = Variable('W_{fuel-tank}', 4, 'lbf', 'fuel tank weight')

        constraints.extend([
            W_wing >= m_skin*g + 1.2*m_cap*g,
            W_fuse >= m_fuse*g + m_rib*g,
            W_fueltot >= W_fuel.sum(),
            W_cent >= (W_fueltot + W_pay + W_engtot + W_fuse + W_avionics +
                       W_skid + W_fueltank),
            W_zfw >= (W_pay + W_engtot + W_fuse + W_wing + m_tail*g +
                      W_avionics + W_skid + W_fueltank)
            ])

        #----------------------------------------------------
        # Structural model

        # Structural parameters
        rho_skin = Variable('\\rho_{skin}', 0.1, 'g/cm^2',
                            'Wing Skin Density')
        rho_cap = Variable('\\rho_{cap}', 1.76, 'g/cm^3', 'Density of CF cap')
        E_cap = Variable('E_{cap}', 2e7, 'psi', 'Youngs modulus of CF cap')
        sigma_cap = Variable('\\sigma_{cap}', 475e6, 'Pa', 'Cap stress')

        # Structural lengths
        h_spar = Variable('h_{spar}', 'm', 'Spar height')
        t_cap = Variable('t_{cap}', 0.028, 'in', 'Spar cap thickness')
        #arbitrarily placed based on available cf
        w_cap = Variable('w_{cap}', 'in', 'Spar cap width')
        c = Variable('c', 'ft', 'Wing chord')
        #assumes straight, untapered wing

        # Structural ratios
        tau = Variable('\\tau', 0.12, '-', 'Airfoil thickness ratio')
        #find better number
        LoverA = Variable('LoverA', 'lbf/ft^2', 'Wing loading')
        lambda_c = Variable('\\lambda_c', '-', 'Taper ratio')

        # Structural areas
        A_capcent = Variable('A_{capcent}', 'm**2', 'Cap area at center')
        A_cap = Variable('A_{cap}', 'm**2', 'Cap area')
        #currently assumes constant area

        # Structural volumes
        Vol_cap = Variable('Vol_{cap}', 'm**3', 'Cap volume')

        # Structural evaluation parameters
        M_cent = Variable('M_cent', 'N*m', 'Center bending moment')
        F = Variable('F', 'N', 'Load on wings')
        SL = Variable('SL', 'Pa', 'Shear load') #need to add constraint
        N_Max = Variable('N_{Max}', 5, '-', 'Load factor')
        #load rating for max number of g's
        P_cap = Variable('P_{cap}', 'N', 'Cap load')
        delta_tip = Variable('\\delta_{tip}', 'ft', 'Tip deflection')
        delta_tip_max = Variable('\\delta_{tip-max}', 0.2, '-',
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
        rho_fuel = Variable('\\rho_{fuel}', 6.01, 'lbf/gallon',
                            'density of 100LL')

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
                            l_cent*w_cent**2 >= Vol_fuel+Vol_avionics+Vol_pay
                           ])

        #----------------------------------------------------
        # wind speed model

        V_wind = VectorVariable(2, 'V_{wind}', 'm/s', 'wind speed')
        wd_cnst = Variable('wd_{cnst}', 0.001077, 'm/s/ft',
                           'wind speed constant predicted by model')
                            #0.002 is worst case, 0.0015 is mean at 45d
        wd_ln = Variable('wd_{ln}', 8.845, 'm/s',
                         'linear wind speed variable')
                       #13.009 is worst case, 8.845 is mean at 45deg
        h_max = Variable('h_{max}', 20866, 'ft', 'maximum height')

        constraints.extend([V_wind[0] >= wd_cnst*h_station + wd_ln,
                            V_wind[1] >= wd_cnst*h_cruise + wd_ln,
                            V[iLoiter] >= V_wind[0],
                            V[iCruise] >= V_wind[1]
                            #h[NCruise] >= 11800*units('ft')
                           ])

        objective = 1/t_station

        Model.__init__(self,objective,constraints, **kwargs)

if __name__ == '__main__':
    M = GasPoweredHALE()
    sol = M.solve('mosek')

    if PLOT:
        M.substitutions.update({'MTOW': ('sweep', np.linspace(70, 500, 15))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

        MTOW = sol('MTOW')
        t_station = sol('t_{station}')

        plt.close()
        plt.plot(MTOW, t_station)
        plt.xlabel('MTOW [lbf]')
        plt.ylabel('t_station [days]')
        plt.grid()
        plt.axis([70, 500, 0, 15])
        plt.savefig('tvsMTOW.png')
