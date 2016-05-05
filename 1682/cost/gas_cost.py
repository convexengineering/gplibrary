''' GasPoweredHALE cost model '''
import numpy as np
import matplotlib.pyplot as plt
from gpkit import VectorVariable, Variable, Model, units
from gpkit.tools import te_exp_minus1
from gpkit.interactive import plotting
import gpkit
gpkit.settings['latex_modelname'] = False
#plt.rc('font', plt.rc('font', family='serif') 
#plt.rc('font', serif='Times New Roman')

plotMTOW = False
plotPayload = False
plotQ = False
plotT = True

class GasPoweredHALE(Model):
    '''
    Model built for a medium altitude aircraft (15,000 ft)
    Long endurance 6 days
    '''
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

        # ----------------------------------------------------
        # altitude constraints
        h_station = Variable('h_{station}', 15000, 'ft',
                             'minimum altitude at station')
        h_cruise = Variable('h_{cruise}', 5000, 'ft',
                            'minimum cruise altitude')
        h = VectorVariable(NSeg, 'h',
                           [h_cruise.value, h_cruise.value, h_station.value,
                            h_station.value, h_station.value, h_station.value,
                            h_station.value, h_station.value, h_cruise.value],
                           'ft', 'altitude')
        deltah = Variable(r'\delta_h', h_station.value-h_cruise.value, 'ft',
                          'delta height')
        t = VectorVariable(NSeg, 't', 'days', 'time per flight segment')
        h_dot = VectorVariable(NClimb, 'h_{dot}', 'ft/min', 'Climb rate')
        h_dotmin = Variable('h_{dot-min}', 100, 'ft/min',
                            'minimum necessary climb rate')

        constraints.extend([t[iClimb[0]]*h_dot[0] >= h_cruise,
                            t[iClimb[1]]*h_dot[1] >= deltah,
                            h_dot >= h_dotmin,
                           ])

        #----------------------------------------------------
        # Atmosphere model
        gamma = Variable(r'\gamma', 1.4, '-', 'Heat capacity ratio of air')
        g = Variable('g', 9.81, 'm/s^2', 'Gravitational acceleration')
        p_sl = Variable('p_{sl}', 101325, 'Pa', 'Pressure at sea level')
        rho_sl = Variable(r'\rho_{sl}', 1.225, 'kg/m^3',
                          'density at sea level')
        T_sl = Variable('T_{sl}', 288.15, 'K', 'Temperature at sea level')
        mu_sl = Variable(r'\mu_{sl}', 1.789*10**-5, 'N*s/m^2',
                         'Dynamic viscosity at sea level')
        L_atm = Variable('L_{atm}', 0.0065, 'K/m', 'Temperature lapse rate')
        T_atm = VectorVariable(NSeg, 'T_{atm}',
                               [T_sl.value - L_atm.value*v.value for v in h],
                               'K', 'Air temperature')
        a_atm = VectorVariable(NSeg, 'a_{atm}', 'm/s',
                               'Speed of sound at altitude')
        mu_atm = VectorVariable(NSeg, r'\mu', 'N*s/m^2', 'Dynamic viscosity')

        R_spec = Variable('R_{spec}', 287.058, 'J/kg/K',
                          'Specific gas constant of air')
        M_atm = Variable("M_{atm}", 0.0289644, "kg/mol",
                         "Molar mass of dry air")
        h_ref = Variable('h_{ref}', 15000, 'ft', 'ref altitude')
        R_atm = Variable("R_{atm}", 8.31447, "J/mol/K")
        rho = VectorVariable(NSeg, r'\rho',
                             [p_sl.value*x.value**(5.257-1)/R_spec.value/
                              (T_sl.value**5.257) for x in T_atm],
                             'kg/m^3', 'air density')

        # Atmospheric variation with altitude (valid from 0-7km of altitude)
        constraints.extend([#(rho/rho_sl)**0.1 == 0.954*(h/h_ref)**(-0.0284),
                            #(T_atm/T_sl)**0.1 == 0.989*(h/h_ref)**(-0.00666),
            (mu_atm/mu_sl)**0.1 == 0.991*(h/h_station)**(-0.00529)
            ])

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
        Vol_avionics = Variable('Vol_{avionics}', 0.125, 'ft^3',
                                'Avionics volume')

        # end of first segment weight + first segment fuel weight must be
        # greater than MTOW.  Each end of segment weight must be greater than
        # the next end of segment weight + the next segment fuel weight. The
        # last end segment weight must be greater than the zero fuel weight

        constraints.extend([MTOW >= W_end[0] + W_fuel[0],
                            W_end[:-1] >= W_end[1:] + W_fuel[1:],
                            W_end[-1] >= W_zfw])

        #----------------------------------------------------
        # Steady level flight model
        CD = VectorVariable(NSeg, 'C_D', '-', 'Drag coefficient')
        CL = VectorVariable(NSeg, 'C_L', '-', 'Lift coefficient')
        V = VectorVariable(NSeg, 'V', 'm/s', 'cruise speed')
        S = Variable('S', 'ft^2', 'wing area')
        eta_prop = VectorVariable(NSeg, r'\eta_{prop}', '-',
                                  'propulsive efficiency')
        eta_propCruise = Variable(r'\eta_{prop-cruise}', 0.6, '-',
                                  'propulsive efficiency in cruise')
        eta_propClimb = Variable(r'\eta_{prop-climb}', 0.5, '-',
                                 'propulsive efficiency in climb')
        eta_propLoiter = Variable(r'\eta_{prop-loiter}', 0.7, '-',
                                  'propulsive efficiency in loiter')
        P_shaft = VectorVariable(NSeg, 'P_{shaft}', 'hp', 'Shaft power')
        T = VectorVariable(NSeg, 'T', 'lbf', 'Thrust')
        P_shaftmaxV = Variable('P_{shaft-maxV}', 'hp', 'Max shaft power at max V')
        V_max = Variable('V_{max}', 'm/s', 'max velocity')

        # Climb model
        # Currently climb rate is being optimized to reduce fuel consumption.
        # In future, could implement min climb rate.

        constraints.extend([
            P_shaft == T*V/eta_prop,
            P_shaftmaxV == V_max*W_zfw/10/eta_propLoiter,
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

        W_eng = Variable('W_{eng}', 'lbf', 'engine weight')
        W_engtot = Variable('W_{eng-tot}', 'lbf', 'Installed engine weight')
        W_engref = Variable('W_{eng-ref}', 4.4107, 'lbf',
                            'Reference engine weight')
        FuelOilFrac = Variable('FuelOilFrac', 0.98, '-', 'Fuel-oil fraction')
        P_shaftref = Variable('P_{shaft-ref}', 2.295, 'hp',
                              'reference shaft power')
        BSFC_min = Variable('BSFC_{min}', 0.32, 'kg/kW/hr', 'Minimum BSFC')
        BSFC = VectorVariable(NSeg, 'BSFC', 'lb/hr/hp',
                              'brake specific fuel consumption')
        RPM_max = Variable('RPM_{max}', 9000, 'rpm', 'Maximum RPM')
        RPM = VectorVariable(NSeg, 'RPM', 'rpm', 'Engine operating RPM')
        P_shaftmaxMSL = Variable('P_{shaft-maxMSL}', 'hp',
                                 'Max shaft power at MSL')
        P_shaftmax = VectorVariable(NSeg, 'P_{shaft-max}', 'hp',
                                    'Max shaft power at altitude')
        Lfactor = VectorVariable(NSeg, 'L_factor', '-',
                                 'Max shaft power loss factor')
        P_avn = Variable('P_{avn}', 40, 'watts', 'avionics power')
        P_pay = Variable('P_{pay}', 10, 'watts', 'payload power')
        P_shafttot = VectorVariable(NSeg, 'P_{shaft-tot}', 'hp',
                                    'total power, including avionics')
        eta_alternator = Variable(r'\eta_{alternator}', 0.8, '-',
                                  'alternator efficiency')
        P_shafttotmaxV = Variable('P_{shaft-tot-maxV}', 'hp', 
                                  'total shaft power at max velocity')

        # Engine Weight Constraints
        constraints.extend([
            Lfactor == 0.906**(1/0.15)*(h/h_ref)**0.92,
            P_shaftmax/P_shaftmaxMSL + Lfactor <= 1,
            P_shaftmax >= P_shafttot,
            W_eng/W_engref >= 0.5538*(P_shaftmaxMSL/P_shaftref)**1.075,
            W_engtot >= 2.572*W_eng**0.922*units('lbf')**0.078,
            P_shafttot[iCruise] >= P_shaft[iCruise] + P_avn/eta_alternator,
            P_shafttot[iClimb] >= P_shaft[iClimb] + P_avn/eta_alternator,
            P_shafttot[iLoiter] >= (P_shaft[iLoiter] +
                                    (P_avn+P_pay)/eta_alternator),
            (BSFC/BSFC_min)**0.129 >= (2*.486*(RPM/RPM_max)**-0.141 +
                                       0.0268*(RPM/RPM_max)**9.62),
            (P_shafttot/P_shaftmax)**0.1 == 0.999*(RPM/RPM_max)**0.292,
            RPM <= RPM_max,
            P_shaftmax[iCruise] >= P_shaftmaxMSL*.81,
            P_shaftmax[iClimb[0]] >= P_shaftmaxMSL*.81,
            P_shaftmax[iClimb[1]] >= P_shaftmaxMSL*0.481,
            P_shaftmax[iLoiter] >= P_shaftmaxMSL*0.481,
            P_shaftmaxV == P_shaftmaxMSL,
            ])

        #----------------------------------------------------
        # Breguet Range
        z_bre = VectorVariable(NSeg, 'z_{bre}', '-', 'breguet coefficient')
        t_cruise = Variable('t_{cruise}', 1, 'days', 'time to station')
        t_station = Variable('t_{station}', 6, 'days', 'time on station')
        R = Variable('R', 200, 'nautical_miles', 'range to station')
        R_cruise = Variable('R_{cruise}', 180, 'nautical_miles',
                            'range to station during climb')

        constraints.extend([
            z_bre >= P_shafttot*t*BSFC*g/W_end,
            R_cruise <= V[iCruise[0]]*t[iCruise[0]],
            R <= V[iCruise[1]]*t[iCruise[1]],
            t[iLoiter] >= t_station/NLoiter,
            sum(t[[0, 1, 2]]) <= t_cruise,
            FuelOilFrac*W_fuel/W_end >= te_exp_minus1(z_bre, 3)])

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
        Re_ref = Variable('Re_{ref}', 3e5, '-', 'Reference Re for cdp')
        cdp = VectorVariable(NSeg, 'c_{dp}', '-', 'wing profile drag coeff')

        constraints.extend([
            CD >= (CDfuse + cdp + CL**2/(np.pi*e*AR))*1.3,
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
        #   CD >= CDfuse + 2*Cf*Kwing + CL**2/(pi*e*AR)
        #   + cl_16*CL**16 + CDAland,
        #   CDAland >= (2*CDland*A_rearland + CDland*A_frontland)/S])

        #----------------------------------------------------
        # Weight breakdown

        W_cent = Variable('W_{cent}', 'lbf', 'Center aircraft weight')
        W_fuse = Variable('W_{fuse}', 'lbf', 'fuselage weight')
        W_wing = Variable('W_{wing}', 'lbf', 'Total wing structural weight')
        W_fueltot = Variable('W_{fuel-tot}', 'lbf', 'total fuel weight')
        #W_skid = Variable('W_{skid}', 3, 'lbf', 'skid weight')
        m_fuse = Variable('m_{fuse}', 'kg', 'fuselage mass')
        m_cap = Variable('m_{cap}', 'kg', 'Cap mass')
        m_skin = Variable('m_{skin}', 'kg', 'Skin mass')
        m_tail = Variable('m_{tail}', 'kg', 'tail mass')
        m_rib = Variable('m_{rib}', 'kg', 'rib mass')

        constraints.extend([
            W_wing >= m_skin*g + 1.2*m_cap*g,
            m_tail*g == W_wing*0.13,
            m_rib*g == W_fuse*0.55,
            W_fuse >= m_fuse*g + m_rib*g,
            W_fueltot >= W_fuel.sum(),
            W_cent >= (W_fueltot + W_pay + W_engtot + W_fuse + W_avionics),
            W_zfw >= (W_pay + W_engtot + W_fuse + W_wing + m_tail*g +
                      W_avionics)])

        #----------------------------------------------------
        # Structural model

        # Structural parameters
        rho_skin = Variable(r'\rho_{skin}', 0.1, 'g/cm^2',
                            'Wing Skin Density')
        rho_cap = Variable(r'\rho_{cap}', 1.76, 'g/cm^3', 'Density of CF cap')
        E_cap = Variable('E_{cap}', 2e7, 'psi', 'Youngs modulus of CF cap')
        sigma_cap = Variable(r'\sigma_{cap}', 475e6, 'Pa', 'Cap stress')

        # Structural lengths
        h_spar = Variable('h_{spar}', 'm', 'Spar height')
        t_cap = Variable('t_{cap}', 0.028, 'in', 'Spar cap thickness')
        #arbitrarily placed based on available cf
        w_cap = Variable('w_{cap}', 'in', 'Spar cap width')
        c = Variable('c', 'ft', 'Wing chord')
        #assumes straight, untapered wing

        # Structural ratios
        tau = Variable(r'\tau', 0.12, '-', 'Airfoil thickness ratio')
        LoverA = Variable('LoverA', 'lbf/ft^2', 'Wing loading')
        lambda_c = Variable(r'\lambda_c', '-', 'Taper ratio')

        # Structural areas
        A_capcent = Variable('A_{capcent}', 'm**2', 'Cap area at center')
        A_cap = Variable('A_{cap}', 'm**2', 'Cap area')
        #currently assumes constant area

        # Structural volumes
        Vol_cap = Variable('Vol_{cap}', 'm**3', 'Cap volume')

        # Structural evaluation parameters
        M_cent = Variable('M_cent', 'N*m', 'Center bending moment')
        F = Variable('F', 'N', 'Load on wings')
        SL = Variable('SL', 'Pa', 'Shear load')
        N_Max = Variable('N_{Max}', 5, '-', 'Load factor')
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

        rho_fuel = Variable(r'\rho_{fuel}', 6.01, 'lbf/gallon',
                            'density of 100LL')
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

        constraints.extend([V_wind[0] >= wd_cnst*h_station + wd_ln,
                            V_wind[1] >= wd_cnst*h_cruise + wd_ln,
                            V[iLoiter] >= V_wind[0],
                            V[iCruise] >= V_wind[1]
                           ])

        #----------------------------------------------------
        # Cost Model
        Q = Variable('Q', 10, 'count', 'Number produced in 5 years')
        FTA = Variable('FTA', 1, 'count', 'Number of Flight Test Aircraft')
        N_eng = Variable('N_{eng}', 'count',
                    'Number of engines or Q*num of engines per aircraft')
        R_avn = Variable('R_{avn}', 5000, 'USD2012/lbf', 'Avionics Cost')

        # Constants Hourly Rates
        R_E = Variable('R_E', 115, 'USD2012/hr', 'Enginering Hourly Rate')
        R_T = Variable('R_T', 118, 'USD2012/hr', 'Tooling Hourly Rate')
        R_M = Variable('R_M', 108, 'USD2012/hr',
                       'Manufacturing Hourly Rate')
        R_Q = Variable('R_Q', 98, 'USD2012/hr', 'Quality Check Hourly Rate')

        # Free Variableiables

        # Hourly costs in hrs. Multiply by hourly work rate to get cost
        H_E = Variable('H_E', 'hrs',
                  'Engineering hours of airframe and component integration')
        # Eng Hours does not include cost of avionics or engine engineering
        H_T = Variable('H_T', 'hrs',
                  'Tooling hours for production preparation')
        # H_T also covers tooling cost through ongoing production
        H_M = Variable('H_M', 'hrs',
                  'Manufacturing hours of main and subcontractors')
        H_Q = Variable('H_Q', 'hrs', 'Quality control hours for inspection')

        # Cost Variableibles
        C_Dev = Variable('C_{Dev}', 'USD2012', 'Development Cost')
        C_F = Variable('C_F', 'USD2012',
                       'Flight Test Cost to prove airworthiness')
        C_M = Variable('C_M', 'USD2012',
                       'Materials cost of aluminum airframe')
        C_avn = Variable('C_{avn}', 'USD2012', 'Cost of avionics')
        C_fly = Variable('C_{fly}', 'USD2012', 'Flyaway cost')
        C_plane = Variable('C_{plane}', 'USD2012', 'cost per plane')
        #Raymer suggests C_avn = 0.15*C_fly

        C_eng = Variable('C_{eng}', 24000, 'USD2012', 'Engine Cost')
        H_E_const = Variable('H_{E_const}', 4.86,
                             units.hr/(units.lbf**0.777 * units.knots**0.894),
                             'H_E power law proportionality constant')
        H_T_const = Variable('H_{T_const}', 5.99,
                             units.hr/(units.lbf**0.777 * units.knots**0.696),
                             'H_T power law proportionality constant')
        H_M_const = Variable('H_{M_const}', 7.37,
                             units.hr/(units.lbf**0.82 *  units.knots**0.484),
                             'H_M power law proportionality constant')
        H_Q_const = Variable('H_{Q_const}', 0.133, '-',
                             'H_Q power law proportionality constant')
        C_Dev_const = Variable('C_{Dev_const}', 91.3,
                        units.USD2012/(units.lbf**0.630 * units.knots**1.3),
                        'C_Dev power law proportionality constant')
        C_F_const = Variable('C_{F_const}', 2498,
                        units.USD2012/(units.lbf**0.325 * units.knots**0.822),
                        'C_F power law proportionality constant')
        C_M_const = Variable('C_{M_const}', 22.1,
                        units.USD2012/(units.lbf**0.921 * units.knots**0.621),
                        'C_M power law proportionality constant')

        constraints.extend([
            C_fly >= (H_E*R_E + H_T*R_T + H_M*R_M + H_Q*R_Q + C_Dev + C_F + C_M +
                     C_eng*N_eng + C_avn),
            C_plane == C_fly/Q,
            H_E == H_E_const*W_zfw**0.777*V_max**0.894*Q**0.163,
            H_T == H_T_const*W_zfw**0.777*V_max**0.696*Q**0.263,
            H_M == H_M_const*W_zfw**0.82*V_max**0.484*Q**0.641,
            H_Q == H_Q_const*H_E,
            C_Dev == C_Dev_const*W_zfw**0.630*V_max**1.3,
            C_F == C_F_const*W_zfw**0.325*V_max**0.822*FTA**1.21,
            C_M == C_M_const*W_zfw**0.921*V_max**0.621*Q**0.799,
            C_avn == R_avn*W_zfw,
            N_eng == Q*2
            ])

        objective = C_fly
        Model.__init__(self, objective, constraints)

if __name__ == '__main__':
    M = GasPoweredHALE()
    sol = M.solve('mosek')

    plt.rcParams.update({'font.size':16})
    plt.rc('font', family='serif') 
    plt.rc('font', serif='Times New Roman')

    if plotMTOW:
        M.substitutions.update({'MTOW':('sweep', np.linspace(100, 500, 50))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

        MTOW = sol('MTOW')
        C_fly = sol('C_{fly}')
        C_plane = sol('C_{plane}')
        t_station = sol('t_{station}')
        
        plt.close()
        plt.plot(MTOW, C_fly/1e6)
        #plt.title('Aircraft Weight vs Flyaway Cost')
        plt.xlabel('Mass Take Off Weight [lbs]')
        plt.ylabel('Flyaway Cost [Millions of $]')
        plt.axis([100,500,0,15])
        plt.grid()
        plt.savefig('MTOWvsC_fly.pdf')
        
        plt.close()
        plt.plot(MTOW, C_plane/1e6)
        #plt.rcParams.update({'font.size':16})
        #plt.title('Aircraft Weight vs Cost per Plane')
        plt.xlabel('Mass Take Off Weight [lbs]')
        plt.ylabel('Cost per plane [Millions of $]')
        plt.axis([100,500,0,1])
        plt.grid()
        plt.savefig('MTOWvsC_plane.pdf')
        
        plt.close()
        plt.plot(MTOW, t_station)
        #plt.rcParams.update({'font.size':16})
        plt.xlabel('MTOW [lbs]')
        plt.ylabel('time on station [days]')
        plt.grid()
        plt.savefig('MTOWvst_station.pdf')

    if plotPayload:

        M.substitutions.update({'t_{station}': 6})
        M.cost = M["MTOW"] + M["C_{fly}"]*units('lbf/USD2012')
        M.substitutions.update({'W_{pay}':('sweep', np.linspace(5, 40, 10))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)
    
        W_pay = sol('W_{pay}')
        C_fly = sol('C_{fly}')
        C_plane = sol('C_{plane}')
        
        plt.close()
        plt.rcParams.update({'font.size':16})
        plt.plot(W_pay, C_fly/1e6)
        #plt.title('Payload Weight vs Flyaway Cost')
        plt.xlabel('Payload Weight [lbs]')
        plt.ylabel('Flyaway Cost [Millions of $]')
        plt.axis([5,40,0,15])
        plt.grid()
        plt.savefig('W_payvsC_fly.pdf')
        
        plt.close()
        plt.rcParams.update({'font.size':16})
        plt.plot(W_pay, C_plane/1e6)
        #plt.title('Payload Weight vs Cost per Plane')
        plt.xlabel('Payload weight [lbs]')
        plt.ylabel('Cost per plane [Millions of $]')
        plt.axis([5,40,0,1])
        plt.grid()
        plt.savefig('W_payvsC_plane.pdf')

    if plotQ:

        M.substitutions.update({'t_{station}': 6})
        M.cost = M["MTOW"] + M["C_{fly}"]*units('lbf/USD2012')
        M.substitutions.update({'Q':('sweep', np.linspace(1, 50, 50))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

        Q = sol('Q')
        C_fly = sol('C_{fly}')
        C_plane = sol('C_{plane}')
        
        plt.close()
        plt.rcParams.update({'font.size':16})
        plt.plot(Q, C_fly/1e6)
        #plt.title('Number of Aircraft vs Flyaway Cost')
        plt.xlabel('Number of Aircraft')
        plt.ylabel('Flyaway Cost [Millions of $]')
        plt.axis([1,20,0,15])
        plt.grid()
        plt.savefig('QvsC_fly.pdf')
        
        plt.close()
        plt.rcParams.update({'font.size':16})
        plt.plot(Q, C_plane/1e6)
        #plt.title('Number of Aircraft vs Cost per Plane')
        plt.xlabel('Number of Aircraft')
        plt.ylabel('Cost per plane [Millions of $]')
        plt.axis([1,20,0,5])
        plt.grid()
        plt.savefig('QvsC_plane.pdf')
    
    if plotT:

        M.substitutions.update({'t_{station}':('sweep', np.linspace(0.05, 10, 30))})
        sol = M.solve(solver='mosek', verbosity=0, skipsweepfailures=True)

        t_station = sol('t_{station}')
        C_fly = sol('C_{fly}')
        C_plane = sol('C_{plane}')
        
        plt.close()
        plt.rcParams.update({'font.size':16})
        plt.plot(t_station, C_fly/1e6)
        #plt.title('Time on Station vs Flyaway Cost')
        plt.xlabel('Time on Station [days]')
        plt.ylabel('Flyaway Cost [Millions of $]')
        plt.axis([0,10,0,15])
        plt.grid()
        plt.savefig('t_stationvsC_fly.pdf')
        
        plt.close()
        plt.rcParams.update({'font.size':16})
        plt.plot(t_station, C_plane/1e6)
        #plt.title('Time on Station vs Cost per Plane')
        plt.xlabel('Time on Station [days]')
        plt.ylabel('Cost per plane [Millions of $]')
        plt.axis([0,10,0,1])
        plt.grid()
        plt.savefig('t_stationvsC_plane.pdf')
