import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.small_scripts import mag
from gpkit.constraints.tight import TightConstraintSet as TCS

class FanAndLPC(Model):
    """
    Free Stream, Fan, and LPC Calcs for the Engine Model
    """
    def __init__(self, **kwargs):
        #air properties
        Rair = Variable('R_{air}', 287, 'J/kg/K', 'R for Air')
        gammaAir = Variable('gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Air')
        Cpair = Variable('Cp_{air}', 1068, 'J/kg/K', "Cp Value for Air at 673K") #http://www.engineeringtoolbox.com/air-properties-d_156.html
        gm1 = Variable('gm1', 0.4, '-', 'Gamma for Air Minus 1')
        c1 = Variable('c1', 1.128, '-', 'Constant in Stagnation Eqn')

        #free stream
        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
        P0 = Variable('P_0', 'kPa', 'Free Stream Static Pressure')
        M0 = Variable('M_0', '-', 'Free Stream Mach Number')
        
        #new Vars for free stream
        a0 = Variable('a_0', 'm/s', 'Speed of Sound in Freestream')
        u0 = Variable('u_0', 'm/s', 'Free Stream Speed')
        Pt0 = Variable('P_{t_0}', 'kPa', 'Free Stream Stagnation Pressure')
        Tt0 = Variable('T_{t_0}', 'K', 'Free Stream Stagnation Temperature')
        ht0 = Variable('h_{t_0}', 'J/kg', 'Free Stream Stagnation Enthalpy')

        #new vars for the diffuser exit
        Pt18 = Variable('P_{t_1.8}', 'kPa', 'Stagnation Pressure at the Diffuser Exit (1.8)')
        Tt18 = Variable('T_{t_1.8}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')
        ht18 = Variable('h_{t_1.8}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
        
        #new vars for the fan inlet
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')

        #new vars for the fan exit (station 2.1)
        Pt21 = Variable('P_{t_2.1}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2.1)')
        Tt21 = Variable('T_{t_2.1}', 'K', 'Stagnation Temperature at the Fan Inlet (2.1)')
        ht21 = Variable('h_{t_2.1}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2.1)')

        #new vars for the LPC exit (station 2.5)
        Pt25 = Variable('P_{t_2.5}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_2.5}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')

        #HPC exit state variables (station 3)
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')

        #new vars for the fan nozzle exit (station 7)
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
        ht7 = Variable('h_{t_7}', 'J/kg', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')
        
        #fan efficiency variable
        etaFan = Variable('\eta_{fan}', 0.9,'-', 'Fan Polytropic Efficiency')

        #compressor efficiency variables
        etalc = Variable('\eta_{lc}', 0.887,'-', 'Low Pressure Compressor Polytropic Efficiency')
        etahc = Variable('\eta_{hc}', 0.887,'-', 'High Pressure Compressor Polytropic Efficiency')

        #pressure ratios
        pid = Variable('\pi_{d}', '-', 'Diffuser Pressure Ratio')
        pif = Variable('\pi_f', '-', 'Fan Pressure Ratio')
        pifn = Variable('\pi_{fn}', '-', 'Fan Duct Pressure Loss Ratio')
        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')
        
        #constraints
        constraints = [
            #free stream speed of sound and free stream velocity
            a0 == (gammaAir * Rair*  T0)**.5,        #traceable to B.109
            u0 == M0 * a0,                           #traceable to B.110
            
            #free stream stagnation values
            Pt0 == P0 / (c1 ** -3.5), #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            Tt0 == T0 / (c1) ** (-1),             #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            ht0 == Cpair * Tt0,

            #diffuser exit stagnation calues (station 1.8)
            Pt18 == pid * Pt0,  #B.113
            Tt18 == Tt0,        #B.114
            ht18 == ht0,        #B.115

            #fan inlet constraints (station 2)
            Tt2 == Tt18,    #B.120
            ht2 == ht18,    #B.121
            Pt2 == Pt18,
                        
            #fan exit constraints (station 2.1)
            Pt21 == pif * Pt2,  #16.50
            Tt21 == Tt2 * pif ** (.31746),   #16.50
            ht21 == Cpair * Tt21,   #16.50

            #fan nozzle exit (station 7)
            Pt7 == pifn * Pt21,     #B.125
            Tt7 == Tt21,    #B.126
            ht7 == ht21,    #B.127

            #LPC exit (station 2.5)
            Pt25 == pilc * Pt2,
            Tt25 == Tt2 * pilc ** (.3221),
            ht25 == Tt25 * Cpair,

            #HPC Exit
            Pt3 == pihc * Pt25,
            Tt3 == Tt25 * pihc ** (.3221),
            ht3 == Cpair * Tt3
            ]
        Model.__init__(self, ht3, constraints, **kwargs)

class CombustorCooling(Model):
    """
    class to represent the engine's combustor and perform calculations
    on engine cooling bleed flow...cooling flow is currently not implemented
    """
    def __init__(self, **kwargs):
        #new vars
        #gas properties
        gammaC =  Variable('gamma_{C}', 1.4, '-', 'Specific Heat Ratio for Gas in Combustor')
        Cpc = Variable('Cp_c', 1068, 'J/kg/K', "Cp Value for Fuel/Air Mix in Combustor")
        
        #HPC exit state variables (station 3)
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')
        
        #combustor exit state variables..recall Tt4 is already set in Engine class
        Pt4 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Combustor Exit (4)')
        ht4 = Variable('h_{t_4}', 'J/kg', 'Stagnation Enthalpy at the Combustor Exit (4)')
        Tt4 = Variable('T_{t_4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')

        #Turbine inlet state variables (station 4.1)
        Pt41 = Variable('P_{t_4.1}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_4.1}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')

        #burner pressure ratio
        pib = Variable('\pi_{b}', '-', 'Burner Pressure Ratio')

        #flow faction f
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        #heat of combustion of jet fuel
        hf = Variable('h_f', 42.8, 'MJ/kg', 'Heat of Combustion of Jet Fuel')     #http://hypertextbook.com/facts/2003/EvelynGofman.shtml...prob need a better source

        with SignomialsEnabled():
        
            constraints = [
                #flow through combustor
                Pt4 == pib * Pt3,   #B.145
                ht4 == Cpc * Tt4,

                #fuel flow fraction f
                f*hf >= ht4 - ht3,
                
                #flow at turbine inlet
                Tt41 == Tt4,
                Pt41 == Pt4,
                ht41 == ht4
                ]
            
        Model.__init__(self, f, constraints, **kwargs)

class Turbine(Model):
    """
    classs to represent the engine's turbine
    """
    def __init__(self, **kwargs):
        #gas properties
        gammaT = Variable('gamma_{T}', 1.4, '-', 'Specific Heat Ratio for Gas in Turbine')
        Cpt =Variable('Cp_t', 1068, 'J/kg/K', "Cp Value for Combustion Products in Turbine")
        
        #new variables
        ht45 = Variable('h_{t_4.5}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (4.5)')
        Pt45 = Variable('P_{t_4.5}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_4.5}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

        #enthalpies used in shaft power balances
        ht18 = Variable('h_{t_1.8}', 'J/kg', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')
        ht21 = Variable('h_{t_2.1}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2.1)')
        ht25 = Variable('h_{t_2.5}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')
        ht3 = Variable('h_{t_3}', 'J/kg', 'Stagnation Enthalpy at the HPC Exit (3)')

        #Turbine inlet state variables (station 4.1)
        Pt41 = Variable('P_{t_4.1}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4.1}', 'K', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_4.1}', 'J/kg', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')

        #HPT rxit states
        Pt49 = Variable('P_{t_4.9}', 'kPa', 'Stagnation Pressure at the HPTExit (49)')
        Tt49 = Variable('T_{t_4.9}', 'K', 'Stagnation Temperature at the HPT Exit (49)')
        ht49 = Variable('h_{t_4.9}', 'J/kg', 'Stagnation Enthalpy at the HPT Exit (49)')
        
        #turbine nozzle exit states
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J/kg', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')

        #pressure ratios
        pihpt = Variable('\pi_{HPT}', '-', 'HPT Pressure Ratio')
        pilpt = Variable('\pi_{LPT}', '-', 'LPT Pressure Ratio')
        
        #turbine efficiences
        etaht = Variable('\eta_{ht}', 0.9, '-', 'Polytropic Efficiency of HPT')
        etalt = Variable('\eta_{lt}', 0.9, '-', 'Polytropic Efficiency of LPT')

        #flow faction f
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        #BPR
        alpha = Variable('alpha', '-', 'By Pass Ratio')

        #relavent pressure ratio
        pitn = Variable('\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')

        with SignomialsEnabled():
            constraints = [
                #HPT shafter power balance

                #SIGNOMIAL
                TCS([(1+f)*(ht41-ht45) >= ht3 - ht25]),    #B.161
                
                #LPT shaft power balance

                #SIGNOMIAL  
                TCS([(1+f)*(ht49 - ht45) <= -((ht25-ht18)+alpha*(ht21 - ht2))]), #B.165

                
                #HPT Exit states (station 4.5)
                Pt45 == pihpt * Pt41,
                pihpt == (Tt45/Tt41)**(3.889),
                ht45 == Cpt * Tt45,

                #LPT Exit States
                Pt49 == pilpt * Pt45,
                pilpt == (Tt49/Tt45)**(3.889),
                ht49 == Cpt * Tt49,

                #turbine nozzle exit states
                Pt5 == pitn * Pt49, #B.167
                Tt5 == Tt49,    #B.168
                ht5 == ht49     #B.169
                ]
        Model.__init__(self, 1/ht49, constraints, **kwargs)

class ExhaustAndThrust(Model):
    """
    Class to calculate the exhaust quantities as well as
    the overall engine thrust and on design TSFC & ISP
    """
    def __init__(self, **kwargs):
        
        a0 = Variable('a_0', 'm/s', 'Speed of Sound in Freestream')
        u0 = Variable('u_0', 'm/s', 'Free Stream Speed')
        
        #external atmosphere properties
        P0 = Variable('P_0', 'kPa', 'Free Stream Stagnation Pressure')
        
        #gas properties
        Cpt =Variable('Cp_t', 1068, 'J/kg/K', "Cp Value for Combustion Products in Turbine")
        gammaT = Variable('gamma_{T}', 1.4, '-', 'Specific Heat Ratio for Gas in Turbine')
        gammaAir = Variable('gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Air')
        Cpair = Variable('Cp_{air}', 1068, 'J/kg/K', "Cp Value for Air at 673K") #http://www.engineeringtoolbox.com/air-properties-d_156.html
        
        #gravity
        g = Variable('g', 9.81, 'm/(s^2)', 'Gravitational Acceleration')
        
        #fan exhaust variables
        P8 = Variable('P_8', 'kPa', 'Fan Exhaust Static Pressure')
        Pt8 = Variable('P_{t_8}', 'kPa', 'Fan Exhaust Stagnation Pressure')
        ht8 = Variable('h_{t_8}', 'J/kg', 'Fan Exhaust Stagnation Enthalpy')
        h8 = Variable('h_8', 'J/kg', 'Fan Exhasut Static Enthalpy')
        Tt8 = Variable('T_{t_8}', 'K', 'Fan Exhaust Stagnation Temperature (8)')
        T8 = Variable('T_{8}', 'K', 'Fan Exhaust Sttic Temperature (8)')

        #turbine nozzle exit states
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'K', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J/kg', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')
        
        #core exit variables
        P6 = Variable('P_6', 'kPa', 'Core Exhaust Static Pressure')
        Pt6 = Variable('P_{t_6}', 'kPa', 'Core Exhaust Stagnation Pressure')
        Tt6 = Variable('T_{t_6}', 'K', 'Core Exhaust Stagnation Temperature (6)')
        T6 = Variable('T_{6}', 'K', 'Core Exhaust Static Temperature (6)')
        ht6 = Variable('h_{t_6}', 'J/kg', 'Core Exhaust Stagnation Enthalpy')
        h6 = Variable('h_6', 'J/kg', 'Core Exhasut Static Enthalpy')

        #new vars for the fan nozzle exit (station 7)
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
        ht7 = Variable('h_{t_7}', 'J/kg', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')

        #thrust variables
        F8 = Variable('F_8', 'N', 'Fan Thrust')
        F6 = Variable('F_6', 'N', 'Core Thrust')
        F = Variable('F', 'N', 'Total Thrust')
        Fsp = Variable('F_{sp}', '-', 'Specific Net Thrust')
        Isp = Variable('I_{sp}', 's', 'Specific Impulse')
        TSFC = Variable('TSFC', '1/hr', 'Thrust Specific Fuel Consumption')

        #exhaust speeds
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')

        #BPR
        alpha = Variable('alpha', '-', 'By Pass Ratio')
        alphap1 = Variable('alphap1', '-', '1 plus BPR')

        #core mass flow
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')

        #flow faction f
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        with SignomialsEnabled():
            constraints = [
                #fan exhaust
                Pt8 == Pt7, #B.179
                Tt8 == Tt7, #B.180
                P8 == P0,
                h8 == Cpair * T8,
                TCS([u8**2 + 2*h8 <= 2*ht8]),
                (P8/Pt8)**(.2857) == T8/Tt8,
                ht8 == Cpair * Tt8,
                
                #core exhaust
                P6 == P0,   #B.4.11 intro
                Pt6 == Pt5, #B.183
                Tt6 == Tt5, #B.184
                (P6/Pt6)**(.2857) == T6/Tt6,
                TCS([u6**2 + 2*h6 <= 2*ht6]),
                h6 == Cpt * T6,
                ht6 == Cpt * Tt6,

                #overall thrust values
                TCS([F8/(alpha * mCore) + u0 <= u8]),  #B.188
                TCS([F6/mCore + u0 <= (1+f)*u6]),      #B.189

                #SIGNOMIAL
                TCS([F <= F6 + F8]),

                Fsp == F/((alphap1)*mCore*a0),   #B.191

                #ISP
                Isp == Fsp*a0*(alphap1)/(f*g),  #B.192

                #TSFC
                TSFC == 1/Isp                   #B.193
                ]
        Model.__init__(self, TSFC, constraints, **kwargs)
        
class OnDesignSizing(Model):
    """
    class to perform the on design sizing of the engine.
    Nozzle areas are calcualted in post processing due to dependence on
    M6 and M8 being greater than or less than 1
    """
    def __init__(self, **kwargs):
        #new variables
        a0 = Variable('a_0', 'm/s', 'Speed of Sound in Freestream')
        
        #air properties
        Cpair = Variable('Cp_{air}', 1068, 'J/kg/K', "Cp Value for Air at 673K") #http://www.engineeringtoolbox.com/air-properties-d_156.html
        Rair = Variable('R_{air}', 287, 'J/kg/K', 'R for Air')
        Cpt =Variable('Cp_t', 1068, 'J/kg/K', "Cp Value for Combustion Products in Turbine")
        Rt = Variable('R_t', 287, 'J/kg/K', 'R for the Turbine Gas')
        
        #fan face variables
        hold2 = Variable('hold_{2}', '-', '1+(gamma-1)/2 * M_2**2')
        rho2 = Variable('\rho_2', 'kg/m^3', 'Air Static Density at Fan Face')
        T2 = Variable('T_2', 'K', 'Air Static Temperature at Fan Face')
        P2 = Variable('P_2', 'kPa', 'Air Static Pressure at Fan Face')
        u2 = Variable('u_2', 'm/s', 'Air Speed at Fan Face')
        A2 = Variable('A_2', 'm^2', 'Fan Area')

        #HPC face variables
        hold25 = Variable('hold_{2.5}', '-', '1+(gamma-1)/2 * M_2.5**2')
        rho25 = Variable('\rho_2.5', 'kg/m^3', 'Static Air Density at HPC Face')
        T25 = Variable('T_{2.5}', 'K', 'Static Air Temperature at HPC Face')
        P25 = Variable('P_{2.5}', 'kPa', 'Static Air Pressure at HPC Face')
        Pt25 = Variable('P_{t_2.5}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        u25 = Variable('u_{2.5}', 'm/s', 'Air Speed at HPC Face')
        A25 = Variable('A_{2.5}', 'm^2', 'HPC Area')
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        ht25 = Variable('h_{t_2.5}', 'J/kg', 'Stagnation Enthalpy at the LPC Exit (2.5)')
        h25 = Variable('h_{2.5}', 'J/kg', 'Static Enthalpy at the LPC Exit (2.5)')

        #mach numbers
        M8 = Variable('M_8', '-', 'Fan Exhaust Mach Number')
        M6 = Variable('M_6', '-', 'Core Exhaust Mach Number')

        #core mass flow
        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')

        #variables specified for sizing purposes
        Fd = Variable('F_D', 'N', 'Design Thrust')
        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')
        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')

        #thrust variables
        Fsp = Variable('F_{sp}', '-', 'Specific Net Thrust')

        #BPR
        alphap1 = Variable('alphap1', '-', '1 plus BPR')

        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')
        ht2 = Variable('h_{t_2}', 'J/kg', 'Stagnation Enthalpy at the Fan Inlet (2)')
        T2 = Variable('T_{2}', 'K', 'Static Temperature at the Fan Inlet (2)')
        h2 = Variable('h_{2}', 'J/kg', 'Static Enthalpy at the Fan Inlet (2)')

        #exhaust speeds
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')

        #exhaust temperatures
        T6 = Variable('T_{6}', 'K', 'Core Exhaust Static Temperature (6)')
        T8 = Variable('T_{8}', 'K', 'Fan Exhaust Sttic Temperature (8)')
        
        constraints = [
            #mass flow sizing
            mCore == Fd/(Fsp*a0*(alphap1)),  #B.194

            #component area sizing
            #fan area
            P2 == Pt2*(hold2)**(-3.5),
            T2 == Tt2 * hold2**-1,
            h2 == Cpair * T2,
            rho2 == P2/(Rair * T2),  #B.196
            u2 == M2*(Cpair*Rair*T2/(781*units('J/kg/K')))**.5,  #B.197
            A2 == (alphap1)*mCore/(rho2*u2),     #B.198

            #HPC area
            P25 == Pt25*(hold25)**(-3.5),
            T25 == Tt25 * hold25**-1,
            h25 == Cpair * T25,
            rho25 == P25/(Rair*T25),
            u25 == M25*(Cpair*Rair*T25/(781*units('J/kg/K')))**.5,   #B.202
            A25 == (alphap1)*mCore/(rho25*u25),     #B.203

            #mach nubmers for post processing of the data
            M8 == u8/((T8*Cpair*Rair/(781*units('J/kg/K')))**.5),
            M6 == u6/((T6*Cpt*Rt/(781*units('J/kg/K')))**.5),
            ]
        #objective is None because all constraints are equality so feasability region is a
        #single point which likely will not solve
        Model.__init__(self, None, constraints, **kwargs)

class NozzleSizing(Model):
    """
    class for post processing of solution data, computes the nozzle sizes

    Input: solution array from engine on design solve
    Output: Nozzle Sizes
    """
    def __init__(self, sol):
        self.sizing(self, sol)

    def sizing(self, sol):
        #first size the fan nozzle

        #size the core nozzle

        #return the two sizes
        return A6, A8
