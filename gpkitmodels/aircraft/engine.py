#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS

#TODO
#figure out how to compute f
#determine the three different Cp, gamma, and R values to be used
#implement the constraints on the nozzle sizing
#get realisitc R and Cp values for the different engine sections

class Engine(Model):
    """
    Engine Sizing

    References:
    1. TASOPT Volume 2: Appendicies dated 31 March 2010
    (comments B.### refer to eqn numbers in TASOPT Appendix
    section B)
    2. Flight Vehicle Aerodynamics (Prof Drela's 16.110 text book)
    3. https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
    4. Prof Barret's Spring 2015 16.50 Notes

    All engine station definitions correltate to TASOPT figure B.1

    Fan efficiency is taken from TASOPT table B.1
    Compressor efficiency is taken from TASOPT table B.2
    """

    def __init__(self, **kwargs):
        #physical constants
        #Fix the R values
        Rair = Variable('R_{air}', 287, 'J/kg/K', 'R for Air')
        Rt = Variable('R_t', 287, 'J/kg/K', 'R for the Turbine Gas')
        Rc = Variable('R_c', 287, 'J/kg/K', 'R for the Combustor Gas')
        gammaAir = Variable('gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Air')
        gammaT = Variable('gamma_{T}', 1.4, '-', 'Specific Heat Ratio for Gas in Turbine')
        gammaC =  Variable('gamma_{C}', 1.4, '-', 'Specific Heat Ratio for Gas in Combustor')
        g = Variable('g', 9.81, 'm/s/s', 'Gravitational Acceleration')

        #NEED TO FIX THESE VALUES
        Cpair = Variable('Cp_{air}', 1068, "Cp Value for Air at 673K") #http://www.engineeringtoolbox.com/air-properties-d_156.html
        Cpc = Variable('Cp_c', 1068, "Cp Value for Fuel/Air Mix in Combustor")
        Cpt =Variable('Cp_t', 1068, "Cp Value for Combustion Products in Turbine")

        #variables linking subclasses
        #enthalpies used in shaft power balances
        ht18 = Variable('h_{t_1.8}', 'J', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
        ht2 = Variable('h_{t_2}', 'J', 'Stagnation Enthalpy at the Fan Inlet (2)')
        ht21 = Variable('h_{t_2.1}', 'J', 'Stagnation Enthalpy at the Fan Inlet (2.1)')
        ht25 = Variable('h_{t_2.5}', 'J', 'Stagnation Enthalpy at the LPC Exit (2.5)')
        
        #HPC exit state variables (station 3)
        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
        ht3 = Variable('h_{t_3}', 'J', 'Stagnation Enthalpy at the HPC Exit (3)')

        #Turbine inlet state variables (station 4.1)
        Pt41 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
        Tt41 = Variable('T_{t_4}', 'J', 'Stagnation Temperature at the Turbine Inlet (4.1)')
        ht41 = Variable('h_{t_4}', 'J', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')

        #HPC exit states
        Pt45 = Variable('P_{t_4.5}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
        Tt45 = Variable('T_{t_4.5}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')

        #turbine nozzle exit states
        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
        Tt5 = Variable('T_{t_5}', 'J', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
        ht5 = Variable('h_{t_5}', 'J', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')

        #new vars for the fan nozzle exit (station 7)
        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
        ht7 = Variable('h_{t_7}', 'J', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')

        #fuel
        #flow faction f
        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')

        #exhaust speeds
        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')
        
        #Variables taken as known, substitutded for later
        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
        P0 = Variable('P_0', 'kPa', 'Free Stream Stagnation Pressure')
        M0 = Variable('M_0', '-', 'Free Stream Mach Number')
        Tt4 = Variable('T_{t4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')
        pif = Variable('\pi_f', '-', 'Fan Pressure Ratio')
        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')
        pitn = Variable('\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')
        pid = Variable('\pi_{d}', '-', 'Diffuser Pressure Ratio')
        pib = Variable('\pi_{b}', '-', 'Burner Pressure Ratio')
        pifn = Variable('\pi_{fn}', '-', 'Fan Duct Pressure Loss Ratio')
        alpha = Variable('alpha', '-', 'By Pass Ratio')
        M4a = Variable('M_{4a}', '-', 'Rep Mach # at Start of GPT Cooling Flow')


        #known variables needed only for sizing purposes
        Fd = Variable('F_D', 'N', 'Design Thrust')
        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')
        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')

        #set up the overeall model for an on design solve
        with SignomialsEnabled():
            
            lpc = FanAndLPC()
            combustor = CombustorCooling()
            turbine = Turbine()
            thrust = ExhaustAndThrust()
            size = OnDesignSizing()

            substitutions = {
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            'T_{t_4}': 1400,
            '\pi_f': 1.5,
            '\pi_{lc}': 3,
            '\pi_{hc}': 10,
            '\pi_{tn}': 1,
            '\pi_{d}': 1,
            '\pi_{b}': 1,
            '\pi_{fn}': 1,
            'alpha': 10,
            'M_{4a}': 1,    #choked turbines
            'F_D': 121436.45, #737 max thrust in N
            'M_2': .4,
            'M_{2.5}': .5 
            }

            #temporary objective is to minimize the fan size
            model = Model(A2, LinkedConstraintSet([lpc.constriants, combustor.constraints, turbine.constraints,
                                                   thrust.constraints, size.constraints]))

            design_sol = model.localsolve(verbosity = 1)

            

class FanAndLPC(Engine):
    """
    Free Stream, Fan, and LPC Calcs for the Engine Model
    """
    def __init__(self):
        #new Vars for free stream
        a0 = Variable('a_0', 'm/s', 'Speed of Sound in Freestream')
        U0 = Variable('U_0', 'm/s', 'Free Stream Speed')
        Pt0 = Variable('P_{t_0}', 'kPa', 'Free Stream Stagnation Pressure')
        Tt0 = Variable('T_{t_0}', 'K', 'Free Stream Stagnation Temperature')
        ht0 = Variable('h_{t_0}', 'J', 'Free Stream Stagnation Enthalpy')

        #new vars for the diffuser exit
        Pt18 = Variable('P_{t_1.8}', 'kPa', 'Stagnation Pressure at the Diffuser Exit (1.8)')
        Tt18 = Variable('T_{t_1.8}', 'K', 'Stagnation Temperature at the Diffuser Exit (1.8)')

        #new vars for the fan inlet
        Pt2 = Variable('P_{t_2}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2)')
        Tt2 = Variable('T_{t_2}', 'K', 'Stagnation Temperature at the Fan Inlet (2)')

        #new vars for the fan exit (station 2.1)
        Pt21 = Variable('P_{t_2.1}', 'kPa', 'Stagnation Pressure at the Fan Inlet (2.1)')
        Tt21 = Variable('T_{t_2.1}', 'K', 'Stagnation Temperature at the Fan Inlet (2.1)')

        #new vars for the LPC exit (station 2.5)
        Pt25 = Variable('P_{t_2.5}', 'kPa', 'Stagnation Pressure at the LPC Exit (2.5)')
        Tt25 = Variable('T_{t_2.5}', 'K', 'Stagnation Temperature at the LPC Exit (2.5)')
        
        #fan efficiency variable
        etaFan = Variable('\eta_{fan}', 0.9,'-', 'Fan Polytropic Efficiency')

        #compressor efficiency variables
        etalc = Variable('\eta_{lc}', 0.887,'-', 'Low Pressure Compressor Polytropic Efficiency')
        etahc = Variable('\eta_{hc}', 0.887,'-', 'High Pressure Compressor Polytropic Efficiency')
        
        #constraints
        constraints = [
            #free stream speed of sound and free stream velocity
            a0 == (gammaAir * Rair*  T0)**.5,        #traceable to B.109
            u0 == M0 * a0,                           #traceable to B.110
            
            #free stream stagnation values
            Pt0 == P0 / (1+.5*(gamma-1)*M0**2) ** (1/(gammaAir-1)), #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
            Tt0 == T0 / (1+.5*(gamma-1)*M0**2) ** (-1),             #https://www.grc.nasa.gov/www/k-12/airplane/isentrop.html
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
            Tt21 == Tt2 * pif ** ((gammaAir-1)/(gammaAir * etaFan)),   #16.50
            ht21 == Cpair * Tt21,   #16.50

            #fan nozzle exit (station 7)
            Pt7 == pifn * Pt21,     #B.125
            Tt7 == Tt21,    #B.126
            ht7 == ht21,    #B.127

            #LPC exit (station 2.5)
            Pt25 == pilc * Pt2,
            Tt25 == Tt2 ** pilc ** ((gammaAir-1)/(etalc*gammaAir)),
            ht25 == Tt25 * Cpair,

            #HPC Exit
            Pt3 == pihc * Pt25,
            Tt3 == Tt25 * pihc ** ((gammaAir-1)/(etalc*gammaAir)),
            ht3 == Cpair * Tt3
            ]

class CombustorCooling(Engine):
    """
    class to represent the engine's combustor and perform calculations
    on engine cooling bleed flow...cooling flow is currently not implemented
    """
    def __init__(self):
        #new vars
        #combustor exit state variables..recall Tt4 is already set in Engine class
        Pt4 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Combustor Exit (4)')
        ht4 = Variable('h_{t_4}', 'J', 'Stagnation Enthalpy at the Combustor Exit (4)')
        
        constraints = [
            #flow through combustor
            Pt4 == pib * Pt3,   #B.145
            ht4 == Cpc * Tt4,

            #fuel flow fraction f

            #flow at turbine inlet
            Tt41 == Tt4,
            Pt41 == Pt4,
            ht41 == ht4
            ]

class Turbine(Engine):
    """
    classs to represent the engine's turbine
    """
    def __init__(self):
        with SignomialsEnabled():
            #new variables
            ht45 = Variable('h_{t_4.5}', 'J', 'Stagnation Enthalpy at the HPT Exit (4.5)')

            #pressure ratios
            pihpt = Variable('\pi_{HPT}', '-', 'HPT Pressure Ratio')
            pilpt = Variable('\pi_{LPT}', '-', 'LPT Pressure Ratio')
            
            #turbine efficiences
            etaht = Variable('\eta_{ht}', '-', 'Polytropic Efficiency of HPT')
            
            etalt = Variable('\eta_{lt}', '-', 'Polytropic Efficiency of LPT')
            contraints = [
                #HPT shafter power balance

                #SIGNOMIAL
                (1+f)*(ht41-ht45) == ht3 - ht25,    #B.161

                #LPT shaft power balance

                #SIGNOMIAL  
                ht49 - ht45 == -((ht25-ht19)+alpha*(ht21 - ht2))/(1+f), #B.165

                #HPT Exit states (station 4.5)
                Pt45 == pihpt * Pt41,
                pihpt == (Tt45/Tt41)**(gammaT/(pihpt*(gammaT-1))),
                ht45 == Cpt * Tt45,

                #LPT Exit States
                Pt49 == pilpt * Pt45,
                pilpt == (Tt49/Tt45)**(gammaT/(pihpt*(gammaT-1))),
                ht49 == Cpt * Tt49,

                #turbine nozzle exit states
                Pt5 == pitn * Pt49, #B.167
                Tt5 == Tt49,    #B.168
                ht5 == ht49     #B.169
                ]

class ExhaustAndThrust(Engine):
    """
    Class to calculate the exhaust quantities as well as
    the overall engine thrust and on design TSFC & ISP
    """
    def __init__(self):
        #new variables
        #fan exhaust variables
        P8 = Variable('P_8', 'kPa', 'Fan Exhaust Static Pressure')
        Pt8 = Variable('P_{t_8}', 'kPa', 'Fan Exhaust Stagnation Pressure')
        ht8 = Variable('h_{t_8}', 'J', 'Fan Exhaust Stagnation Enthalpy')
        h8 = Variable('h_8', 'J', 'Fan Exhasut Static Enthalpy')
        
        #core exit variables
        P6 = Variable('P_6', 'kPa', 'Core Exhaust Static Pressure')
        Pt6 = Variable('P_{t_6}', 'kPa', 'Core Exhaust Stagnation Pressure')
        ht6 = Variable('h_{t_6}', 'J', 'Core Exhaust Stagnation Enthalpy')
        h6 = Variable('h_6', 'J', 'Core Exhasut Static Enthalpy')

        #thrust variables
        F8 = Variable('F_8', 'N', 'Fan Thrust')
        F6 = Variable('F_6', 'N', 'Core Thrust')
        F = Variable('F', 'N', 'Total Thrust')
        Fsp = Variable('F_{sp}', 'N/kg', 'Specific Net Thrust')

        with SignomialsEnabled():
            constraints = [
                #fan exhaust
                Pt8 == Pt7, #B.179
                Tt8 == Tt7, #B.180
                P8 == P0,
                h8 == Cpt * T8,
                u8**2 + 2*h8 <= 2*ht8,
                (P8/Pt8)**((gammaT-1)/gammaT) == T8/Tt8,
                
                #core exhaust
                P6 == P0,   #B.4.11 intro
                Pt6 == Pt5, #B.183
                Tt6 == Tt5, #B.184
                (P6/Pt6)**((gammaT-1)/gammaT) == T6/Tt6,
                u6**2 + 2*h6 <= 2*ht6,
                h6 == Cpt * T6,

                #overall thrust values
                F8/(alpha * mCore) + u0 <= u8,  #B.188
                F6/mCore + u0 <= (1+f)*u6,      #B.189

                #SIGNOMIAL
                F == F6 + F8,

                Fsp == F/((1+alpha)*mCore*a0),   #B.191

                #ISP
                Isp == Fsp*a0*(1+alpha)/(f*g),  #B.192

                #TSFC
                TSFC == 1/ISP                   #B.193
                ]

class OnDesignSizing(Engine):
    """
    class to perform the on design sizing of the engine.
    Nozzle areas are calcualted in post processing due to dependence on
    M6 and M8 being greater than or less than 1
    """
    def __init__(self):
        #new variables
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
        u25 = Variable('u_{2.5}', 'm/s', 'Air Speed at HPC Face')
        A25 = Variable('A_{2.5}', 'm^2', 'HPC Area')

        #mach numbers
        M8 = Variable('M_8', '-', 'Fan Exhaust Mach Number')
        M6 = Variable('M_6', '-', 'Core Exhaust Mach Number')
        
        constraints = [
            #mass flow sizing
            mCore == Fd/(Fsp*a0*(1+alpha)),  #B.194

            #component area sizing
            #fan area
            P2 == Pt2*(hold2)**(-gammaAir/(gammaAir-1)),
            T2 == Tt2 * hold2**-1,
            h2 == Cpair * T2,
            rho2 == P2/(Rair * T2),  #B.196
            u2 == M2*(Cpair*Rair*T2/(Cpair-Rair))**.5,  #B.197
            A2 == (1+alpha)*mCore/(rho2*u2),     #B.198

            #HPC area
            P25 == Pt25*(hold25)**(-gammaAir/(gammaAir-1)),
            T2 == Tt25 * hold25**-1,
            h25 == Cpair * T25,
            rho25 == P25/(Rair*T25),
            u25 == M25*(Cpair*Rair*T25/(Cpair-Rair))**.5,   #B.202
            A25 == (1+alpha)*mCore/(rho25*u25),     #B.203

            #mach nubmers for post processing of the data
            M8 == u8/((CpAir*Rair/(CpAir-Rair))**.5),
            M8 == u6/((Cpt*Rt/(Cpt-Rt))**.5),
            ]

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


