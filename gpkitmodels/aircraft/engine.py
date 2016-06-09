#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing

#TODO
#figure out how to compute f
#determine the three different Cp, gamma, and R values to be used
#implement the constraints on the nozzle sizing
#get realisitc R and Cp values for the different engine sections
#verify the values of all constants...and fix where differnces of constants are hard coded

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
        #set up the overeall model for an on design solve
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing()

        mCore = Variable('m_{core}', 'kg/s', 'Core Mass Flow')

        self.submodels = [lpc, combustor, turbine, thrust, size]

        #print [self.submodels, constraints]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

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
            'alphap1': 11,
            'M_{4a}': 1,    #choked turbines
            'F_D': 121436.45, #737 max thrust in N
            'M_2': .4,
            'M_{2.5}': .5
            }

            #temporary objective is to minimize the core mass flux
            objective = mCore
            Model.__init__(self, objective, lc, substitutions)

if __name__ == "__main__":
    engine = Engine()
    sol = engine.localsolve()
            

                #physical constants
        #Fix the R values
##        Rair = Variable('R_{air}', 287, 'J/kg/K', 'R for Air')
##        Rt = Variable('R_t', 287, 'J/kg/K', 'R for the Turbine Gas')
##        Rc = Variable('R_c', 287, 'J/kg/K', 'R for the Combustor Gas')
##        gammaAir = Variable('gamma_{air}', 1.4, '-', 'Specific Heat Ratio for Air')
##        gammaT = Variable('gamma_{T}', 1.4, '-', 'Specific Heat Ratio for Gas in Turbine')
##        gammaC =  Variable('gamma_{C}', 1.4, '-', 'Specific Heat Ratio for Gas in Combustor')
##        g = Variable('g', 9.81, 'm/s/s', 'Gravitational Acceleration')
##
##        #NEED TO FIX THESE VALUES
##        Cpair = Variable('Cp_{air}', 1068, 'J/kg/K', "Cp Value for Air at 673K") #http://www.engineeringtoolbox.com/air-properties-d_156.html
##        Cpc = Variable('Cp_c', 1068, 'J/kg/K', "Cp Value for Fuel/Air Mix in Combustor")
##        Cpt =Variable('Cp_t', 1068, 'J/kg/K', "Cp Value for Combustion Products in Turbine")
##
##        #variables linking subclasses
##        #enthalpies used in shaft power balances
##        ht18 = Variable('h_{t_1.8}', 'J', 'Stagnation Enthalpy at the Diffuser Exit (1.8)')
##        ht2 = Variable('h_{t_2}', 'J', 'Stagnation Enthalpy at the Fan Inlet (2)')
##        ht21 = Variable('h_{t_2.1}', 'J', 'Stagnation Enthalpy at the Fan Inlet (2.1)')
##        ht25 = Variable('h_{t_2.5}', 'J', 'Stagnation Enthalpy at the LPC Exit (2.5)')
##        
##        #HPC exit state variables (station 3)
##        Pt3 = Variable('P_{t_3}', 'kPa', 'Stagnation Pressure at the HPC Exit (3)')
##        Tt3 = Variable('T_{t_3}', 'K', 'Stagnation Temperature at the HPC Exit (3)')
##        ht3 = Variable('h_{t_3}', 'J', 'Stagnation Enthalpy at the HPC Exit (3)')
##
##        #Turbine inlet state variables (station 4.1)
##        Pt41 = Variable('P_{t_4}', 'kPa', 'Stagnation Pressure at the Turbine Inlet (4.1)')
##        Tt41 = Variable('T_{t_4}', 'J', 'Stagnation Temperature at the Turbine Inlet (4.1)')
##        ht41 = Variable('h_{t_4}', 'J', 'Stagnation Enthalpy at the Turbine Inlet (4.1)')
##
##        #HPC exit states
##        Pt45 = Variable('P_{t_4.5}', 'kPa', 'Stagnation Pressure at the HPT Exit (4.5)')
##        Tt45 = Variable('T_{t_4.5}', 'K', 'Stagnation Temperature at the HPT Exit (4.5)')
##
##        #turbine nozzle exit states
##        Pt5 = Variable('P_{t_5}', 'kPa', 'Stagnation Pressure at the Turbine Nozzle Exit (5)')
##        Tt5 = Variable('T_{t_5}', 'J', 'Stagnation Temperature at the Turbine Nozzle Exit (5)')
##        ht5 = Variable('h_{t_5}', 'J', 'Stagnation Enthalpy at the Turbine Nozzle Exit (5)')
##
##        #new vars for the fan nozzle exit (station 7)
##        Pt7 = Variable('P_{t_7}', 'kPa', 'Stagnation Pressure at the Fan Nozzle Exit (7)')
##        Tt7 = Variable('T_{t_7}', 'K', 'Stagnation Temperature at the Fan Nozzle Exit (7)')
##        ht7 = Variable('h_{t_7}', 'J', 'Stagnation Enthalpy at the Fan Nozzle Exit (7)')
##
##        #fuel
##        #flow faction f
##        f = Variable('f', '-', 'Fuel Air Mass Flow Fraction')
##
##        #exhaust speeds
##        u6 = Variable('u_6', 'm/s', 'Core Exhaust Velocity')
##        u8 = Variable('u_8', 'm/s', 'Fan Exhaust Velocity')
##        
##        #Variables taken as known, substitutded for later
##        T0 = Variable('T_0', 'K', 'Free Stream Stagnation Temperature')
##        P0 = Variable('P_0', 'kPa', 'Free Stream Stagnation Pressure')
##        M0 = Variable('M_0', '-', 'Free Stream Mach Number')
##        Tt4 = Variable('T_{t4}', 'K', 'Combustor Exit (Station 4) Stagnation Temperature')
##        pif = Variable('\pi_f', '-', 'Fan Pressure Ratio')
##        pilc = Variable('\pi_{lc}', '-', 'LPC Pressure Ratio')
##        pihc = Variable('\pi_{hc}', '-', 'HPC Pressure Ratio')
##        pitn = Variable('\pi_{tn}', '-', 'Turbine Nozzle Pressure Ratio')
##        pid = Variable('\pi_{d}', '-', 'Diffuser Pressure Ratio')
##        pib = Variable('\pi_{b}', '-', 'Burner Pressure Ratio')
##        pifn = Variable('\pi_{fn}', '-', 'Fan Duct Pressure Loss Ratio')
##        alpha = Variable('alpha', '-', 'By Pass Ratio')
##        M4a = Variable('M_{4a}', '-', 'Rep Mach # at Start of GPT Cooling Flow')
##
##
##        #known variables needed only for sizing purposes
##        Fd = Variable('F_D', 'N', 'Design Thrust')
##        M2 = Variable('M_2', '-', 'Fan Face/LPC Face Axial Mach Number')
##        M25 = Variable('M_{2.5}', '-', 'HPC Face Axial Mach Number')


