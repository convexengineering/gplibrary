#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing, OffDesign, FanMap, LPCMap, HPCMap

#TODO
#determine the three different Cp, gamma, and R values to be used
#get realisitc R and Cp values for the different engine sections
#verify the values of all constants...and fix where differnces of constants are hard coded
#replace the gammaAirs with gammaTs in post processing for nozzle at station 5

class EngineOnDesign(Model):
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
        m6opt = 1
        m8opt = 1
        
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing(m6opt, m8opt)

        self.submodels = [lpc, combustor, turbine, thrust, size]

        M2 = .6
        M25 = .6
        M4a = .6
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            'T_{t_4}': 1400,
            '\pi_f': 1.5,
            '\pi_{lc}': 4,
            '\pi_{hc}': 10,
            '\pi_{d}': 1,
            '\pi_{fn}': 1,
            '\pi_{tn}': 1,
            '\pi_{b}': .94,
            'alpha': 8,
            'alphap1': 9,
            'M_{4a}': M4a,    #choked turbines
            'F_D': 10.67*1000, #737 max thrust in N
            'M_2': M2,
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,

            #new subs for cooling flow losses
            'T_{t_f}': 1200,
            '\alpca_c': 0.1,
            'hold_{4a}': 1+.5*(1.313-1)*M25**2,
            'r_{uc}': 0.5
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, thrust.cost, lc, substitutions)

class EngineOffDesign(Model):
    """
    Engine Sizing for off design operation

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

    Off design model takes fan pressure ratio, LPC pressure ratio,
    HPC pressure ratio, fan corrected mass flow, LPC corrected mass flow,
    HPC corrected mass flow, Tt4, and Pt5 as uknowns that are solved for
    """
    def __init__(self, sol):
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 1
        m5opt = 0
        m7opt = 1
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': sol('T_0'),   #36K feet
                'P_0': sol('P_0'),    #36K feet
                'M_0': sol('M_0'),
                'M_{4a}': sol('M_{4a}'),
                '\pi_{tn}': sol('\pi_{tn}'),
                '\pi_{b}': sol('\pi_{b}'),
                '\pi_{d}': sol('\pi_{d}'),
                '\pi_{fn}': sol('\pi_{fn}'),
                
                'A_5': sol('A_5'),
                'A_7': sol('A_7'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': sol('m_{htD}'),
                'm_{ltD}': sol('m_{ltD}'),
                
                'G_f': 1,
                'alpha': 8,
                'alphap1': 9,
                
                'F_{spec}': 8.0e+04 ,
                'T_{t_{4spec}}': 1500,
                
                'm_{fan_D}': sol('alpha')*sol('m_{core}'),
                'N_{{bar}_Df}': 1,
                '\pi_{f_D}': sol('\pi_f'),
                'm_{core_D}': sol('m_{core}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
                'm_{lc_D}': sol('m_{lc_D}'),
                'm_{fan_bar_D}': sol('m_{fan_bar_D}'),
                'm_{hc_D}': sol('m_{hc_D}'),
                '\pi_{hc_D}': sol('\pi_{hc}')
            }
        
        Model.__init__(self, thrust.cost, lc, substitutions)
   
if __name__ == "__main__":
    engineOnD = EngineOnDesign()
    
    solOn = engineOnD.localsolve(verbosity = 4, kktsolver="ldl")
    
##    engineOffD = EngineOffDesign(solOn)
##    
##    solOff = engineOffD.localsolve(verbosity = 4, kktsolver="ldl",iteration_limit=200)
##    print solOff('F')
