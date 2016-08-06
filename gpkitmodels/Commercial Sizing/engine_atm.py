#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units, ConstraintSet
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
        m8opt = 0
        
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing(m6opt, m8opt)

        self.submodels = [lpc, combustor, turbine, thrust, size]
            
        with SignomialsEnabled():

            substitutions = {
                'P_0': 19.8,  
                'M_0': 0.8,
                'T_{t_4}': 1400,
                '\pi_f': 1.5,
                '\pi_{lc}': 3,
                '\pi_{hc}': 10,
                'alpha': 8,
                'alphap1': 9,
                'M_{4a}': 1,    #choked turbines
                'M_2': .4,
                'M_{2.5}': .5,
                'hold_{2}': 1+.5*(1.398-1)*.4**2,
                'hold_{2.5}': 1+.5*(1.354-1)*.5**2,
                }
                    
            lc = LinkedConstraintSet([self.submodels])

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, None, lc, substitutions, **kwargs)
        

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
    def __init__(self, **kwargs):
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 1
        m5opt = 0
        m7opt = 0
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
            
        with SignomialsEnabled():

            constraints = ConstraintSet([self.submodels])

            constraints.subinplace({'M_0': 'M_0_1'})

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 230,   #36K feet
                'P_0': 19.8,  
##                'M_0': 0.8,
            }
        
        Model.__init__(self, None, lc, substitutions, **kwargs)

class EngineOffDesign2(Model):
    def __init__(self, **kwargs):
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 1
        m5opt = 0
        m7opt = 0
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
            
        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_E': 'TSFC_E2', 'F': 'F_2', 'M_0': 'M_0_2'})
    

        lc = LinkedConstraintSet(constraints)

        substitutions = {
            'T_0': 230,   #36K feet
            'P_0': 19.8,  
##            'M_0_2': 0.8,
        }
        
        Model.__init__(self, None, lc, substitutions, **kwargs)

class EngineOffDesign3(Model):
    def __init__(self, **kwargs):
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 1
        m5opt = 0
        m7opt = 0
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
            
        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_E': 'TSFC_E3', 'F': 'F_3', 'M_0': 'M_0_3'})
    

        lc = LinkedConstraintSet(constraints)

        substitutions = {
            'T_0': 230,   #36K feet
            'P_0': 19.8,  
            'M_0_3': 0.8,
        }
        
        Model.__init__(self, None, lc, substitutions, **kwargs)

class EngineOffDesign4(Model):
    def __init__(self, **kwargs):
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
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]

        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_E': 'TSFC_E4', 'F': 'F_4', 'M_0': 'M_0_4'})

        lc = LinkedConstraintSet(constraints)

        substitutions = {
            'T_0': 230,   #36K feet
            'P_0': 19.8,  
            'M_0_4': 0.8,
        }
        
        Model.__init__(self, None, lc, substitutions, **kwargs)

class EngineOffDesign5(Model):
    def __init__(self, **kwargs):
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 0
        m5opt = 0
        m7opt = 1
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_E': 'TSFC_E5', 'F_{spec}': 'F_{spec5}', 'M_0': 'M_0_5'})
    

        lc = LinkedConstraintSet(constraints)

        substitutions = {
            'T_0': 230,   #36K feet
            'P_0': 19.8,  
            'M_0_5': 0.8,
        }
        Model.__init__(self, None, lc, substitutions, **kwargs)

class EngineOffDesign6(Model):
     def __init__(self, **kwargs):
        lpc = FanAndLPC()
        combustor = CombustorCooling()
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap()
        lpcmap = LPCMap()
        hpcmap = HPCMap()

        res7 = 0
        m5opt = 0
        m7opt = 1
        
        offD = OffDesign(res7, m5opt, m7opt)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        constraints = ConstraintSet([self.submodels])

        constraints.subinplace({'TSFC_E': 'TSFC_E6', 'F_{spec}': 'F_{spec6}', 'M_0': 'M_0_6'})
    

        lc = LinkedConstraintSet(constraints)

        substitutions = {
            'T_0': 230,   #36K feet
            'P_0': 19.8,  
            'M_0_6': 0.8,
        }
        Model.__init__(self, None, lc, substitutions, **kwargs)
   
if __name__ == "__main__":
    engineOnD = EngineOnDesign()
    
    solOn = engineOnD.localsolve(verbosity = 4, kktsolver="ldl")
    
    engineOffD = EngineOffDesign(solOn)
    
    solOff = engineOffD.localsolve(verbosity = 4, kktsolver="ldl",iteration_limit=200)
    print solOff('F')
