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

        self.submodels = [lpc, combustor, turbine, thrust, size]
            
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
            '\pi_{d}': 1,
            '\pi_{fn}': 1,
            '\pi_{tn}': 1,
            '\pi_{b}': 1,
            'alpha': 10,
            'alphap1': 11,
            'M_{4a}': 1,    #choked turbines
            'F_D': 121436.45, #737 max thrust in N
            'M_2': .4,
            'M_{2.5}': .5,
            'hold_{2}': 1.032,
            'hold_{2.5}': 1.05
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, thrust.cost, lc, substitutions)
            

    def sizing(self, sol):
        #first size the fan nozzle
        if sol('M_8') >= 1:
           M7 = (sol('P_7')/sol('P_{t_7}')**((sol('gamma_{air}')-1)/(-sol('gamma_{air}')))-1)*2/(sol('gamma_{air}')-1)
        #size the core nozzle

        #return the two sizes
        print actual

if __name__ == "__main__":
    engine = Engine()
    sol = engine.localsolve(verbosity = 4, kktsolver="ldl")
    engine.sizing(sol)
