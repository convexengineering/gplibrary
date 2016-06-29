#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing
from engine_components_test import FanAndLPCTEST, CombustorCoolingTEST, TurbineTEST, ExhaustAndThrustTEST, OffDesignTEST, FanMapTEST, LPCMapTEST
from collections import defaultdict
from gpkit.small_scripts import mag
#from gpkit.tools import determine_unbounded_variables

#TODO
#determine the three different Cp, gamma, and R values to be used
#get realisitc R and Cp values for the different engine sections
#verify the values of all constants...and fix where differnces of constants are hard coded
#figure out Tref and Pref and make them the right value in the sizing post processing
#replace the gammaAirs with gammaTs in post processing for nozzle at station 5
#verify that the Tref nad Pref are correct

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
            'hold_{2.5}': 1.05,
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, thrust.cost, lc, substitutions)
            

    def sizing(self, sol):
        #find the design normalized speeds and mass flow
        mhtD = (1+sol('f'))*sol('m_{core}')*((sol('T_{t_4.1}')/(288.15*units('K')))**.5)/(sol('P_{t_4.1}')/(101.325*units('kPa'))) #B.225
        mltD = (1+sol('f'))*sol('m_{core}')*((sol('T_{t_4.5}')/(288.15*units('K')))**.5)/(sol('P_{t_4.5}')/(101.325*units('kPa'))) #B.226
        #normalized on design compressor mass flows for the maps
        mFanBarD = sol('alpha')*sol('m_{core}')*((sol('T_{t_2}')/(288.15*units('K')))**.5)/(sol('P_{t_2}')/(101.325*units('kPa'))) #B.226
        mlcD = sol('m_{core}')*((sol('T_{t_2}')/(288.15*units('K')))**.5)/(sol('P_{t_2}')/(101.325*units('kPa'))) #B.226
        #noting that 
        NlpcD = 1/(sol('T_{t_1.8}')/(288.15*units('K')))**.5    #B.223
        NhpcD = 1/(sol('T_{t_2.5}')/(288.15*units('K')))**.5    #B.224

        #first size the fan nozzle
        if sol('M_8') >= 1:
            M7 = 1
            P7 = sol('P_{t_7}')*(1+.5*(sol('gamma_{air}')-1))**(-sol('gamma_{air}')/(sol('gamma_{air}')-1))
            T7 = sol('T_{t_7}')*(1+.5*(sol('gamma_{air}')-1))**-1
        else:
            P7 = sol('P_0')
            T7 = sol('T_{t_7}')*(P5/sol('P_{t_7}'))**((sol('gamma_{air}')-1)/sol('gamma_{air}'))
            
        h7 = sol('Cp_{air}')*T7    
        #compute the fan stream velocity and nozzle area
        u7 = (2*(sol('h_{t_7}')-h7)**.5)
        rho7 = P7/(sol('R_t')*T7)
        A7 = sol('alpha')*sol('m_{core}')/(rho7*u7)
        
        #core nozzle area calcualtion
        if sol('M_6') >=1:
            M5 = 1
            P5 = sol('P_{t_5}')*(1+.5*(sol('gamma_{air}')-1))**(-sol('gamma_{air}')/(sol('gamma_{air}')-1))
            T5 = sol('T_{t_5}')*(1+.5*(sol('gamma_{air}')-1))**-1  
        else:
            P5 = sol('P_0')
            T5 = sol('T_{t_5}')*(P5/sol('P_{t_5}'))**((sol('gamma_{air}')-1)/sol('gamma_{air}'))
            
        h5 = sol('Cp_{air}')*T5
        #compute the core stream velocity and nozzle area    
        u5 = (2*(sol('h_{t_5}')-h5)**.5)
        rho5 = P5/(sol('R_t')*T5)
        A5 = sol('m_{core}')/(rho5*u5)
       
        return mhtD, mltD, mFanBarD, mlcD, NlpcD, NhpcD, A5, A7

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
    def bound_all_variables(self, model, eps=1e-30, lower=None, upper=None):
        "Returns model with additional constraints bounding all free variables"
        lb = lower if lower else eps
        ub = upper if upper else 1/eps
        constraints = []
        freevks = tuple(vk for vk in model.varkeys if "value" not in vk.descr)
        for varkey in freevks:
            units = varkey.descr.get("units", 1)
            constraints.append([ub*units >= Variable(**varkey.descr),
                                Variable(**varkey.descr) >= lb*units])
        m = Model(model.cost, [constraints, model], model.substitutions)
        m.bound_all = {"lb": lb, "ub": ub, "varkeys": freevks}
        return m


    # pylint: disable=too-many-locals
    def determine_unbounded_variables(self, model, solver=None, verbosity=0,
                                      eps=1e-30, lower=None, upper=None, **kwargs):
        "Returns labeled dictionary of unbounded variables."
        m = self.bound_all_variables(model, eps, lower, upper)
        sol = m.localsolve(solver, verbosity, **kwargs)
        lam = sol["sensitivities"]["la"][1:]
        out = defaultdict(list)
        for i, varkey in enumerate(m.bound_all["varkeys"]):
            lam_gt, lam_lt = lam[2*i], lam[2*i+1]
            if abs(lam_gt) >= 1e-7:  # arbitrary threshold
                out["sensitive to upper bound"].append(varkey)
            if abs(lam_lt) >= 1e-7:  # arbitrary threshold
                out["sensitive to lower bound"].append(varkey)
            value = mag(sol["variables"][varkey])
            distance_below = np.log(value/m.bound_all["lb"])
            distance_above = np.log(m.bound_all["ub"]/value)
            if distance_below <= 3:  # arbitrary threshold
                out["value near lower bound"].append(varkey)
            elif distance_above <= 3:  # arbitrary threshold
                out["value near upper bound"].append(varkey)
        return out

    def __init__(self, sol, mhtD, mltD, mFanBarD, mlcD, NlpcD, NhpcD, A5, A7, **kwargs):
        lpc = FanAndLPCTEST()
        combustor = CombustorCoolingTEST()
        turbine = TurbineTEST()
        thrust = ExhaustAndThrustTEST()
        fanmap = FanMapTEST()
        lpcmap = LPCMapTEST()

        res7 = 1
        m5opt = 1
        m7opt = 1
        
        offD = OffDesignTEST(res7, m5opt, m7opt)

        self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': sol('T_0'),   #36K feet
                'P_0': sol('P_0'),    #36K feet
                'M_0': sol('M_0'),
                '\pi_{tn}': sol('\pi_{tn}'),
                '\pi_{b}': sol('\pi_{b}'),
                '\pi_{d}': sol('\pi_{d}'),
                '\pi_{fn}': sol('\pi_{fn}'),
                
                'A_5': A5.to('m^2'),
                'A_7': A7.to('m^2'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': mhtD,
                'm_{ltD}': mltD,
                
                #'N_1': 1,
                #'\pi_{lc}': sol('\pi_{lc}'),
                #'\pi_{hc}': sol('\pi_{hc}'),
                #'\pi_f': sol('\pi_f'),
                #'m_{core}': sol('m_{core}'),
                #might need Tt4 in full reversion
                
                'G_f': 1,
                'alpha': 10,
                'alphap1': 11,
                
                'F_{spec}': sol('F_D') ,
                'T_{t_{4spec}}': 1400,
                
                'm_{fan_D}': sol('alpha')*sol('m_{core}'),
                'N_{{bar}_Df}': 1,
                '\pi_{f_D}': sol('\pi_f'),
                'm_{core_D}': sol('m_{core}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
                'N_{{bar}_D_lc}': NlpcD,
                'm_{lc_D}': mlcD,
                'm_{fan_bar_D}': mFanBardD,
            }

        Model.__init__(self, thrust.cost, lc, substitutions)

 

            
if __name__ == "__main__":
    engineOnD = EngineOnDesign()
    
    solOn = engineOnD.localsolve(verbosity = 1, kktsolver="ldl")
    
    mhtD, mltD, mFanBardD, mlcD, NlpcD, NhpcD, A5, A7 = engineOnD.sizing(solOn)
    
    engineOffD = EngineOffDesign(solOn, mhtD, mltD, mFanBardD, mlcD, NlpcD, NhpcD, A5, A7)
    
    solOff = engineOffD.localsolve(verbosity = 4, kktsolver="ldl")
##    bounds = engineOffD.determine_unbounded_variables(engineOffD,verbosity=4)
    
    #print model.program.solver_out
    #print solOff('\rho_5')
