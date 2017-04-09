#Implements the TASOPT engine model, currently disregards BLI
import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components_works import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing, OffDesign, FanMap, LPCMap, HPCMap
from collections import defaultdict
from gpkit.small_scripts import mag
#TODO
#get jet fuel Cp

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
        mixing = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        size = OnDesignSizing()

        self.submodels = [lpc, combustor, turbine, thrust, size]

        M2 = .6
        M25 = .6
        M4a = .6
        Mexit = 1
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            'T_{t_4}': 1400,
            '\pi_f': 1.5,
            '\pi_{lc}': 3.28,
            '\pi_{hc}': 10,
            '\pi_{d}': .98,
            '\pi_{fn}': .98,
            '\pi_{tn}': .98,
            '\pi_{b}': .94,
            'alpha': 5.5,
            'alphap1': 6.5,
            'M_{4a}': M4a,    #choked turbines
            'F_D': 86.7*1000, #737 max thrust in N
            'M_2': M2,
            'M_{2.5}':M25,
            'hold_{2}': 1+.5*(1.398-1)*M2**2,
            'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            'T_{ref}': 288.15,
            'P_{ref}': 101.325,
            '\eta_{HPshaft}': .99,
            '\eta_{LPshaft}': .98,
            'M_{takeoff}': .95,
            'stag41': 1+.5*(.312)*M4a**2,
            '\alpca_c': .05,

            #new subs for mixing flow losses
            'T_{t_f}': 600,
            'hold_{4a}': 1+.5*(1.313-1)*M4a**2,
            'r_{uc}': 0.5,
            'T_{m_TO}': 1000,
            'M_{t_exit}': Mexit,
            'chold_2': (1+.5*(1.318-1)*Mexit**2)**-1,
            'chold_3': (1+.5*(1.318-1)*Mexit**2)**-2,
            'T_{t_4TO}': 1600,
            }

            #temporary objective is to minimize the core mass flux 
            Model.__init__(self, size.cost, lc, substitutions)

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
            solhold = sol
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
            return out, solhold

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
        mixing = False
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        res7 = 1
        
        offD = OffDesign(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        else:
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
                
                'A_5': sol('A_5'),
                'A_7': sol('A_7'),
                'T_{ref}': 288.15,
                'P_{ref}': 101.325,
                'm_{htD}': sol('m_{htD}'),
                'm_{ltD}': sol('m_{ltD}'),
                
                'G_f': 1,
                'alpha': sol('alpha'),
                'alphap1': sol('alphap1'),
                
                'F_{spec}': 8.0e+04 ,
                'T_{t_{4spec}}': 1200,

                '\eta_{HPshaft}': sol('\eta_{HPshaft}'),
                '\eta_{LPshaft}': sol('\eta_{LPshaft}'),
                'M_{takeoff}': sol('M_{takeoff}'),
                '\alpca_c': sol('\alpca_c'),
                'T_{t_f}': sol('T_{t_f}'),
                
                'm_{fan_D}': sol('alpha')*sol('m_{core}'),
                'N_{{bar}_Df}': 1,
                '\pi_{f_D}': sol('\pi_f'),
                'm_{core_D}': sol('m_{core}'),
                '\pi_{lc_D}': sol('\pi_{lc}'),
                'm_{lc_D}': sol('m_{lc_D}'),
                'm_{fan_bar_D}': sol('m_{fan_bar_D}'),
                'm_{hc_D}': sol('m_{hc_D}'),
                '\pi_{hc_D}': sol('\pi_{hc}'),

            }
            if mixing == True:
                substitutions.update({
                    'stag41': 1+.5*(.312)*sol('M_{4a}')**2,
                    'M_{4a}': sol('M_{4a}'),
                    'hold_{4a}': 1+.5*(1.313-1)*.6**2,#sol('hold_{4a}'),
                    'r_{uc}': sol('r_{uc}'),
                })
            
        Model.__init__(self, thrust.cost, lc, substitutions)
   
if __name__ == "__main__":
    engineOnD = EngineOnDesign()
    
    solOn = engineOnD.localsolve(verbosity = 4, solver="mosek")
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOnD, solver="mosek",verbosity=4, iteration_limit=100)
    
##    engineOffD = EngineOffDesign(solOn)
    
##    solOff = engineOffD.localsolve(verbosity = 4, solver="mosek",iteration_limit=100)
##    bounds, sol = engineOnD.determine_unbounded_variables(engineOffD, solver="mosek",verbosity=4, iteration_limit=100)
