import numpy as np
from gpkit import Model, Variable, SignomialsEnabled, units, ConstraintSet
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components_NPSS_GE90_CFM_noOnD_validation import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, Sizing, FanMap, LPCMap, HPCMap
from collections import defaultdict
from gpkit.small_scripts import mag

class SizingTOC(Model):

    def __init__(self):
        mixing = True
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        res7 = 0

        #need to give a Tt4 to run with res7 = 0

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 218,   #36K feet
                'P_0': 23.84,    #36K feet
                'M_0': .8,
                'M_2': M2,
                'M_{2.5}':M25,
                'hold_{2}': 1+.5*(1.398-1)*M2**2,
                'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            }
            
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1450,
                })
            else:
                substitutions.update({
                    'F_{spec}': 19600.4*4.4,
##                    'T_{t_4}': 1450,
                })
            
            
        Model.__init__(self, offD.cost, lc, substitutions)
        
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

        
class SizingOnDRerun(Model):

    def __init__(self):
        mixing = True
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        res7 = 0

        #need to give a Tt4 to run with res7 = 0

        M2 = .8
        M25 = .65
        M4a = .1025
        Mexit = 1
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 218,   #36K feet
                'P_0': 23.84,    #36K feet
                'M_0': .8,
                'M_2': M2,
                'M_{2.5}':M25,
                'hold_{2}': 1+.5*(1.398-1)*M2**2,
                'hold_{2.5}': 1+.5*(1.354-1)*M25**2,
            }
            
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1450,
                })
            else:
                substitutions.update({
                    'F_{spec}': 16408.4 * 4.4,
##                    'T_{t_4}': 1450,
                })
        Model.__init__(self, offD.cost, lc, substitutions)

class SizingTO(Model):

    def __init__(self):
        mixing = True
        SPmaps = False
        
        lpc = FanAndLPC()
        combustor = CombustorCooling(mixing)
        turbine = Turbine()
        thrust = ExhaustAndThrust()
        fanmap = FanMap(SPmaps)
        lpcmap = LPCMap(SPmaps)
        hpcmap = HPCMap(SPmaps)

        res7 = 0

        M2 = .25
        M25 = .15
        M4a = .1025
        Mexit = 1
        
        offD = Sizing(res7, mixing)

        #only add the HPCmap if residual 7 specifies a thrust
        if res7 ==0:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap, hpcmap]
        if res7 == 1 and SPmaps == True:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
        if res7 == 1 and SPmaps == False:
            self.submodels = [lpc, combustor, turbine, thrust, offD, fanmap, lpcmap]
            
        with SignomialsEnabled():

            lc = LinkedConstraintSet([self.submodels])

            substitutions = {
                'T_0': 315,   #36K feet
                'P_0': 101.325,    #36K feet
                'M_0': .25,
                'M_2': M2,
                'M_{2.5}':M25,

            }
            
            if res7 == 1:
               substitutions.update({
                    'T_{t_{4spec}}': 1400,
                })
            else:
                substitutions.update({
                    'F_{spec}': 78399.1*4.4,
                })
            
        Model.__init__(self, offD.cost, lc, substitutions)

if __name__ == "__main__":
    W_engine = Variable('W_{engine}', 'N', 'Weight of a Single Turbofan Engine')
    
    engine1 = SizingTOC()
    engine2 = SizingOnDRerun()
    engine3 = SizingTO()
 
##    sol1 = engine1.localsolve(verbosity = 4, solver="mosek")
##    bounds, sol = engine1.determine_unbounded_variables(engine1, solver="mosek",verbosity=4, iteration_limit=50)

##    sol2 = engine2.localsolve(verbosity = 4, solver="mosek")

    #create the big linked engine model
    submodels = [engine1, engine2]
    constraints = ConstraintSet([submodels])

    lc = LinkedConstraintSet(constraints, include_only = {'A_5', 'A_7', 'A_2', 'A_{2.5}', '\pi_{tn}', '\pi_{b}', '\pi_{d}', '\pi_{fn}',
                                                          'T_{ref}', 'P_{ref}', '\eta_{HPshaft}', '\eta_{LPshaft}', 'eta_{B}',
                                                          'W_{engine}', 'm_{fan_bar_D}', 'm_{lc_D}', 'm_{hc_D}', '\pi_{f_D}',
                                                          '\pi_{hc_D}', '\pi_{lc_D}', 'm_{htD}', 'm_{ltD}', 'm_{coreD}', 'M_{4a}',
                                                          'hold_{4a}', 'r_{uc}', '\\alpha_c', 'T_{t_f}', 'M_{takeoff}', 'G_f'})

    M4a = .10

    fan = 1.58
    lpc  = 1.26
    hpc = 20.033
    
    valsubs = {
##    'A_5': .2727,
##    'A_7': 1.1,
    '\pi_{tn}': .98,
    '\pi_{b}': .94,
    '\pi_{d}': .98,
    '\pi_{fn}': .98,
    'T_{ref}': 288.15,
    'P_{ref}': 101.325,
    '\eta_{HPshaft}': .97,
    '\eta_{LPshaft}': .97,
    'eta_{B}': .9970,

    '\pi_{f_D}': 1.58,
    '\pi_{hc_D}': 20.033,
    '\pi_{lc_D}': 1.26,

    '\\alpha_SizingOnDRerun': 8.7877,

    'M_{4a}': M4a,
    'hold_{4a}': 1+.5*(1.313-1)*M4a**2,#sol('hold_{4a}'),
    'r_{uc}': .01,
    '\\alpha_c': .14,
    'T_{t_f}': 435,

    'M_{takeoff}': .955,

    'G_f': 1,
    }

    m = Model((2*engine2.cost+engine1.cost), constraints, valsubs)
    
    Pt0 = 40
    Tt0 = 220
    Pt3 = Pt0*lpc*fan*hpc
    Pt21 = fan * Pt0
    Pt25 = Pt0 * fan * lpc

    Tt21 = Tt0 * (fan)**(.4/(1.4*.9005))
    Tt25 = Tt21 * (lpc)**(.4/(1.4*.9306))
    Tt3 = Tt25 * (hpc)**(.4/(1.4*.9030))

    Tt41 = 1300

    Tt45 = Tt41 - (Tt3 - Tt25)

    Tt49 = Tt45 - (Tt25 - Tt21)

    piHPT = (Tt45/Tt41)**(.8851*1.4/.4)

    piLPT = (Tt49/Tt45)**(.8731*1.4/.4)

    Pt45 = piHPT * Pt3

    print Tt21, Tt25, Pt21, Pt25, Tt41, Tt45, Pt3, Pt45
    
##    sol = m.localsolve(verbosity = 4, solver="mosek", iteration_limit=100)
    bounds, sol = engine1.determine_unbounded_variables(m, solver="mosek",verbosity=4, iteration_limit=200)
