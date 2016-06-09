from gpkit import Model, Variable, SignomialsEnabled
from gpkit.constraints.linked import LinkedConstraintSet
from gpkit.constraints.tight import TightConstraintSet as TCS
from engine_components import FanAndLPC, CombustorCooling, Turbine, ExhaustAndThrust, OnDesignSizing

if __name__ == "__main__":
    fan = FanAndLPC()
    fan.substitutions.update({
            'T_0': 216.5,   #36K feet
            'P_0': 22.8,    #36K feet
            'M_0': 0.8,
            '\pi_f': 1.5,
            '\pi_{lc}': 3,
            '\pi_{hc}': 10,
            '\pi_{d}': 1,
            '\pi_{fn}': 1
            })
    fansol = fan.solve(verbosity = 0)

    combustor = CombustorCooling()
    combustor.substitutions.update({
        'T_{t_4}': 1400,
        '\pi_{b}': 1,
        'h_{t_3}': fansol('h_{t_3}'),
        'P_{t_3}': fansol('P_{t_3}')
        })
    combsol = combustor.localsolve(verbosity = 0)
    print combsol('h_{t_4.1}')-(fansol('h_{t_3}')-fansol('h_{t_2.5}'))/(1+combsol('f'))
    print combsol('h_{t_4.1}')
    print (fansol('h_{t_2.5}')-fansol('h_{t_1.8}')+10*(fansol('h_{t_2.1}')-fansol('h_{t_2}')))
    
##    turbine = Turbine()
##    turbine.substitutions.update({
##        'alpha': 10,
##        'h_{t_3}': fansol('h_{t_3}'),
##        'h_{t_2.5}': fansol('h_{t_2.5}'),
##        'h_{t_1.8}': fansol('h_{t_1.8}'),
##        'h_{t_2.1}': fansol('h_{t_2.1}'),
##        'h_{t_2}': fansol('h_{t_2}'),
##        'h_{t_4.1}': combsol('h_{t_4.1}'),
##        'P_{t_4.1}': combsol('P_{t_4.1}'),
##        '\pi_{tn}': 1,
##        'h_{t_4.5}': 1093420.94025,
##        'f': combsol('f'),
##        'T_{t_4.1}': combsol('T_{t_4.1}')
##        })
##    turbinesol = turbine.solve()
##    #turbinesol = turbine.localsolve(kktsolver="ldl",verbosity = 1)

##    exhaust = ExhaustAndThrust()
##    exhaust.substitutions.update({
##        'a_0': fansol('a_0'),
##        'u_0': fansol('u_0'),
##        'T_0': 216.5,   #36K feet
##        'P_0': 22.8,    #36K feet
##        'M_0': 0.8,
##        'alpha': 10,
##        'alphap1': 11,
##        'M_{4a}': 1,    #choked turbines
##        'P_{t_7}': fansol('P_{t_7}'),
##        'T_{t_7}': fansol('T_{t_7}'),
##        'h_{t_7}': fansol('h_{t_7}'),
##        'P_{t_5}': 2026.5,
##        'T_{t_5}': 836.53,
##        'h_{t_5}': 893420.94025
##        })
##
##    exhaustsol = exhaust.solve()

    design = OnDesignSizing()
    design.substitutions.update({
        'a_0': fansol('a_0'),
        'u_0': fansol('u_0'),
        'T_0': 216.5,   #36K feet
        'P_0': 22.8,    #36K feet
        'M_0': 0.8,
        'alpha': 10,
        'alphap1': 11,
        'M_{4a}': 1,    #choked turbines
        'F_D': 121436.45, #737 max thrust in N
        'M_2': .4,
        'M_{2.5}': .5,
        'm_{core}': 40,
        'u_8': 600,
        'u_6': 625,
        'T_{t_2}': fansol('T_{t_2}'),
        'T_6': 700,
        'T_8': 500,
        'T_{t_2.5}': fansol('T_{t_2.5}')
        })

    designsol = design.solve()

    
