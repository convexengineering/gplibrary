import numpy as np
import numpy as np
from gasmale import GasMALEDiscritizedClimb
import plot
from plotting import fix_vars, plot_mission_var, plot_altitude_sweeps

fixed = False
# PLOTS #
fixedPLOTS = False
missionPLOTS = False
wind = False
W_pay = False
P_pay = False
altitude = False

M = GasMALEDiscritizedClimb()
sol_fix = M.solve('mosek')

vars_to_fix = {"S":0.0, "b":0.0, "Vol_{fuse}":0.0001}
fix_vars(M, sol_fix, vars_to_fix)

# plot mission profiles
fig, ax = plot_mission_var(M, "V", [0, 40])
fig.savefig("profile_t_vs_velocity.pdf")

fig, ax = plot_mission_var(M, "\\eta_{prop}", [0, 1])
fig.savefig("profile_t_vs_etaprop.pdf")

fig, ax = plot_mission_var(M, "BSFC", [0, 2])
fig.savefig("profile_t_vs_BSFC.pdf")

fig, ax = plot_mission_var(M, "P_{shaft-max}", [0, 5])
fig.savefig("profile_t_vs_Pshaftmax.pdf")

fig, ax = plot_mission_var(M, "P_{shaft-tot}", [0, 5])
fig.savefig("profile_t_vs_Pshafttot.pdf")

fig, ax = plot_mission_var(M, "RPM", [0, 9000])
fig.savefig("profile_t_vs_RPM.pdf")

fig, ax = plot_mission_var(M, "W_{end}", [0, 150], "aircraft weight [lbf]")
fig.savefig("profile_t_vs_weight.pdf")

# plot altitude sweeps
plot_altitude_sweeps(np.linspace(14000, 24000, 20), {"t_{station}"},
                     vars_to_fix, CLIMB=True)




#numplots = sum([wind, P_pay, W_pay, altitude])
#if numplots > 1:
#    raise ValueError("more than one is true, but there can only be one!")
#elif numplots == 1:
#    if "t_{station}" in M.substitutions:
#        del M.substitutions["t_{station}"]
#    M.cost = 1/M["t_{station}"]
#    if wind:
#        plot.wind(M)
#    elif P_pay:
#        plot.P_pay(M)
#    elif W_pay:
#        plot.W_pay(M)
#    elif altitude:
#        plot.altitude(M)
