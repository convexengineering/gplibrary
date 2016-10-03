# Gas MALE using 100 W payload

This paper evaluates using a payload that draws 100 W instead of 10W. 

## Sizing

This model was created and then a sweep was done to determine the MTOW required to meet 5 days. 

Using a 100W payload the trend of MTOW vs loiter time changes.  The graph below represents the change and the new design space.
```python
#inPDF: replace with tex/tstation_vs_MTOW_rubber100W.fig.generated.tex
from gasmale import GasMALE
from plotting import plot_sweep, fix_vars, plot_altitude_sweeps
from gen_tex import gen_tex_fig, gen_fixvars_tex
import numpy as np

M = GasMALE()
M.substitutions.update({"MTOW": 150})
M.substitutions.update({"P_{pay}": 100})
fig, ax = plot_sweep(M, "MTOW", np.linspace(70, 500, 15), ["t_{loiter}"])
gen_tex_fig(fig, "tstation_vs_MTOW_rubber100W")
```

## With DF70 engine

Using the DF70 engine, the following solution is the best possible, given the 100W payload. Note, I am maximizing time on station. 

\input{tex/maxt_loiter100W.generated.tex}

In order to keep under 150 lbs we have to sacrifice endurance to 5.6 days. The solution to this aircraft is found below:  
```python
#inPDF: replace with tex/sol100W.generated.tex
M = GasMALE(DF70=True)
M.substitutions.update({"P_{pay}": 100})
sol = M.solve("mosek")

with open("tex/maxt_loiter100W.generated.tex", "w") as f:
    f.write("The maximum time on station possible with the DF70 engine is %0.3f days" % sol("t_{loiter}").magnitude)
    f.write("This is achieved at a max take of weight of %0.3f lbs and a wing span of %0.3f feet" % (sol("MTOW").magnitude, sol("b").magnitude))

M.substitutions.update({"t_{loiter}": 5.6})
M.cost = M["MTOW"]
with open("tex/sol100W.generated.tex", "w") as f:
    f.write(sol.table(latex=True))
```

By fixing the following variables to their respective values we were also able to generate performance curves. 

```python
#inPDF: replace with tex/fixvars.table.generated.tex

vars_to_fix = {"S":0.0, "b": 0.0, "Vol_{fuse}":0.0, "W_Wing, GasMALE":0.0}
gen_fixvars_tex(M, sol, vars_to_fix)

fix_vars(M, sol, vars_to_fix)
sol = M.solve("mosek") # check for solving errors
```

## Sweeps

```python
#inPDF: skip
del M.substitutions["t_{loiter}"]
M.cost = 1/M["t_{loiter}"]

# payload weight vs time on station
fig, ax = plot_sweep(M, "W_{pay}", np.linspace(5, 20, 15), ["t_{loiter}"], ylim=[0, 10])
gen_tex_fig(fig, "t_vs_Wpay100W")

# wind speed vs time on station
M = GasMALE(wind=True, DF70=True)
fix_vars(M, sol, vars_to_fix)
fig, ax = plot_sweep(M, "V_{wind}_Wind, Loiter, Mission, GasMALE", np.linspace(5, 40, 15), ["t_{loiter}"], ylim=[0,10])
gen_tex_fig(fig, "t_vs_Vwind100W")

# altitude vs time on loiter
altitude_vars = {"t_{loiter}"}
figs, axs = plot_altitude_sweeps(np.linspace(14000, 20000, 10), altitude_vars,
                     vars_to_fix)
axs[0].set_ylim([0,10])
for var, fig in zip(altitude_vars, figs):
    gen_tex_fig(fig, "%s_vs_altitude100W" % var.replace("{", "").replace("}", "").replace("_", ""))
```
\input{tex/t_vs_Wpay100W.fig.generated.tex}
\input{tex/t_vs_Vwind100W.fig.generated.tex}
\input{tex/tloiter_vs_altitude100W.fig.generated.tex}

## Flight Profile

By further discretizing the climb, cruise, and loiter mission segments the following figures were generated to follow the performance over the duration of the mission. 

```python
#inPDF: skip
from plotting import plot_mission_var

Mprof = GasMALE(DF70=True)
Mprof.substitutions.update({"P_{pay}": 100})
M.substitutions.update({"t_{loiter}": 5.6})
M.cost = M["MTOW"]
sol = Mprof.solve("mosek")

# plot mission profiles
fig, ax = plot_mission_var(Mprof, sol, "V", [0, 40])
gen_tex_fig(fig, "profile_velocity100W")

fig, ax = plot_mission_var(Mprof, sol, "h", [0, 20000])
gen_tex_fig(fig, "profile_altitude100W")

fig, ax = plot_mission_var(Mprof, sol, "\\eta_{prop}", [0, 1])
gen_tex_fig(fig, "profile_etaprop100W")

fig, ax = plot_mission_var(Mprof, sol, "BSFC", [0, 2])
gen_tex_fig(fig, "profile_BSFC100W")

fig, ax = plot_mission_var(Mprof, sol, "P_{shaft-max}", [0, 5])
gen_tex_fig(fig, "profile_Pshaftmax100W")

fig, ax = plot_mission_var(Mprof, sol, "P_{shaft-tot}", [0, 5])
gen_tex_fig(fig, "profile_Pshafttot100W")

fig, ax = plot_mission_var(Mprof, sol, "RPM", [0, 9000])
gen_tex_fig(fig, "profile_RPM100W")

fig, ax = plot_mission_var(Mprof, sol, "W_{N+1}", [0, 150], "aircraft weight [lbf]")
gen_tex_fig(fig, "profile_weight100W")
```
\input{tex/profile_velocity100W.fig.generated.tex}
\input{tex/profile_etaprop100W.fig.generated.tex}
\input{tex/profile_BSFC100W.fig.generated.tex}
\input{tex/profile_Pshaftmax100W.fig.generated.tex}
\input{tex/profile_Pshafttot100W.fig.generated.tex}
\input{tex/profile_RPM100W.fig.generated.tex}
\input{tex/profile_weight100W.fig.generated.tex}
