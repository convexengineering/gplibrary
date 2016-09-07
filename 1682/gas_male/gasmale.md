# GAS-POWERED, MEDIUM-ALTITUDE, LONG-ENDURANCE AIRCRAFT

This paper presents the designs achieved in the 82 project. This also presents all of the models used in the design of this aircraft.

## Aerodynamic Model
\input{Aerodynamics.vars.generated.tex}
\input{Aerodynamics.cnstrs.generated.tex}

## Atmospheric Model
\input{Atmosphere.vars.generated.tex}
\input{Atmosphere.cnstrs.generated.tex}

## Breguet Endurance Model
\input{BreguetEndurance.vars.generated.tex}
\input{BreguetEndurance.cnstrs.generated.tex}

## Breguet Range Model
\input{BreguetRange.vars.generated.tex}
\input{BreguetRange.cnstrs.generated.tex}

## Climb Model
\input{Climb.vars.generated.tex}
\input{Climb.cnstrs.generated.tex}

## Cruise Model
\input{Cruise.vars.generated.tex}
\input{Cruise.cnstrs.generated.tex}

## Engine Model
\input{Engine.vars.generated.tex}
\input{Engine.cnstrs.generated.tex}

## Fuel Model
\input{Fuel.vars.generated.tex}
\input{Fuel.cnstrs.generated.tex}

## Fuselage Model
\input{Fuselage.vars.generated.tex}
\input{Fuselage.cnstrs.generated.tex}

## Loiter Model
\input{Loiter.vars.generated.tex}
\input{Loiter.cnstrs.generated.tex}

## Mission Profile Model
\input{Mission.vars.generated.tex}
\input{Mission.cnstrs.generated.tex}

## Steady Level Flight Model
\input{SteadyLevelFlight.vars.generated.tex}
\input{SteadyLevelFlight.cnstrs.generated.tex}

## Structural Model
\input{Structures.vars.generated.tex}
\input{Structures.cnstrs.generated.tex}

## Weight Model
\input{Weight.vars.generated.tex}
\input{Weight.cnstrs.generated.tex}

## Wind Model
\input{Wind.vars.generated.tex}
\input{Wind.cnstrs.generated.tex}

## Overall Model
\input{GasMALE.vars.generated.tex}
\input{GasMALE.cnstrs.generated.tex}


```python
#inPDF: skip
from gasmale import GasMALE
```

```python
#inPDF: replace with sol.generated.tex

if __name__ == "__main__":
    M = GasMALE()
    sol = M.solve("mosek")
    from gpkit.small_scripts import unitstr

    def gen_model_tex(model, modelname):
        with open('%s.vars.generated.tex' % modelname, 'w') as f:
            f.write("\\begin{longtable}{llll}\n \\toprule\n")
            f.write("\\toprule\n")
            f.write("Variables & Value & Units & Description \\\\ \n")
            f.write("\\midrule\n")
            #f.write("\\multicolumn{3}{l}\n")
            varnames = ["firstname"]
            for var in model.varkeys:
                name = var.name
                if name in varnames:
                    pass
                else:
                    varnames.append(name)
                    unitstr = var.unitstr()[1:]
                    unitstr = "$[%s]$" % unitstr if unitstr else ""
                    val = "%0.3f" % var.value if var.value else ""
                    f.write("$%s$ & %s & %s & %s \\\\\n" % 
                                (var.name, val, unitstr, var.label))
            f.write("\\bottomrule\n")
            f.write("\\end{longtable}\n")

        with open('%s.cnstrs.generated.tex' % modelname, 'w') as f:
            lines = model.latex(excluded=["models"]).replace("[ll]", "{ll}")..split("\n")
            modeltex = "\n".join(lines[:1] + lines[3:])
            f.write("$$ %s $$" % modeltex)

    def find_base_submodel(models, modelnames):
        runAgain = 0
        for m in models:
            if "submodels" in m.__dict__.keys():
                for sub in m.submodels:
                    if sub.__class__.__name__ not in modelnames:
                        models.append(sub)
                        modelnames.append(sub.__class__.__name__)
                        runAgain += 1
                    else:
                        pass
            else:
                pass
        if runAgain > 0:
            return find_base_submodel(models, modelnames)
        else:
            return models, modelnames

    models, modelnames = find_base_submodel([M], [])
    for m in models: 
        gen_model_tex(m, m.__class__.__name__)

    with open("sol.generated.tex", "w") as f:
        f.write(sol.table(latex=True))

```




