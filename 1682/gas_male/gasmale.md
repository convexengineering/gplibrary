# GAS-POWERED, MEDIUM-ALTITUDE, LONG-ENDURANCE AIRCRAFT

This paper presents the designs achieved in the 82 project. 

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
            for var in model.varkeys:
                unitstr = var.unitstr()[1:]
                unitstr = "$[%s]$" % unitstr if unitstr else ""
                val = "%0.3f" % var.value if var.value else ""
                f.write("$%s$ & %s & %s & %s \\\\\n" % (var.name, val, unitstr, var.label))
            f.write("\\bottomrule\n")
            f.write("\\end{longtable}\n")

        with open('%s.cnstrs.generated.tex' % modelname, 'w') as f:
            lines = model.latex(excluded=["models"]).replace("[ll]", "{ll}").split("\n")
            modeltex = "\n".join(lines[:1] + lines[3:])
            f.write("$$ %s $$" % modeltex)

    for model in M.submodels:
        gen_model_tex(model, model.__class__.__name__)
    with open("sol.generated.tex", "w") as f:
        f.write(sol.table(latex=True))

```




