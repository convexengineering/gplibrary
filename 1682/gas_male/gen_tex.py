from gpkit.small_scripts import unitstr

def gen_model_tex(model, modelname, texname=None):
    if texname:
        filename = texname
    else:
        filename = modelname
    with open('%s.vars.generated.tex' % filename, 'w') as f:
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
                if var.models[0] == modelname:
                    varnames.append(name)
                    unitstring = var.unitstr()[1:]
                    unitstring = "$[%s]$" % unitstring if unitstring else ""
                    val = "%0.3f" % var.value if var.value else ""
                    f.write("$%s$ & %s & %s & %s \\\\\n" %
                            (var.name, val, unitstring, var.label))
                else:
                    pass
        f.write("\\bottomrule\n")
        f.write("\\end{longtable}\n")

    with open('%s.cnstrs.generated.tex' % texname, 'w') as f:
        lines = model.latex(excluded=["models"]). \
                replace("[ll]", "{ll}").split("\n")
        modeltex = "\n".join(lines[:1] + lines[3:])
        f.write("$$ %s $$" % modeltex)

def find_submodels(models, modelnames):
    runagain = 0
    for m in models:
        if "submodels" in m.__dict__.keys():
            for sub in m.submodels:
                if sub.__class__.__name__ not in modelnames:
                    models.append(sub)
                    modelnames.append(sub.__class__.__name__)
                    runagain += 1
                else:
                    pass
        else:
            pass
    if runagain > 0:
        return find_submodels(models, modelnames)
    else:
        return models, modelnames

def gen_tex_fig(fig, filename, caption=None):
    fig.savefig("%s.pdf" % filename)
    with open("%s.fig.generated.tex" % filename, "w") as f:
        f.write("\\begin{figure}[H]")
        f.write("\\label{f:%s}" % filename)
        f.write("\\begin{center}")
        f.write("\\includegraphics[scale=0.5]{%s}" % filename)
        if caption:
            f.write("\\caption{%s}" % caption)
        f.write("\\end{center}")
        f.write("\\end{figure}")

def gen_fixvars_tex(model, solution, fixvars):
    with open('fixvars.table.generated.tex', 'w') as f:
        f.write("\\begin{longtable}{llll}\n \\toprule\n")
        f.write("\\toprule\n")
        f.write("Variables & Value & Units & Description \\\\ \n")
        f.write("\\midrule\n")
        for varname in fixvars:
            val = "%0.3f" % solution(varname).magnitude
            unitstring = unitstr(model[varname].units)
            label = model[varname].descr["label"]
            f.write("$%s$ & %s & %s & %s \\\\\n" %
                    (varname, val, unitstring, label))
        f.write("\\bottomrule\n")
        f.write("\\end{longtable}\n")
