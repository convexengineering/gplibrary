import numpy as np
import pandas as pd
from gasmale import GasMALE
from gpkit.small_scripts import unitstr

def output_csv(PATH, M, sol, varnames):
    """
    This ouputs variables relevant accross a mission
    """

    fseg = {"name": [], "shape": [], "index": []}
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            for fs in subm.submodels:
                fseg["name"].append(fs.__class__.__name__)
                fseg["shape"].append(fs.N)
                fseg["index"].append(fs.num)

    fseg["start"] = [0]
    for i in range(0, len(fseg["shape"])-1):
        fseg["start"].append(fseg["start"][i] + fseg["shape"][i])

    colnames = ["Units"]
    for n, s in zip(fseg["name"], fseg["shape"]):
        for i in range(s):
            colnames.append(n)
    colnames.append("Label")

    data = {}
    for vname in varnames:
        data[vname] = [0]*(sum(fseg["shape"]) + 2)

    for var in M.varkeys:
        if var.descr["name"] in data:
            mname = list(set(fseg["name"]) & set(var.descr["models"]))[0]
            ind = var.descr["modelnums"][var.descr["models"].index(mname)]
            if "shape" in var.descr:
                shape = var.descr["shape"][0]
            for n, i, s, l in zip(fseg["name"], fseg["index"],
                                  fseg["shape"], fseg["start"]):
                if i == ind and n == mname and s == shape:
                    data[var.descr["name"]][0] = unitstr(var.units)
                    data[var.descr["name"]][-1] = var.label
                    if "idx" in var.descr:
                        data[var.descr["name"]][var.descr["idx"][0] + \
                                l + 1] = sol(var).magnitude
                    else:
                        data[var.descr["name"]][l + 1] = [sol(var).magnitude]

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    df.to_csv("%soutput1.csv" % PATH)

def bd_csv_output(PATH, sol, varname):

    if varname in sol["sensitivities"]["constants"]:
        colnames = ["Value", "Units", "Sensitivitiy", "Label"]
    else:
        colnames = ["Value", "Units", "Label"]

    data = {}
    for sv in sol(varname):
        name = sv
        data[name] = [sol(sv).magnitude]
        data[name].append(unitstr(sv.units))
        if varname in sol["sensitivities"]["constants"]:
            data[name].append(sol["sensitivities"]["constants"][sv])
        data[name].append(sv.label)

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    df.to_csv("%s%s_breakdown.csv" %
              (PATH, varname.replace("{", "").replace("}", "")))

if __name__ == "__main__":
    M = GasMALE()
    M.substitutions.update({"t_{loiter}": 6})
    M.cost = M["MTOW"]
    Sol = M.solve("mosek")
    PATH = "/Users/mjburton11/Dropbox (MIT)/16.82GasMALE/GpkitReports/csvs/"

    Mission_vars = ["RPM", "BSFC", "V", "P_{shaft}", "P_{shaft-tot}",
                    "h_{dot}", "h", "T_{atm}", "\\mu", "\\rho", "W_{fuel}",
                    "W_{N}", "W_{N+1}", "C_D", "C_L", "\\eta_{prop}", "T",
                    "h_{loss}", "P_{shaft-max}", "t", "Re", "C_{f-fuse}",
                    "C_{D-fuse}", "c_{dp}", "V_{wind}"]
    output_csv(PATH, M, Sol, Mission_vars)
    bd_csv_output(PATH, Sol, "W")
    bd_csv_output(PATH, Sol, "m_{fac}")
