import numpy as np
import pandas as pd
from gasmale import GasMALE
from gpkit.small_scripts import unitstr
from gpkit import Variable
from gen_tex import find_submodels
import xlsxwriter

def mission_vars(M, sol, varnames, margins):
    """
    This ouputs variables relevant accross a mission
    """
    sens = sol["sensitivities"]["constants"]

    data = {}
    n = []
    colnames = ["Units"]
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            mission = subm
            for fs in mission.submodels:
                n.append(fs.N)
                for i in range(fs.N):
                    colnames.append(fs.__class__.__name__ +
                                    "%s.%s" % (fs.num, i))
    colnames.append("Label")

    for varname in varnames:
        data[varname] = [""]
        for flightseg in mission.submodels:
            if varname in flightseg.varkeys:
                if flightseg[varname] in sens:
                    data[varname + " sens"] = [""]
        for i, fs in enumerate(mission.submodels):
            if varname not in fs.varkeys:
                data[varname].append([""]*n[i])
                if varname+" sens" in data:
                    data[varname+" sens"].append([""]*n[i])
                continue
            nonvector = isinstance(fs[varname], Variable)
            units = (unitstr(fs[varname][0].descr["units"]) if not nonvector
                     else unitstr(fs[varname].descr["units"]))
            data[varname].append(sol(fs[varname]).magnitude if not nonvector
                                 else [sol(fs[varname]).magnitude]*n[i])
            if fs[varname] in sens:
                if not nonvector:
                    data[varname+" sens"].append(sens[fs[varname]])
                else:
                    data[varname+" sens"].append(list(sens[fs[varname]])*n[i])
            if i == len(mission.submodels)-1:
                data[varname].append(fs[varname][0].descr["label"] if not
                                     nonvector else fs[varname].descr["label"])
        data[varname][0] = units
        data[varname] = np.hstack(data[varname])

    for d in data:
        if "sens" in d or len(data[d]) < sum(n)+2:
            data[d] = np.hstack(data[d])
            data[d] = np.append(data[d], [""])

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    return df

def bd_vars(M, sol, varname, morevars):

    colnames = ["Value", "Units", "Margin", "Margin Sens", "Label"]

    data = {}
    for sv in sol(varname):
        name = max(list(sv.keys), key=len)
        data[name] = [sol(sv).magnitude]
        data[name].append(unitstr(sv.units))
        for mfac in sol("m_{fac}"):
            if not sv.models == mfac.models:
                continue
            data[name].append(sol(mfac).magnitude)
            data[name].append(sol["sensitivities"]["constants"][mfac])
        if len(data[name]) != 4:
            data[name] += [""]*2
        data[name].append(sv.label)

    for name in morevars:
        sv = sol(name)
        data[name] = [sv.magnitude]
        data[name].append(unitstr(M[name].descr["units"]))
        for mfac in sol("m_{fac}"):
            if not M[name].descr["models"] == mfac.models:
                continue
            data[name].append(sol(mfac).magnitude)
            data[name].append(sol["sensitivities"]["constants"][mfac])
        if len(data[name]) != 4:
            data[name] += [""]*2
        data[name].append(M[name].descr["label"])

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    return df

def sketch_params(M, sol, varnames, othervars=None, pointmasses=None):

    data = {}
    for vname in varnames:
        data[vname] = [sol(vname).magnitude, unitstr(M[vname].descr["units"]),
                       M[vname].descr["label"]]

    if othervars:
        data.update(othervars)

    if hasattr(M, "get_cgs"):
        xnp, xcg, SM = M.get_cgs()
        data["x_{np}"] = [xnp.magnitude, xnp.units, "neutral point"]
        data["x_{cg}"] = [xcg.magnitude, xcg.units, "center of gravity"]
        data["SM"] = [SM.magnitude, "-", "static margin"]

    if pointmasses:
        for pm in pointmasses:
            data[pm] = []

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = ["Value", "Units", "Label"]
    return df

def write_to_excel(path, filename, df, sens_formatting):

    coldepth = df.count()[0]
    rowdepth = len(df.columns)

    colind = []
    for colname in df.columns:
        if "Sens" in colname:
            colind.append(list(df.columns).index(colname))

    rowind = []
    for rowname in df.index:
        if "Sens" in rowname:
            rowind.append(list(df.index).index(rowname))

    alp = ["A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K", "L", "M"]

    writer = pd.ExcelWriter("%s%s" % (path, filename), engine="xlsxwriter")
    df.to_excel(writer, sheet_name="Sheet1")

    workbook = writer.book
    ws = writer.sheets["Sheet1"]
    # Light red fill
    format1 = workbook.add_format({'bg_color': '#FFC7CE'})

    # Green fill
    format2 = workbook.add_format({'bg_color': '#FFCC99'})

    # Green fill
    format3 = workbook.add_format({'bg_color': '#C6EFCE'})

    for i in colind:
        ws.conditional_format("%s2:%s%d" % (alp[i+1], alp[i+1], coldepth+1),
                              {"type": "cell",
                               "criteria": ">=",
                               "value": sens_formatting["bad"],
                               "format": format1})

        ws.conditional_format("%s2:%s%d" % (alp[i+1], alp[i+1], coldepth+1),
                              {"type": "cell",
                               "criteria": "between",
                               "minimum": sens_formatting["bad"],
                               "maximum": sens_formatting["good"],
                               "format": format2})

        ws.conditional_format("%s2:%s%d" % (alp[i+1], alp[i+1], coldepth+1),
                              {"type": "cell",
                               "criteria": "<",
                               "value": sens_formatting["good"],
                               "format": format3})

    for i in rowind:
        ws.conditional_format("%s%d:%s%d" % (alp[2], i+2, alp[rowdepth-1], i+2),
                              {"type": "cell",
                               "criteria": ">=",
                               "value": sens_formatting["bad"],
                               "format": format1})

        ws.conditional_format("%s%d:%s%d" % (alp[2], i+2, alp[rowdepth-1], i+2),
                              {"type": "cell",
                               "criteria": "between",
                               "minimum": sens_formatting["bad"],
                               "maximum": sens_formatting["good"],
                               "format": format2})

        ws.conditional_format("%s%d:%s%d" % (alp[2], i+2, alp[rowdepth-1], i+2),
                              {"type": "cell",
                               "criteria": "<",
                               "value": sens_formatting["good"],
                               "format": format3})

    writer.save()

def model_params(subM, sol):

    data = {}
    for v in subM.varkeys:
        if "Mission" not in v.descr["models"]:
            data[v] = [sol(v).magnitude]
            data[v].append(unitstr(M[v].units))
            data[v].append(v.descr["label"])

    if data:
        df = pd.DataFrame(data)
        df = df.transpose()
        df.columns = ["Value", "Units", "Label"]
    else:
        df = None
    return df

if __name__ == "__main__":
    M = GasMALE(DF70=True)
    M.substitutions.update({"t_{loiter}": 6})
    M.cost = M["MTOW"]
    Sol = M.solve("mosek")
    PATH = "/Users/mjburton11/Dropbox (MIT)/16.82GasMALE/Management/GpkitReports/"

    Mission_vars = ["RPM", "BSFC", "V", "P_{shaft}",
                    "P_{shaft-tot}", "h_{dot}", "h", "T_{atm}", "\\mu",
                    "\\rho", "W_{fuel}", "W_{N}", "W_{N+1}", "C_D", "C_L",
                    "\\eta_{prop}", "T", "h_{loss}", "P_{shaft-max}", "t",
                    "C_{f-fuse}", "C_{D-fuse}", "c_{dp}", "V_{wind}"]
    Margins = ["BSFC", "c_{dp}"]
    Sens_boundaries = {"bad": 0.8, "good": 0.2}
    DF = mission_vars(M, Sol, Mission_vars, Margins)
    DF.to_csv("test.csv")
    write_to_excel(PATH, "Mission_params.xlsx", DF, Sens_boundaries)
    DF = bd_vars(M, Sol, "W", ["MTOW", "W_{fuel-tot}", "W_{zfw}"])
    write_to_excel(PATH, "W_breakdown.xlsx", DF, Sens_boundaries)

    m, mn = find_submodels([M], [])
    mn = [""] + mn
    for subm, name in zip(m, mn):
        df = model_params(subm, Sol)
        if df is not None:
            df.to_csv(PATH + "%s.csv" % name)

    # df = sketch_params(
    #     M, Sol, ["S", "b", "l_{fuel}", "d", "L", "S_h", "S_v", "b_h", "b_v", "d_0"],
    #     othervars={"lambda":[0.5, "-", "taper ratio"],
    #                "eta":[3, "-", "tail boom separation/fuselage diameter"]}
    #     )
    # df.to_csv("sketch_params.csv")
