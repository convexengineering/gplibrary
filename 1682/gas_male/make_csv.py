import numpy as np
import pandas as pd
from gasmale import GasMALE
from gpkit.small_scripts import unitstr
import xlsxwriter

def mission_vars(M, sol, varnames, margins):
    """
    This ouputs variables relevant accross a mission
    """
    sens = sol["sensitivities"]["constants"]

    fseg = {}
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            for fs in subm.submodels:
                fseg[fs.name] = {"index": [], "shape": [], "start": []}

    start = [1]
    colnames = ["Units"]
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            for i, fs in enumerate(subm.submodels):
                fseg[fs.name]["index"].append(fs.num)
                fseg[fs.name]["shape"].append(fs.N)
                start.append(start[i] + fs.N)
                fseg[fs.name]["start"].append(start[i])
                colnames += [fs.name + "%d.%d" %
                             (fs.num, n) for n in range(fs.N)]
    colnames.append("Label")

    data = {}
    for vname in varnames:
        data[vname] = [0]*(start[-1] + 1)
        if vname in sens:
            data[vname + " Sens"] = [""] + [0]*(start[-1]-1) + [""]
        if vname in margins:
            data[vname + " Margin"] = [""] + [0]*(start[-1]-1) + [""]
            data[vname+" Margin Sens"] = [""] + [0]*(start[-1]-1) + [""]

    i = 0
    for vname in varnames:
        for sv in sol(vname):
            for fs in fseg:
                if fs not in sv.models:
                    continue
                ind = sv.models.index(fs)
                ifs = fseg[fs]["index"].index(sv.modelnums[ind])
                st = fseg[fs]["start"][ifs]
                data[vname][0] = unitstr(sv.units)
                data[vname][-1] = sv.label
                if "shape" in sv.descr:
                    data[vname][st:st + sv.shape[0]] = sol(sv).magnitude[0:]
                    if vname in sens:
                        data[vname+" Sens"][st:st+sv.shape[0]] = sens[sv][0:]
                else:
                    data[vname][st:st + fseg[fs]["shape"][ifs]] = \
                            [sol(sv).magnitude]*fseg[fs]["shape"][ifs]
                    if vname in sens:
                        data[vname+" Sens"][st:st+fseg[fs]["shape"][ifs]] = \
                                [sens[sv]]*fseg[fs]["shape"][ifs]
                if vname not in margins:
                    continue
                for mfac in sol("m_{fac}"):
                    # the varname must be in the m_fac label or it won't work
                    if not vname in mfac.label:
                        continue
                    data[vname+" Margin"][st:st+sv.shape[0]] = \
                        [sol(mfac).magnitude]*sv.shape[0]
                    data[vname+" Margin Sens"][st:st+sv.shape[0]] = \
                        [sens[mfac]]*sv.shape[0]


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

if __name__ == "__main__":
    M = GasMALE(DF70=True)
    M.substitutions.update({"t_{loiter}": 6})
    M.cost = M["MTOW"]
    Sol = M.solve("mosek")
    PATH = "/Users/mjburton11/Dropbox (MIT)/16.82GasMALE/GpkitReports/csvs/"

    Mission_vars = ["RPM", "RPM_{max}", "BSFC", "V", "P_{shaft}",
                    "P_{shaft-tot}", "h_{dot}", "h", "T_{atm}", "\\mu",
                    "\\rho", "W_{fuel}", "W_{N}", "W_{N+1}", "C_D", "C_L",
                    "\\eta_{prop}", "T", "h_{loss}", "P_{shaft-max}", "t",
                    "Re", "C_{f-fuse}", "C_{D-fuse}", "c_{dp}", "V_{wind}"]
    Margins = ["BSFC", "c_{dp}"]
    Sens_boundaries = {"bad": 0.8, "good": 0.2}
    DF = mission_vars(M, Sol, Mission_vars, Margins)
    DF.to_csv("test.csv")
    write_to_excel(PATH, "Mission_params.xlsx", DF, Sens_boundaries)
    DF = bd_vars(M, Sol, "W", ["MTOW", "W_{fuel-tot}", "W_{zfw}"])
    write_to_excel(PATH, "W_breakdown.xlsx", DF, Sens_boundaries)
