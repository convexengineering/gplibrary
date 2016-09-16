import numpy as np
import pandas as pd
from gasmale import GasMALE
from gpkit.small_scripts import unitstr
import xlsxwriter

def output_csv(path, M, sol, varnames, margins):
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
    for subm in M.submodels:
        if subm.__class__.__name__ == "Mission":
            for i, fs in enumerate(subm.submodels):
                fseg[fs.name]["index"].append(fs.num)
                fseg[fs.name]["shape"].append(fs.N)
                start.append(start[i] + fs.N)
                fseg[fs.name]["start"].append(start[i])

    colnames = ["Units"]
    for fs in fseg:
        for i in fseg[fs]["shape"]:
            colnames += [fs]*i
    colnames.append("Label")

    data = {}
    for vname in varnames:
        data[vname] = [0]*(start[-1] + 1)
        if vname in sens:
            data[vname + " sensitivity"] = [""] + [0]*(start[-1]-1) + [""]
        if vname in margins:
            data[vname + " margin"] = [""] + [0]*(start[-1]-1) + [""]
            data[vname+" margin sensitivity"] = [""] + [0]*(start[-1]-1) + [""]

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
                data[vname][st:st + sv.shape[0]] = sol(sv).magnitude[0:]
                if vname in sens:
                    data[vname+" sensitivity"][st:st+sv.shape[0]] = sens[sv][0:]
                if vname not in margins:
                    continue
                for mfac in sol("m_{fac}"):
                    # the varname must be in the m_fac label or it won't work
                    if not vname in mfac.label:
                        continue
                    data[vname+" margin"][st:st+sv.shape[0]] = \
                        [sol(mfac).magnitude]*sv.shape[0]
                    data[vname+" margin sensitivity"][st:st+sv.shape[0]] = \
                        [sens[mfac]]*sv.shape[0]


    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    df.to_csv("%soutput1.csv" % path)

def bd_csv_output(path, sol, varname):

    colnames = ["Value", "Units", "Margin", "Margin Sensitivity", "Label"]
    if varname in sol["sensitivities"]["constants"]:
        colnames.insert(2, "Sensitivitiy")

    data = {}
    for sv in sol(varname):
        name = max(list(sv.keys), key=len)
        data[name] = [sol(sv).magnitude]
        data[name].append(unitstr(sv.units))
        if varname in sol["sensitivities"]["constants"]:
            data[name].append(sol["sensitivities"]["constants"][sv])
        for mfac in sol("m_{fac}"):
            if not sv.models == mfac.models:
                continue
            data[name].append(sol(mfac).magnitude)
            data[name].append(sol["sensitivities"]["constants"][mfac])
        data[name].append(sv.label)

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    writer = pd.ExcelWriter("%s%s_breakdown.xlsx" % (path, varname.replace("{", "").replace("}", "")), engine="xlsxwriter")
    df.to_excel(writer, sheet_name="Sheet1")

    workbook = writer.book
    worksheet = writer.sheets["Sheet1"]
    # Add a format. Light red fill with dark red text.
    format1 = workbook.add_format({'bg_color': '#FFC7CE'})

    # Add a format. Green fill with dark green text.
    format2 = workbook.add_format({'bg_color': '#C6EFCE',
                                   'font_color': '#006100'})

    worksheet.conditional_format("E2:E%d" % (len(df["Margin Sensitivity"])+1),
                                 {"type": "cell",
                                  "criteria": ">=",
                                  "value": 0.8,
                                  "format": format1})
    # df.to_csv("%s%s_breakdown.csv" %
    #           (path, varname.replace("{", "").replace("}", "")))

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
    Margins = ["BSFC", "c_{dp}"]
    output_csv(PATH, M, Sol, Mission_vars, Margins)
    bd_csv_output(PATH, Sol, "W")
