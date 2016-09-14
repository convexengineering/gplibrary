import numpy as np
import pandas as pd
from gasmale import GasMALE

def output_csv(M, sol):

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

    colnames = []
    for n, s in zip(fseg["name"], fseg["shape"]):
        for i in range(s):
            colnames.append(n)

    data = {}
    for var in M.varkeys:
        for mname in var.descr["models"]:
            if mname in fseg["name"] and var.descr["name"] not in data:
                data[var.descr["name"]] = np.zeros(sum(fseg["shape"]))

    for var in M.varkeys:
        if var.descr["name"] in data:
            mname = list(set(fseg["name"]) & set(var.descr["models"]))[0]
            ind = var.descr["modelnums"][var.descr["models"].index(mname)]
            if "shape" in var.descr:
                shape = var.descr["shape"][0]
            for n, i, s, l in zip(fseg["name"], fseg["index"],
                                  fseg["shape"], fseg["start"]):
                if i == ind and n == mname and s == shape:
                    if "idx" in var.descr:
                        data[var.descr["name"]][var.descr["idx"][0] + \
                                l] = sol(var).magnitude
                    else:
                        data[var.descr["name"]] = sol(var).magnitude

    df = pd.DataFrame(data)
    df = df.transpose()
    df.columns = colnames
    df.to_csv("csv/output1.csv")

if __name__ == "__main__":
    M = GasMALE()
    M.substitutions.update({"t_{loiter}": 6})
    M.cost = M["MTOW"]
    sol = M.solve("mosek")

    output_csv(M, sol)
