import numpy as np
import pandas as pd

def output_csv(M, sol):

    data = {}
    values = []
    units = []
    variables = []
    for var in sol["variables"]:
        values.append(sol(var).magnitude)
        units.append(sol(var).units)
        variables.append(var)

    data["Variables"] = variables
    data["Value"] = values
    data["units"] = units

    columns = ["Variables", "Value", "units"]

    df = pd.DataFrame(data, columns=columns)
    df.to_csv("csv/output1.csv")
