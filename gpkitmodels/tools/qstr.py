from gpkit.small_scripts import mag

def qstr(qty, fmtstr="%-.4g"):
    'Converts a number/pint quantity to a nice string for display'
    if hasattr(qty, "magnitude"):
        magstr = fmt % qty.magnitude
        unitstr = u" {:~}".format(qty.units))
        if "dimensionless" in unitstr:  # boring and long way to put it
            unitstr = ""
        return magstr + unitstr
    else:
        return fmt % qty  # it's just a number
