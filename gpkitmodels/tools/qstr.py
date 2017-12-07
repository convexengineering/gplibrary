from gpkit.small_scripts import mag

def qstr(qty):
    'Converts a pint quantity to a nice string for display'
    return "%-.4g%s" % (qty.magnitude,  u" {:~}".format(qty.units)) if hasattr(qty, "magnitude") else (qty, "")
