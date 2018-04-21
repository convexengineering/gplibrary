" tail boom flexibility "
from numpy import pi
from gpkit import Model, parse_variables, SignomialsEnabled

class TailBoomFlexibility(Model):
    """ Tail Boom Flexibility Model

    Variables
    ---------
    Fne                 [-]     tail boom flexibility factor
    deda                [-]     wing downwash derivative
    SMcorr      0.55    [-]     corrected static margin
    sph1                [-]     flexibility helper variable 1
    sph2                [-]     flexibility helper variable 2

    LaTex Strings
    -------------
    Fne         F_{\mathrm{NE}}
    deda        d\\epsilon/d\\alpha
    SMcorr      SM_{\\mathrm{corr}}

    """
    def setup(self, htail, hbending, wing):
        exec parse_variables(TailBoomFlexibility.__doc__)

        mh = htail.mh
        mw = wing.mw
        Vh = htail.Vh
        th = hbending.th
        CLhmin = htail.CLhmin
        CLwmax = wing.planform.CLmax
        Sw = wing.planform.S
        bw = wing.planform.b
        lh = htail.lh
        CM = wing.planform.CM

        constraints = [
            Fne >= 1 + mh*th,
            sph1*(mw*Fne/mh/Vh) + deda <= 1,
            sph2 <= Vh*CLhmin/CLwmax,
            # (sph1 + sph2).mono_lower_bound({"sph1": .48, "sph2": .52}) >= (
            #     SMcorr + wing["C_M"]/wing["C_{L_{max}}"]),
            deda >= mw*Sw/bw/4/pi/lh]

        with SignomialsEnabled():
            constraints.extend([sph1 + sph2 >= SMcorr + CM/CLwmax])

        return constraints
