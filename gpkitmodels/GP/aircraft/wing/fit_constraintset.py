from gpkit import ConstraintSet
from gpkit import Variable, NomialArray
import numpy as np

class FitCS(ConstraintSet):
    def __init__(self, df, ivar, dvars, nobounds=False, err_margin=False):

        K = int(df["K"].iloc[0])
        d = int(df["d"].iloc[0])
        ftype = df["ftype"].iloc[0]
        A = np.array(df[["e%d%d" % (k, i) for k in range(1, K+1) for i in
                         range(1, d+1)]])[0].astype(float)
        B = np.array(df[["c%d" % k for k in range(1, K+1)]])[0].astype(float)

        if np.array(dvars).ndim > 1:
            vvars = np.array(dvars).T
        else:
            vvars = np.array([dvars])
        monos = [B*NomialArray([(dv**A[k*d:(k+1)*d]).prod() for k in
                                range(K)]) for dv in vvars]

        if err_margin:
            maxerr = float(df["max_err"].iloc[0])
            mfac = Variable("m_{fac-fit}", maxerr, "-",
                            "max error of " + ivar.descr["label"] + " fit")
        else:
            mfac = 1

        if ftype == "ISMA":
            # constraint of the form 1 >= c1*u1^exp1*u2^exp2*w^(-alpha) + ....
            alpha = np.array(df[["a%d" % k for k in
                                 range(1, K+1)]])[0].astype(float)
            lhs = 1
            rhs = NomialArray([(mono/(ivar/mfac)**alpha).sum() for mono
                               in monos])
        elif ftype == "SMA":
            # constraint of the form w^alpha >= c1*u1^exp1 + c2*u2^exp2 +....
            alpha = float(df["a1"].iloc[0])
            lhs = (ivar/mfac)**alpha
            rhs = NomialArray([mono.sum() for mono in monos])
        elif ftype == "MA":
            # constraint of the form w >= c1*u1^exp1, w >= c2*u2^exp2, ....
            lhs, rhs = (ivar/mfac), NomialArray(monos)

        if K == 1:
            # when possible, return an equality constraint
            cstrt = (lhs == rhs)
        else:
            cstrt = (lhs >= rhs)

        constraints = [cstrt]

        if not nobounds:
            self.boundingvars = []
            for i, v in enumerate(dvars):
                if hasattr(v, "__len__"):
                    desv = v[0]
                else:
                    desv = v
                if "units" in desv.descr:
                    unt = desv.descr["units"]
                else:
                    unt = "-"
                low = Variable(desv.descr["name"] + "_{low-bound}",
                               float(df["lb%d" % (i+1)].iloc[0]), unt,
                               desv.descr["label"] + " lower bound")
                up = Variable(desv.descr["name"] + "_{up-bound}",
                              float(df["ub%d" % (i+1)].iloc[0]), unt,
                              desv.descr["label"] + " upper bound")

                self.boundingvars.extend(np.hstack([low, up]))

                constraints.extend([v >= low,
                                    v <= up])
            self.boundingvars = np.hstack(self.boundingvars)

        else:
            self.boundingvars = None

        ConstraintSet.__init__(self, constraints)

    def process_result(self, result, TOL=1e-5):
        super(FitCS, self).process_result(result)

        if not self.boundingvars is None:
            for var in self.boundingvars:
                sen = result["sensitivities"]["constants"][var]

                if abs(sen) > TOL:
                    msg = ("Variable %.100s could cause inaccurate result"
                           " because sensitivity to " % var
                           + var.descr["label"] + " of value %.3f is greater"
                           " than 0, with a sensitivity of %.5f" %
                           (var.descr["value"], sen))
                    print "Warning: " + msg

