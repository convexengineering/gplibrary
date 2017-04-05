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

        monos = B*NomialArray([(dvars**A[k*d:(k+1)*d]).prod() for k in
                               range(K)])

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
            lhs, rhs = 1, (monos/(ivar/mfac)**alpha).sum()
        elif ftype == "SMA":
            # constraint of the form w^alpha >= c1*u1^exp1 + c2*u2^exp2 +....
            alpha = float(df["a1"].iloc[0])
            lhs, rhs = (ivar/mfac)**alpha, monos.sum()
        elif ftype == "MA":
            # constraint of the form w >= c1*u1^exp1, w >= c2*u2^exp2, ....
            lhs, rhs = (ivar/mfac), monos

        if K == 1:
            # when possible, return an equality constraint
            cstrt = (lhs == rhs)
        else:
            cstrt = (lhs >= rhs)

        constraints = [cstrt]

        if not nobounds:
            self.boundingvars = []
            for i, v in enumerate(dvars):
                if "units" in v.descr:
                    unt = v.descr["units"]
                else:
                    unt = "-"
                low = Variable(v.descr["name"] + "_{low-bound}",
                               float(df["lb%d" % (i+1)].iloc[0]), unt,
                               v.descr["label"] + " lower bound")
                up = Variable(v.descr["name"] + "_{up-bound}",
                              float(df["ub%d" % (i+1)].iloc[0]), unt,
                              v.descr["label"] + " upper bound")

                self.boundingvars.extend([low, up])

                constraints.extend([v >= low,
                                    v <= up])

        else:
            self.boundingvars = None

        ConstraintSet.__init__(self, constraints)

    def process_result(self, result):
        super(FitCS, self).process_result(result)

        if self.boundingvars:
            for var in self.boundingvars:
                sen = result["sensitivities"]["constants"][var]
                if hasattr(sen, "__len__"):
                    sen = max(np.abs(sen.values()))

                if sen > 0.01:
                    msg = ("Sensitive to " + var.descr["label"] +
                           " of value %.3f with a sensitivity of %.3f\n" %
                           (var.descr["value"], sen))
                    print "Warning: " + msg

