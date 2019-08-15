from __future__ import print_function
from builtins import range
import subprocess
import numpy as np
import scipy.optimize as spo
import math
import sys
from gpkit.tests.helpers import NullFile

def blind_call(topline, cl, Re, M, max_iter = 100,
               pathname = "/usr/local/bin/xfoil"):

    if '.dat' in topline:
        tl=topline.split()
        afile=tl[1]
    if '.txt' in topline:
        tl=topline.split()
        afile=tl[1]

    proc = subprocess.Popen([pathname], stdout=subprocess.PIPE, stdin=subprocess.PIPE)
    proc.stdin.write(topline +
                     'oper \n' +
                     "iter %d\n" %(max_iter)+
                     'visc \n' +
                     "%.2e \n" %(Re) +
                     "M \n" +
                     "%.2f \n" %(M) +
                     "a 2.0 \n" +
                     "cl %.4f \n" %(cl) +
                     '\n' +
                     'quit \n')
    stdout_val = proc.communicate()[0]
    proc.stdin.close()

    if ("VISCAL:  Convergence failed\n" in stdout_val):
        return stdout_val

    res = {}
    if ("VISCAL:  Convergence failed\n" not in stdout_val):
        ostr = stdout_val.split()
        ctr = 0
        for i in range(0,len(ostr)):
            ix = len(ostr)-(i+1)
            vl = ostr[ix]
            if vl in ['a','CL','CD','Cm']:
                res[vl] = ostr[ix + 2]
                ctr += 1
            if ctr >= 4:
                break
        cd = res['CD']
        cl = res['CL']
        # alpha_ret = res['a']
        cm = res['Cm']
        return float(cd), float(cl), float(cm), stdout_val

def single_cl(CL, Re = 1e7, M = 0.0, airfoil=[], pathname = "/home/ckarcher/Xfoil/bin/./xfoil",
                    number_of_samples = 51, sampling_min=-10, sampling_max=20, fitting_fraction = 1.4):

    num_samples = number_of_samples
    sample_min=sampling_min
    sample_max=sampling_max

    ls_res = subprocess.check_output(["ls -a"], shell = True)
    ls = ls_res.split()

    remove_kulfan = False
    if list(airfoil):
        if ('.dat' in airfoil) or ('.txt' in airfoil):
            topline = 'load ' + airfoil + ' \n afl \n'
        elif ('naca' == airfoil.lower()[0:4]) and (len(airfoil)==8):
            topline = airfoil + ' \n'
    else:
        print("Error: Invalid airfoil passed into XFOIL.  Defaulting to a NACA0012.")
        topline = 'naca0012 \n'

    initial_list = np.linspace(sample_min,sample_max,num_samples).tolist()
    cd_calcl = []
    cl_calcl = []
    alpha_calcl = []
    cm_calcl = []
    for alpha in initial_list:
        # x = blind_call(topline, alpha, Re, M)
        try:
            x = blind_call(topline, alpha, Re, M)
        except:
            x=[1,1]
        if len(x)==5:
            cd_calcl.append(x[0])
            cl_calcl.append(x[1])
            alpha_calcl.append(x[2])
            cm_calcl.append(x[3])
        elif len(x)>10:
            pass

    cd_calc = np.asarray(cd_calcl)
    cl_calc = np.asarray(cl_calcl)
    alpha_calc = np.asarray(alpha_calcl)
    cm_calc = np.asarray(cm_calcl)

    vldidx = np.where(cl_calc<=max(cl_calc))
    cd_calc=cd_calc[vldidx[0].tolist()]
    cl_calc=cl_calc[vldidx[0].tolist()]
    alpha_calc=alpha_calc[vldidx[0].tolist()]
    cm_calc=cm_calc[vldidx[0].tolist()]

    vldidx2 = np.where(cl_calc>=0.0)
    cd_calc=cd_calc[vldidx2[0].tolist()]
    cl_calc=cl_calc[vldidx2[0].tolist()]
    alpha_calc=alpha_calc[vldidx2[0].tolist()]
    cm_calc=cm_calc[vldidx2[0].tolist()]

    p_cd  = np.polyfit(np.append(-cl_calc,cl_calc), np.append(cd_calc,cd_calc),int(len(cl_calc)/fitting_fraction))
    p_alpha = np.polyfit(np.append(-cl_calc,cl_calc), np.append(alpha_calc,alpha_calc),int(len(cl_calc)/fitting_fraction))
    p_cm = np.polyfit(np.append(-cl_calc,cl_calc), np.append(cm_calc,cm_calc),int(len(cl_calc)/fitting_fraction))

    cd_guess = np.polyval(p_cd,CL)
    alpha_guess = np.polyval(p_alpha,CL)
    cm_guess = np.polyval(p_cm,CL)

    return cd_guess, CL, alpha_guess, cm_guess
