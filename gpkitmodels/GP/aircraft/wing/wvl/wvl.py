" Weissenger vortex lattice method "
import numpy as np
from numpy import arange, array, interp, pi, sqrt, zeros, flip, cos, sin, dot
from numpy import diff
from numpy.linalg import lstsq
from vorvel import vorvel

def wvl(geom, N, ispace, Sref, bref, cref, itmax, toler, alspec, bespec, pbspec, rbspec, CLspec, Crspec, Cnspec, ialspec, ibespec, ipbspec, irbspec, iCLspec, iCrspec, iCnspec):

    ns = geom.shape[0]

    # calculate running length along span
    ssec = zeros([ns, 1])
    for i in range(ns-1):
        dy = geom[i+1, 1] - geom[i, 1]
        dz = geom[i+1, 2] - geom[i, 2]
        ssec[i+1] = ssec[i] + sqrt(dy**2 + dz**2)

    # number of flow parameters
    nflow = 1
    npar = 4

    if ispace == 1:
        dt = ssec[ns-1]/N
        t = arange(0, ssec[ns-1] + dt, dt)
        s = arange(0, ssec[ns-1] + dt, dt)
    else:
        dt = 0.5*pi/N
        t = arange(0, 0.5*pi + dt, dt)
        s = ssec[ns-1]*cos(0.5*pi - t)

    c = interp(s, ssec[:, 0], geom[:, 3])
    x = interp(s, ssec[:, 0], geom[:, 0]) + 0.25*c
    y = interp(s, ssec[:, 0], geom[:, 1])
    z = interp(s, ssec[:, 0], geom[:, 2])

    if ispace == 1:
        tv = t[1:] - 0.5*dt
        sv = s[1:] - 0.5*dt
    else:
        tv = t[1:] - 0.5*dt
        sv = ssec[ns-1]*cos(0.5*pi - tv)

    cv = interp(sv, ssec[:, 0], geom[:, 3])
    xv = interp(sv, ssec[:, 0], geom[:, 0]) + 0.25*cv
    yv = interp(sv, ssec[:, 0], geom[:, 1])
    zv = interp(sv, ssec[:, 0], geom[:, 2])

    xc = xv + 0.5*cv
    yc = yv
    zc = zv

    sind = zeros(N)
    cosd = zeros(N)
    for i in range(N):
        sind[i] = (z[i+1] - z[i])/(s[i+1]-s[i])
        cosd[i] = (y[i+1] - y[i])/(s[i+1]-s[i])

    twv = interp(sv, ssec[:, 0], geom[:, 4])
    a0v = interp(sv, ssec[:, 0], geom[:, 5])
    twa = twv - a0v

    c = array(list(flip(c, 0)) + list(c[1:]))
    x = array(list(flip(x, 0)) + list(x[1:]))
    y = array(list(flip(-y, 0)) + list(y[1:]))
    z = array(list(flip(z, 0)) + list(z[1:]))
    s = array(list(flip(-s, 0)) + list(s[1:]))

    cv = array(list(flip(cv, 0)) + list(cv))
    xv = array(list(flip(xv, 0)) + list(xv))
    yv = array(list(flip(-yv, 0)) + list(yv))
    zv = array(list(flip(zv, 0)) + list(zv))

    xc = array(list(flip(xc, 0)) + list(xc))
    yc = array(list(flip(-yc, 0)) + list(yc))
    zc = array(list(flip(zc, 0)) + list(zc))

    twa = array(list(flip(twa, 0)) + list(twa))
    cost = cos(twa)
    sint = sin(twa)

    cosd = array(list(flip(cosd, 0)) + list(cosd))
    sind = array(list(flip(-sind, 0)) + list(sind))

    wyG = zeros([2*N, 2*N])
    wzG = zeros([2*N, 2*N])

    for i in range(2*N):
        rsqi = (yv-y[i])**2 + (zv-z[i])**2
        rsqp = (yv-y[i+1])**2 + (zv-z[i+1])**2
        wyG[i, :] = ((zv-z[i])/rsqi - (zv-z[i+1])/rsqp)/(4*pi)
        wzG[i, :] = -((yv-y[i])/rsqi - (yv-y[i+1])/rsqp)/(4*pi)

    dx = diff(x)
    dy = diff(y)
    dz = diff(z)

    A = zeros([2*N, 2*N])

    for i in range(2*N):
        vvor = vorvel(array([x[i], y[i], z[i]]),
                      array([x[i+1], y[i+1], z[i+1]]),
                      array([xc, yc, zc]), 1)
        A[i, :] = vvor[0, :]*sint + vvor[2, :]*cost

    alpha = 0.0
    beta = 0.0
    pbar = 0.0
    rbar = 0.0

    dalpha = 1.0
    dbeta = 1.0
    dpbar = 1.0
    drbar = 1.0

    itr = 0
    istop = 0

    while max([abs(dalpha), abs(dbeta), abs(dpbar), abs(drbar)]) > toler:
        print "\n Iteration %5i ...\n" % itr
        print " alpha=%7.3f" % (alpha*180/pi)
        print " beta=%7.3f" % (beta*180/pi)

        if (itr >= itmax):
            print "\n Iteration limit exceeded\n"
            istop = 1
            break

        if (abs(alpha)*180/pi > 45.0):
            print 'Alpha > 45 deg reached. CLspec is probably excessive.\n'
            istop = 2
            break

        if (abs(beta)*180/pi > 45.0):
            print 'Beta > 45 deg, likely due to insufficient dihedral.\n'
            istop = 2
            break

        sina = sin(alpha)
        cosa = cos(alpha)
        sinb = sin(beta)
        cosb = cos(beta)

        # set up rhs vectors for each flow
        R = zeros([nflow, 2*N])
        R_a = zeros([nflow, 2*N])
        R_b = zeros([nflow, 2*N])
        R_p = zeros([nflow, 2*N])
        R_r = zeros([nflow, 2*N])

        Vx = cosa*cosb - 2*yc*rbar
        Vy = -sinb - 2*zc*pbar
        Vz = sina*cosb + 2*yc*pbar

        Vx_a = -sina*cosb
        Vy_a = 0.0
        Vz_a = cosa*cosb

        Vx_b = -cosa*sinb
        Vy_b = -cosb
        Vz_b = -sina*sinb

        Vx_p = 0.0
        Vy_p = -zc
        Vz_p = yc

        Vx_r = -yc
        Vy_r = 0.0
        Vz_r = 0.0

        R = -(Vx*sint + Vz*cost*cosd - Vy*sind)
        R_a = -(Vx_a*sint + Vz_a*cost*cosd - Vy_a*sind)
        R_b = -(Vx_b*sint + Vz_b*cost*cosd - Vy_b*sind)
        R_p = -(Vx_p*sint + Vz_p*cost*cosd - Vy_p*sind)
        R_r = -(Vx_r*sint + Vz_r*cost*cosd - Vy_r*sind)

        G = lstsq(A.T, R.T)[0]
        G_a = lstsq(A.T, R_a.T)[0]
        G_b = lstsq(A.T, R_b.T)[0]
        G_p = lstsq(A.T, R_p.T)[0]
        G_r = lstsq(A.T, R_r.T)[0]

        vi = dot(G, wyG)
        wi = dot(G, wzG)

        vi_a = dot(G_a, wyG)
        vi_b = dot(G_b, wyG)
        vi_p = dot(G_p, wyG)
        vi_r = dot(G_r, wyG)

        wi_a = dot(G_a, wzG)
        wi_b = dot(G_b, wzG)
        wi_p = dot(G_p, wzG)
        wi_r = dot(G_r, wzG)

        cl = zeros(2*N)  # local L'/(0.5 rho V^2 c)        =  cl
        ccl = zeros(2*N)  # local L'/(0.5 rho Vinf^2 cref)
                                # =  cl*(c/cref)*(V/Vinf)^2   (local loading)
        CL = 0.0  # overall CL
        CDi = 0.0  # overall CDi
        Cr = 0.0  # overall Cr  (roll moment coefficient)
        Cn = 0.0  # overall Cn  (yaw  moment coefficient)
        Cb = 0.0  # root bending moment coefficient

        CL_a = 0.0
        CL_b = 0.0
        CL_p = 0.0
        CL_r = 0.0

        Cr_a = 0.0
        Cr_b = 0.0
        Cr_p = 0.0
        Cr_r = 0.0

        Cn_a = 0.0
        Cn_b = 0.0
        Cn_p = 0.0
        Cn_r = 0.0

        Vx = zeros(2*N)

        for i in range(2*N):
            # % normalized local velocity relative to wing station,
            # [Vx,Vy,Vz] = [u,v,w]/Vinf
            Vx[i] = cosa*cosb - yc[i]*rbar
            Vy = -sinb - zc[i]*pbar
            Vz = sina*cosb + yc[i]*pbar

            Vx_a = -sina*cosb
            Vy_a = 0.0
            Vz_a = cosa*cosb

            Vx_b = -cosa*sinb
            Vy_b = -cosb
            Vz_b = -sina*sinb

            Vx_p = 0.0
            Vy_p = -zc[i]
            Vz_p = yc[i]

            Vx_r = -yc[i]
            Vy_r = 0.0
            Vz_r = 0.0

            cl[i] = 2*G[i]/(Vx[i]*cv[i])
            ccl[i] = 2*G[i]*Vx[i]/cref

            CL = CL + 2*G[i]* Vx[i]*dy[i]/Sref
            Cr = Cr - 2*G[i]* Vx[i]*(yv[i]*dy[i]+zv[i]*dz[i])/(bref*Sref)
            Cn = Cn - 2*G[i]*(Vz+wi[i])* yv[i]*dy[i]/(bref*Sref)

            Cb = Cb + 2*G[i]* Vx[i]*abs(yv[i]*dy[i]+zv[i]*dz[i])/(bref*Sref)
            CDi = CDi- 2*G[i]*(wi[i]*dy[i]-vi[i]*dz[i])/Sref

            CL_a = CL_a + 2*(G_a[i]*Vx[i] + G[i]*Vx_a) * dy[i]/Sref
            Cr_a = (Cr_a - 2*(G_a[i]*Vx[i] + G[i]*Vx_a)
                    * (yv[i]*dy[i]+zv[i]*dz[i])/(bref*Sref))
            Cn_a = (Cn_a - 2*(G_a[i]*(Vz+wi[i]) + G[i]*(Vz_a+wi_a[i]))
                    * yv[i]*dy[i]/(bref*Sref))

            CL_b = CL_b + 2*(G_b[i]*Vx[i] + G[i]*Vx_b) * dy[i]/Sref
            Cr_b = (Cr_b - 2*(G_b[i]*Vx[i] + G[i]*Vx_b)
                    * (yv[i]*dy[i]+zv[i]*dz[i])/(bref*Sref))
            Cn_b = (Cn_b - 2*(G_b[i]*(Vz+wi[i]) + G[i]*(Vz_b+wi_b[i]))
                    * yv[i]*dy[i]/(bref*Sref))

            CL_p = CL_p + 2*(G_p[i]*Vx[i] + G[i]*Vx_p) * dy[i]/Sref
            Cr_p = (Cr_p - 2*(G_p[i]*Vx[i] + G[i]*Vx_p)
                    * (yv[i]*dy[i]+zv[i]*dz[i])/(bref*Sref))
            Cn_p = (Cn_p - 2*(G_p[i]*(Vz+wi[i]) + G[i]*(Vz_p+wi_p[i]))
                    * yv[i]*dy[i]/(bref*Sref))

            CL_r = CL_r + 2*(G_r[i]*Vx[i] + G[i]*Vx_r) * dy[i]/Sref
            Cr_r = (Cr_r - 2*(G_r[i]*Vx[i] + G[i]*Vx_r)
                    * (yv[i]*dy[i]+zv[i]*dz[i])/(bref*Sref))
            Cn_r = (Cn_r - 2*(G_r[i]*(Vz+wi[i]) + G[i]*(Vz_r+wi_r[i]))
                    * yv[i]*dy[i]/(bref*Sref))

        nsys = npar
        Asys = zeros([nsys, nsys])
        Rsys = zeros(nsys)

        ksys = -1
        if ialspec == 1:
            ksys = ksys+1
            Rsys[ksys] = alpha - alspec
            Asys[ksys, 0] = 1
        if ibespec == 1:
            ksys = ksys+1
            Rsys[ksys] = beta - bespec
            Asys[ksys, 1] = 1
        if ipbspec == 1:
            ksys = ksys+1
            Rsys[ksys] = pbar - pbspec
            Asys[ksys, 2] = 1
        if irbspec == 1:
            ksys = ksys+1
            Rsys[ksys] = rbar - rbspec
            Asys[ksys, 3] = 1
        if iCLspec == 1:
            ksys = ksys+1
            Rsys[ksys] = CL - CLspec
            Asys[ksys, 0] = CL_a
            Asys[ksys, 1] = CL_b
            Asys[ksys, 2] = CL_p
            Asys[ksys, 3] = CL_r
        if iCrspec == 1:
            ksys = ksys+1
            Rsys[ksys] = Cr - Crspec
            Asys[ksys, 0] = Cr_a
            Asys[ksys, 1] = Cr_b
            Asys[ksys, 2] = Cr_p
            Asys[ksys, 3] = Cr_r
        if iCnspec == 1:
            ksys = ksys+1
            Rsys[ksys] = Cn - Cnspec
            Asys[ksys, 0] = Cn_a
            Asys[ksys, 1] = Cn_b
            Asys[ksys, 2] = Cn_p
            Asys[ksys, 3] = Cn_r

        if (ksys+1 != nsys):
            print ('Error: %3i quantities are specified for flow '
                   'condition %3i .  Must have 4. \n' % (ksys, 1))

        dpar = lstsq(-Asys, Rsys)[0]

        dalpha = dpar[0]
        dbeta = dpar[1]
        dpbar = dpar[2]
        drbar = dpar[3]

        alpha = alpha+ dalpha
        beta = beta + dbeta
        pbar = pbar + dpbar
        rbar = rbar + drbar

        itr = itr + 1

    return A, wzG, dy, twa, Vx, G, yv,zv,cl,ccl,vi,wi,alpha,beta,pbar,rbar,CL,CDi,Cr,Cn,Cb
