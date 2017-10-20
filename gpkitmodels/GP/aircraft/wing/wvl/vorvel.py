import numpy as np
from numpy.matlib import repmat

def vorvel(ra, rb, r, ibound):

    r = r.transpose()
    N = r.shape[0]
    a = r - repmat(ra, N, 1)
    b = r - repmat(rb, N, 1)
    ih = repmat(np.array([1, 0, 0]), N, 1)

    am = np.array([np.sqrt(np.dot(ai, ai)) for ai in a])
    bm = np.array([np.sqrt(np.dot(bi, bi)) for bi in b])
    adb = np.array([np.dot(ai, bi) for ai, bi in zip(a, b)])
    axb = np.cross(a, b)
    axi = np.cross(a, ih)
    bxi = np.cross(b, ih)

    den = am*bm + adb
    dena = am - a[:, 0]
    denb = bm - b[:, 0]

    v = (repmat((1/am)/dena, 3, 1).transpose()*axi
         - repmat((1/bm)/denb, 3, 1).transpose()*bxi)

    if ibound == 1:
        v += repmat((1/am + 1/bm)/den, 3, 1).transpose()*axb

    vel = v/4/np.pi
    vel = vel.transpose()
    return vel
