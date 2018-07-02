" Weissenger vortex lattice method "
from numpy import array, pi, zeros, flip, diff

def wvl(self, c):

    N = self.Nwvl

    y = c[self.eta]
    yv = y[:-1] + diff(y)/2
    y = array(list(flip(-y, 0)) + list(y[1:]))
    yv = array(list(flip(-yv, 0)) + list(yv))

    z = zeros(len(y))
    zv = zeros(len(yv))

    wyG = zeros([N, N])
    wzG = zeros([N, N])

    for i in range(N):
        rsqi = (yv-y[i])**2 + (zv-z[i])**2
        rsqp = (yv-y[i+1])**2 + (zv-z[i+1])**2
        wyG[i, :] = ((zv-z[i])/rsqi - (zv-z[i+1])/rsqp)/(4*pi)
        wzG[i, :] = -((yv-y[i])/rsqi - (yv-y[i+1])/rsqp)/(4*pi)

    dy = diff(y)
    dz = diff(z)

    B = - wzG*dy + wyG*dz

    return B
