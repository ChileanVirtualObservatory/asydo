# Create several cubes in parallel, using random parameters
from asydopy import *
import matplotlib.pyplot as plt
from scipy import signal
import time
import math
import pickle
import numpy as np
from matplotlib import animation

dbpath = "../../ASYDO"
template = factory.IMCConf(0, dbpath,
                           mol_list="all",
                           mol_prob=0.3,
                           x_pos=0.0,
                           y_pos=0.0,
                           f_pos=300000,
                           spa_pix=32,
                           spe_pix=500,
                           fov=500,
                           bw=2000,
                           rvel=(150, 1000),
                           temp=(50, 500),
                           semiaxis=(10, 300),
                           fwhm=(10, 50),
                           angle=(0, math.pi),
                           rot=(50, 500),
                           curtosis=(-5, 5))

inte = 100
cube = pickle.loads(factory.unitary_IMC(template))
x = 0
# fig=plt.figure()
fig, axarr = plt.subplots(3, sharex=True)
cwtmatr = []
specm = []
mmin = float('Inf')
mmax = float('-Inf')
img = []
img2 = []
img3 = []
for x in range(len(cube.alpha_axis)):
    spec = cube.data[:, 0, x]
    img2.append(axarr[1].plot(spec, 'b'))
    wavelet = signal.ricker
    widths = np.arange(2, 5, 0.2)
    cmat = signal.cwt(spec, wavelet, widths)
    for y in range(len(widths)):
      img3.append(axarr[2].plot(cmat[y]))
    cmax = cmat.max()
    if cmax > mmax:
        mmax = cmax
    cmin = cmat.min()
    if cmin < mmin:
        mmin = cmin
    cwtmatr.append(cmat)
for cmat in cwtmatr:
    img.append([axarr[0].imshow(cmat, vmin=mmin, vmax=mmax)])

ani = animation.ArtistAnimation(fig, img, interval=inte, blit=True,
                                repeat=True)
ani = animation.ArtistAnimation(fig, img2, interval=inte, blit=True,
                                repeat=True)
ani = animation.ArtistAnimation(fig, img3, interval=inte/len(widths), blit=True,
                                repeat=True)
plt.show()

