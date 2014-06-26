from asydoCore import *
import matplotlib.pyplot as plt
import sys

log = sys.stdout

univ=Universe(log)
univ.addSource(log,'P-33SO2',0.0,1.0,150,300.0)
univ.addStruct(log,'P-33SO2','33SO2,COv=0',('Gaussian',0.1,0.2,0.6,1.0),('Gaussian',0.4,1.0))
cube=univ.genCube(log,'obs1',0.0,1.0,324.0,0.1,0.8,24500,98000000,'p33SO2-obs1.fits')
plt.plot(cube.getSpectrum(0.0,1.0))
plt.show()
#cube.animate(1)

