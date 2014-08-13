from Universe import *
from IMCM import *
import matplotlib.pyplot as plt
import sys

log = sys.stdout

cspec=CubeSpec(0.0,1.0,324000.0,10,800,24.5,98000)

univ=Universe(log)
univ.createSource('example',0.0,1.0)
model=IMCM(log,'33SO2,COv=0',300,('Normal',0.01,0.02,0.6,1.0),('Normal',500,3),('Linear',math.pi/9,3))
model.setRadialVelocity(150)
univ.addComponent('example',model)
cube=univ.genCube('obs1',cspec, 'p33SO2-obs1.fits')
plt.plot(cube.getSpectrum(0.0,1.0))
plt.show()
cube.animate(1,False)

