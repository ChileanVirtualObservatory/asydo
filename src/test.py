from Universe import *
from IMCM import *
import matplotlib.pyplot as plt
import sys

log = sys.stdout

cspec=CubeSpec(0.0,1.0,324000.0,10,800,2.45,3000)

univ=Universe(log)
univ.createSource('example',0.0,1.0)
model=IMCM(log,'33SO2,COv=0',300,('Normal',100,70,math.pi/9),('Normal',100,3),('Linear',math.pi/9,500.0))
model.setRadialVelocity(150)
univ.addComponent('example',model)
cube=univ.genCube('obs1',cspec, 'p33SO2-obs1.fits')
plt.plot(cube.getSpectrum(0.0,1.0))
plt.show()
cube.animate(10,False)

