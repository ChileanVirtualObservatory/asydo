from Universe import *
from IMCM import *
import matplotlib.pyplot as plt
import sys

log = sys.stdout

cspec=CubeSpec(0.0,1.0,324000.0,10,800,2.45,3000)

univ=Universe(log)
univ.createSource('example',0.0,1.0)
univ.createSource('example2',0.05,0.95)
model=IMCM(log,'33SO2',300,('normal',100,70,0),('skew',100,3),('linear',math.pi/9,500.0))
model2=IMCM(log,'33SO2',300,('exponential',100,70,math.pi/9),('skew',100,0),('linear',math.pi/4,100.0))
model.setRadialVelocity(150)
model2.setRadialVelocity(100)
univ.addComponent('example',model)
univ.addComponent('example2',model2)
cube=univ.genCube('obs1',cspec)
univ.saveCube(cube,'p33SO2-obs1.fits')
plt.plot(cube.getSpectrum(0.0,1.0))
plt.show()
cube.animate(10,False)

