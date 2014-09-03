from asydopy import *
import matplotlib.pyplot as plt
import sys
import math

log = sys.stdout
dbpath = 'ASYDO'

univ=vu.Universe(log)
univ.create_source('example',0.0,1.0)
univ.create_source('example2',0.05,0.95)
model=vu.IMCM(log,dbpath,'33SO2',300,('normal',100,70,0),('skew',100,3),('linear',math.pi/9,500.0))
model2=vu.IMCM(log,dbpath,'33SO2',300,('exp',100,70,math.pi/9),('skew',100,0),('linear',math.pi/4,100.0))
model.set_radial_velocity(150)
model2.set_radial_velocity(100)
univ.add_component('example',model)
univ.add_component('example2',model2)
cube=univ.gen_cube('obs1',0.0,1.0,324000.0,10,800,2.45,3000)
univ.save_cube(cube,'p33SO2-obs1.fits')
plt.plot(cube.get_spectrum(0.0,1.0))
plt.show()
cube.animate(10,True)

