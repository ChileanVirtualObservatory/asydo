from SynUniverse import *
import sys

log=sys.stdout

univ=SynUniverse(log)
univ.addSource(log,'P-33SO2',0.0,1.0,150,300.0)
univ.addStruct(log,'P-33SO2','(33)SO2,CO',('Gaussian',0.1,0.2,0.6,1.0),('Gaussian',0.1,1.0))
cube=univ.genCube(log,'obs1',0.0,1.0,0.1,0.8,324.0,24500,98000000,'p33SO2-obs1.fits')
plot(cube.getSpectrum(0.0,1.0))
show()
#cube.animate(0.01)

