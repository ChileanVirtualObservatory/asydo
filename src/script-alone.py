from Universe import *
from IMCM import *
import random
import math
import sys

rvel=150.0
temp=300.0
molist=('COv=0','13COv=0','C18O','C17O','13C18O','NH2','N2H+v=0','CNv=0','HCNv=0','HNCv=0','H2CN','CSv=0','CCS','H2S','H2CS','SO2v=0','H2CO','HCO+v=0','HC3Nv=0','HC5Nv=0','CH3OHvt=0')

log=open('script-alone.log', 'w')

counter=0.0
univ=Universe(log)
for mol in molist:
   univ.createSource('alone-'+mol,0.0,counter)
   s_x=random.uniform(50, 150)
   s_y=random.uniform(40, 100)
   rot=random.uniform(10, 150)
   s_f=random.uniform(50, 120)
   angle=random.uniform(0,math.pi)
   model=IMCM(log,mol,temp,('normal',s_x,s_y,angle),('skew',s_f,0),('linear',angle,rot))
   model.setRadialVelocity(rvel)
   univ.addComponent('alone-'+mol,model)
   for i in range(49):
      fcenter=(276 + 2*i)*1000
      cspec=CubeSpec(0.0,counter,fcenter,10,800,2,2000)
      cube=univ.genCube('alone-'+mol,cspec)
      univ.saveCube(sube, 'alone-'+mol+'-'+str(fcenter)+'.fits')
   counter+=1

