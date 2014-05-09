from asydoCore import *
import sys

rvel=150.0
temp=300.0
molist=('COv=0','13COv=0','C18O','C17O','13C18O','NH2','N2H+v=0','CNv=0','HCNv=0','HNCv=0','H2CN','CSv=0','CCS','H2S','H2CS','SO2v=0','H2CO','HCO+v=0','HC3Nv=0','HC5Nv=0','CH3OHvt=0')


log=open('script-alone.log', 'w')

counter=0.0
for mol in molist:
   univ=SynUniverse(log)
   univ.addSource(log,'alone-'+mol,0.0,counter,rvel,temp)
   s_x=random.uniform(0.0001, 0.001)
   s_y=random.uniform(0.0001, 0.001)
   rot=random.uniform(0.0, 1.67)
   mov=random.uniform(0.05,0.1)
   univ.addStruct(log,'alone-'+mol,mol,('Gaussian',s_x,s_y,rot,1.0),('Gaussian',mov,1.0))
   for i in range(49):
      fcenter=276 + 2*i
      cube=univ.genCube(log,'alone-'+mol,0.0,counter,float(fcenter),0.1,0.8,2000,2000000,'alone-'+mol+'-'+str(fcenter)+'.fits')
   counter+=1

