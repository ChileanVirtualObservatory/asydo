from SynUniverse import *
import sys

rvel=150.0
temp=300.0
molist=('CO','13CO','C18O','C17O','13C18O','NH3','N2H+','CN','HCN','HNC','CH3CN','H2CN','CS','CCS','H2S','H2CS','SO','SO2','C3H2','H2CO','HCO+','HC3N','HC5N','HC7N','HC9N','HC11N','CH3OH')


log=open('PUC-script.log', 'w')

univ=SynUniverse(log)
for mol in molist:
   univ.addSource(log,'alone-'+mol,0.0,float(i),rvel,temp)
   s_x=random.uniform(0.01, 0.1)
   s_y=random.uniform(0.01, 0.1)
   rot=random.uniform(0.0, 1.67)
   mov=random.uniform(0.01,0.05)
   univ.addStruct(log,'alone-'+mol,mol,('Gaussian',s_x,s_y,rot,1.0),('Gaussian',mov,1.0))
   for i in range(98):
      fcenter=276 + i
		cube=univ.genCube(log,'alone-'+mol,0.0,float(i),0.1,0.8,float(fcenter),2000,2000000,'alone-'+mol+'-'+fcenter+'.fits')

