from multiprocessing import Pool
from Universe import *
from IMCM import *
import random
import numpy as np
import sys
import time
import matplotlib.pyplot as plt


log=open('script-combined.log', 'w')

molist=('COv=0','13COv=0','C18O','C17O','13C18O','NH2','N2H+v=0','CNv=0','HCNv=0','HNCv=0','H2CN','CSv=0','CCS','H2S','H2CS','SO2v=0','H2CO','HCO+v=0','HC3Nv=0','HC5Nv=0','CH3OHvt=0')
siz=len(molist)


target=3
samples=100
spa_pix=5
spe_pix=100
#vect=np.zeros((samples,spa_pix*spa_pix*spe_pix))

def gen_cube(x):
   global vect
   print x
   mask = np.random.randint(2,size=siz)
   #print mask
   if x > samples/2:
      mask[target]=False
   else:
      mask[target]=True
   univ=Universe(log)
   for j in range(siz):
      if mask[j]==0:
          continue
      mol=molist[j]
      univ.createSource('combined-'+mol,0.0,0.0)
      rvel=random.uniform(1, 200)
      temp=random.uniform(100, 400)
      s_x=random.uniform(1, 100)
      s_y=random.uniform(1, 100)
      rot=random.uniform(1, 100)
      s_f=random.uniform(40, 400)
      angle=random.uniform(0,math.pi)
      model=IMCM(log,mol,temp,('normal',s_x,s_y,angle),('skew',s_f,randint(-7,7)),('linear',angle,rot))
      model.setRadialVelocity(rvel)
      univ.addComponent('combined-'+mol,model)
   fcenter=300000
   cspec=CubeSpec(0.0,0.0,fcenter,500/spa_pix,500,20000/spe_pix,20000)
   cube=univ.genCube('combined',cspec)
   return cube.data.flatten()


p = Pool(8)
result=p.map(gen_cube,range(samples))
print result
info = np.asarray(result)
np.save('exp1.npy', info)
#np.savetxt('exp1.txt', vect, newline="\n")
  
#   cube.animate(1,False)
#print vect
#plt.plot(vect.transpose())
#plt.show()

