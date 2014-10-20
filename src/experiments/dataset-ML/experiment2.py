# Create several cubes in parallel, using random parameters
import matplotlib.pyplot as plt
from asydopy import *
import math
import numpy as np

sample_size=10000
dbpath="../../ASYDO"
template=factory.IMCConf(0,dbpath,
                mol_list="all",
                mol_prob=0.3,
                x_pos=0.0,
                y_pos=0.0,
                f_pos=300000,
                spa_pix=5,
                spe_pix=100,
                fov=500,
                bw=2000,
                rvel=(150,1000),
                temp=(50,500),
                semiaxis=(10,300),
                fwhm=(10,50),
                angle=(0,math.pi),
                rot=(50,500),
                curtosis=(-5,5))

template.ban_list.append('Phosphapropynylidyne')
conf_list=list()

for i in range(0,sample_size/2):
    conf=factory.IMCConf(i)
    conf.set_params(template)
    conf_list.append(conf)
cubes=factory.gen_IMC_cubes(conf_list)
data=[]
for cube in cubes:
   data.append(cube.data.flatten())
negative=np.asarray(data)
template.ban_list=list()
template.force_list.append('Phosphapropynylidyne')
conf_list=list()
for i in range(sample_size/2,sample_size):
    conf=factory.IMCConf(i)
    conf.set_params(template)
    conf_list.append(conf)
cubes=factory.gen_IMC_cubes(conf_list)
data=[]
for cube in cubes:
   data.append(cube.data.flatten())
positive=np.asarray(data)
result=np.vstack((negative,positive))
np.save('exp2.npy', result)


#import cProfile
#import re
#cProfile.run('unitaryGen(0)')
