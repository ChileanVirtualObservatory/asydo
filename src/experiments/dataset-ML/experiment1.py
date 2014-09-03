# Create several cubes in parallel, using random parameters
import matplotlib.pyplot as plt
from CubeFactory import *

global CFParams

sample_size=30
cores=8
CFParams.ban_list.append('Phosphapropynylidyne')
cubes=parallelGen(sample_size/2,cores)
data=[]
for cube in cubes:
   data.append(cube.data.flatten())
negative=np.asarray(data)
CFParams.ban_list=list()
CFParams.force_list.append('Phosphapropynylidyne')
cubes=parallelGen(sample_size/2,cores)
data=[]
for cube in cubes:
   data.append(cube.data.flatten())
positive=np.asarray(data)
result=np.vstack((negative,positive))
np.save('test.npy', result)


#import cProfile
#import re
#cProfile.run('unitaryGen(0)')
