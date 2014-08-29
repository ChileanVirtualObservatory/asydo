# Create several cubes in parallel, using random parameters
import matplotlib.pyplot as plt
from CubeFactory import *

global CFParams

sample_size=30000
CFParams.ban_list.append('Phosphapropynylidyne')
negative=np.asarray(parallelGen(sample_size/2,4))
CFParams.ban_list=list()
CFParams.force_list.append('Phosphapropynylidyne')
positive=np.asarray(parallelGen(sample_size/2,4))
result=np.vstack((negative,positive))
np.save('exp1.npy', result)


#import cProfile
#import re
#cProfile.run('unitaryGen(0)')
