# Create several cubes in parallel, using random parameters
from CubeFactory import *

global CFParams

sample_size=10000
CFParams.ban_list.append('Phosphapropynylidyne')
CFParams.spe_pix=1000
CFParams.spa_pix=5
negative=np.asarray(parallelGen(sample_size/2,4))
CFParams.ban_list=list()
CFParams.force_list.append('Phosphapropynylidyne')
positive=np.asarray(parallelGen(sample_size/2,4))
result=np.vstack((negative,positive))
np.save('exp2.npy', result)


#import cProfile
#import re
#cProfile.run('unitaryGen(0)')
