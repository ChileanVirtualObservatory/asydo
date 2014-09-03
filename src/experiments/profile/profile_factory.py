# Create several cubes in parallel, using random parameters
import matplotlib.pyplot as plt
from asydopy import *
import math
import numpy as np
import cProfile
import re

dbpath="../../ASYDO"
template=factory.IMCConf(0,dbpath)
cProfile.run('factory.unitary_IMC(template)')

