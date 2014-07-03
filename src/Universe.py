from Source import *
from Cube import *


import numpy as np

from collections import namedtuple

import string
from astropy.io import fits
from scipy import signal
import matplotlib.animation as animation
import matplotlib.pyplot as plt


UniverseSpec = namedtuple('SynUniverse','profile band_freq band_noise inten_group inten_values iso_abun base_abun base_CO')
CubeSpec = namedtuple('CubeSpec', 'alpha delta freq ang_res ang_fov spe_res spe_bw')


defaultUniverse = UniverseSpec('default', \
                          {'3': [88, 116], '4': [125, 163], '6': [211, 275], '7': [275, 373], '8': [385, 500],'9': [602, 720]},\
                          {'3': 0.01, '4': 0.012, '6': 0.02, '7': 0.04, '8': 0.08, '9': 0.16},\
                          [('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')],\
                          [[0.1, 2], [20, 60], [5, 20], [1, 10]],\
                          {'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': 1.0/120, '34S': 1.0 / 30, '33S': 1.0 / 120,'13N': 1.0 / 30, 'D': 1.0 / 30},\
                          [10 ** -5, 10 ** -6],\
                          1.0)

class Universe:
    """Synthetic Cube Generation Class"""

    def __init__(self, log):
        self.conf = defaultUniverse
        self.log=log
        self.sources = dict()

    def setSpec(self, spec):
        self.conf = spec

    def createSource(self, name, alpha , delta):
        self.sources[name] = Source(self.log, name, alpha, delta)

    def addComponent(self, source_name, model):
        self.sources[source_name].addComponent(model)

    def genCube(self, name, cubespec,filename) 
        cube = Cube(self.log, name, self.conf, cubespec, filename)
        for src in self.sources: 
            self.log.write('   * Source: ' + src + '\n')
            self.sources[src].emission(cube, self.conf.inten_group, self.conf.inten_values)
        self.log.write('   * Adding Noise... \n')
        cube.saveFits(self.sources, filename)
        return cube

    def removeSource(self, name):
        self.log.write('Removing source '+source)
        return sources.remove(name)

