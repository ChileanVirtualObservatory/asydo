from Source import *
from Cube import *


import numpy as np

from collections import namedtuple

import string
from astropy.io import fits
from scipy import signal
import matplotlib.animation as animation
import matplotlib.pyplot as plt


SynConf = namedtuple('SynUniverse','profile band_freq band_noise inten_group inten_values iso_abun base_abun base_CO')

CubeSpec = namedtuple('CubeSpec', 'alpha delta freq ang_res ang_fov spe_res spe_bw')

defaultUniverse = SynConf('default', \
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
        self.sources = dict()

    def setConf(self, log, univ):
        self.conf = univ

    def addSource(self, log, name, alpha , delta, rad_vel, temp):
        self.sources[name] = Source(log, name, alpha, delta, rad_vel, temp)

    def addStruct(self, log, name, mol_list, spa_form, spe_form):
        self.sources[name].addStruct(log, mol_list, spa_form, spe_form, self.conf)

    def genCube(self, log, name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw, filename):
        cube = Cube(log, name, self.conf, CubeSpec(alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw))
        for src in self.sources: 
            log.write('   * Source: ' + src + '\n')
            self.sources[src].emission(log, cube, self.conf.inten_group, self.conf.inten_values)
        log.write('   * Adding Noise... \n')
        cube.saveFits(self.sources, filename)
        return cube

    def removeSource(self, log, name):
        self.sources.remove(name)
        log.write('Removing Source ' + name + '\n')

    def genImage(self, log, name, alpha, delta, ang_res, ang_fov, filename):
        # TODO: must be implemented by the IS group
        return 0




