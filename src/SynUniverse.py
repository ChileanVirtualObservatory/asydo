import numpy as np
from pylab import *
from collections import namedtuple
import random
import string
from astropy.io import fits
from astropy.io.votable import parse_single_table
from scipy import signal
import matplotlib.animation as animation
import matplotlib.pyplot as plt


SynConf = namedtuple('SynUniverse',
                     'profile band_freq band_noise inten_group inten_values iso_abun base_abun base_CO')

SynStruct = namedtuple('SynStruct', 'code intens spa_form spe_form')

CubeSpec = namedtuple('CubeSpec', 'x_center y_center ang_res ang_fov v_center spe_res spe_bw')

defaultUniverse = SynConf('default', \
                          {'3': [88, 116], '4': [125, 163], '6': [211, 275], '7': [275, 373], '8': [385, 500],
                           '9': [602, 720]}, \
                          {'3': 0.1, '4': 0.12, '6': 0.2, '7': 0.4, '8': 0.8, '9': 1.6}, \
                          [('default'), ('CO'), ('13CO'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')], \
                          [[0.1, 2], [20, 60], [5, 20], [1, 10]], \
                          {'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': '1.0/120', '34S': 1.0 / 30, '33S': 1.0 / 120,
                           '13N': 1.0 / 30, 'D': 1.0 / 30}, \
                          [10 ** -5, 10 ** -6], \
                          1.0)

DEG2ARCSEC = 3600.0
S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2)
MAX_CHANNELS = 9000
MAX_BW = 2000000.0 # kHz
SPEED_OF_LIGHT = 299792458.0


class SynCube:
    """ A synthetic ALMA cube."""


    def __init__(self, log, name, conf, spec):
        """ Parameters:
                 log	: descriptor of a log file
                 name : name of the cube
                 conf : the SynConf configuration
                 spec	: a CubeSpec specification
            """
        self.name = name
        self.spec = spec
        self.conf = conf
        log.write('Generating cube ' + name + '\n')
        log.write('  -> Spatial Center (deg): ra=' + str(spec.x_center) + ' dec=' + str(spec.y_center) + '\n')
        fact = spec.ang_fov / (DEG2ARCSEC* spec.ang_res)
        self.x_border = [spec.x_center - fact / 2, spec.x_center + fact / 2]
        self.y_border = [spec.y_center - fact / 2, spec.y_center + fact / 2]
        self.x_axis = np.linspace(self.x_border[0], self.x_border[1], spec.ang_fov / spec.ang_res)
        self.y_axis = np.linspace(self.y_border[0], self.y_border[1], spec.ang_fov / spec.ang_res)
        if spec.x_center > 90 or spec.x_center < -90:
            raise Exception('ERROR: invalid coordinate: ra=' + spec.x_center)
        if spec.y_center > 90 or spec.y_center < -90:
            raise Exception('ERROR: invalid coordinate: dec=' + spec.y_center)
        log.write('  -> FOV (arcsec): ra=' + str(self.x_border) + ' dec=' + str(self.y_border) + '\n')
        self.v_border = [spec.v_center - spec.spe_bw / (2.0 * 1000000), spec.v_center + spec.spe_bw / (2.0 * 1000000)]
        if spec.spe_bw > MAX_BW:
            log.write('WARNING: max ALMA bandwidth exceeded\n')
        self.channels = round(spec.spe_bw / spec.spe_res)
        if self.channels > MAX_CHANNELS:
            log.write('WARNING: max ALMA channels exceeded\n')
        self.v_axis = np.linspace(self.v_border[0], self.v_border[1], self.channels)
        log.write('  -> Spectral (GHz): center=' + str(spec.v_center) + ' bandwidth=' + str(self.v_border) + '\n')
        log.write('  -> Cube size: ' + str(len(self.x_axis)) + ' x ' + str(len(self.y_axis)) + ' x ' + str(
            len(self.v_axis)) + ' \n')
        self.band = 'NO_BAND'
        for bnd in conf.band_freq:
            freqs = conf.band_freq[bnd]
            if self.v_border[0] >= freqs[0] and self.v_border[1] <= freqs[1]:
                self.band = bnd
                log.write('  -> Band: ' + bnd + '\n')
        if self.band == 'NO_BAND':
            log.write('WARNING: not in a valid ALMA band\n')
        self.data = zeros((len(self.x_axis), len(self.y_axis), len(self.v_axis)))
        self.hdu = fits.PrimaryHDU()


    def freqWindow(self, freq, fwhm):
        """ Frequency window.
                 Given a freq and a fwhm compute returns the range of channels that are affectet
            """
        factor = 2.0;
        ldiff = freq - self.v_border[0] - factor * fwhm;
        udiff = freq - self.v_border[0] + factor * fwhm;
        l_w = 0
        if (ldiff > 0):
            l_w = ldiff * 1000000 / self.spec.spe_res
        l_u = udiff * 1000000 / self.spec.spe_res
        if l_u > self.channels:
            l_u = channels - 1
        return (int(l_w), int(l_u))


    def _updatefig(self, j):
        """ Animate helper function """
        self.im.set_array(self.data[:, :, j])
        return self.im,


    def animate(self, inte, rep=True):
        """ Animate the cube.
                inte       : time interval between frames
                 rep[=True] : boolean to repeat the animation
          """
        fig = plt.figure()
        self.im = plt.imshow(self.data[:, :, 0], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
                             extent=(self.x_border[0], self.x_border[1], self.y_border[0], self.y_border[1]))
        ani = animation.FuncAnimation(fig, self._updatefig, frames=range(len(self.v_axis)), interval=inte, blit=True,
                                      repeat=rep)
        plt.show()


    def addNoise(self):
        """ Adds the noise. The noise of the level is defined in the Universe configuration (band_noise) """
        if self.band == 'NO_BAND':
            noise = 0.1
        else:
            noise = self.conf.band_noise[self.band]
        self.data += np.random.random((len(self.x_axis), len(self.y_axis), len(self.v_axis))) * noise


    def getSpectrum(self, x, y):
        """ Return an spectrum in x, y """
        xi = int(round((x - self.x_axis[0]) / self.spec.ang_res))
        yi = int(round((y - self.y_axis[0]) / self.spec.ang_res))
        return self.data[xi][yi]


    def saveFits(self, filename):
        """ Write the final FITS file in filename """
        self.hdu.data = self.data
        self.hdu.writeto(filename, clobber=True)


class SynSource:
    """A source"""

    def __init__(self, log, name, pos, rad_vel, temp):
        log.write('Source \'' + name + '\' added \n')
        self.pos = pos
        self.name = name
        self.rad_vel = rad_vel
        self.temp = temp
        self.structs = list()

    def addStruct(self, log, mol_list, spa_form, spe_form, conf):
        code = self.name + '-' + str(len(self.structs) + 1)
        intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(conf.base_abun[0], conf.base_abun[1])
            if mol in ('CO', '13CO', 'C18O', 'C17O', '13C18O'):
                abun += conf.base_CO
            for iso in conf.iso_abun:
                if iso in mol:
                    abun *= conf.iso_abun[iso]
            intens[mol] = abun
        self.structs.append(SynStruct(code, intens, spa_form, spe_form))
        log.write('Added to \'' + self.name + '\': molecules ' + str(mol_list) + '\n')
        log.write('  -> Spatial form  = ' + str(spa_form) + '\n')
        log.write('  -> Spectral form = ' + str(spe_form) + '\n')

    def genSurface(self, form, cube):
        "Create a gaussian surface over a mesh created by x and y axes"
        sx = form[1]
        sy = form[2]
        theta = form[3]
        r = 3 * sqrt(sx ** 2 + sy ** 2)
        x_mesh, y_mesh = meshgrid(cube.x_axis, cube.y_axis, sparse=False, indexing='xy')
        Xc = x_mesh.flatten() - self.pos[0] * ones(len(cube.x_axis) * len(cube.y_axis))
        Yc = y_mesh.flatten() - self.pos[1] * ones(len(cube.x_axis) * len(cube.y_axis))
        XX = (Xc) * cos(theta) - (Yc) * sin(theta);
        YY = (Xc) * sin(theta) + (Yc) * cos(theta);
        u = (XX / sx) ** 2 + (YY / sy) ** 2;
        sol = sx * sy * exp(-u / 2) / (2 * pi);
        res = transpose(reshape(sol, (len(cube.y_axis), len(cube.y_axis))))
        return res

    def emission(self, log, cube, inten_group, inten_values):
        log.write('Loading visible lines in band ' + cube.band + ' (rad_vel=' + str(self.rad_vel) + ')\n')
        lines = self.loadLines(cube.band, cube.v_border[0], cube.v_border[1], self.rad_vel)
        for struct in self.structs:
            log.write(' --> Struct Name: ' + struct.code + '\n')
            tcub = self.genSurface(struct.spa_form, cube)
            fwhm = struct.spe_form[1]
            shape = struct.spe_form[2]
            for mol in struct.intens:
                mlin = lines.where((lines.molformula == mol))
                rinte = inten_values[0]
                for j in range(len(inten_group)):
                    if mol in inten_group[j]:
                        rinte = inten_values[j]
                rinte = random.uniform(rinte[0], rinte[1])
                for i in range(len(mlin)):
                    try:
                        trans_temp = float(mlin[i]['lowerstateenergyK'])
                    except ValueError:
                        log.write('WARNING: Strange transition temperature, printing line...\n')
                        log.write('WARNING: ' + str(mlin[i]) + '\n')
                        log.write('WARNING: Setting transition temperature to ' + str(self.temp) + '\n')
                        trans_temp = self.temp
                    temp = exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                    log.write('E: ' + mlin[i]['chemicalname'] + ' (' + mlin[i]['molformula'] + ') at ' + str(
                        mlin[i]['frequency']) + ' Mhz, T = exp(-|' \
                              + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(
                        rinte) + ' = ' + str(temp) + ' K\n')
                    freq = mlin[i]['frequency']
                    window = cube.freqWindow(freq / 1000.0, fwhm)
                    log.write('Window:' + str(window) + '\n')
                    sigma = fwhm / S_FACTOR 
                    for idx in range(window[0], window[1]):
                        cube.data[:, :, idx] += tcub * temp * exp(
                            (-0.5 * (cube.v_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape)))
                    #TODO: FIX THIS PROBLEM NOW!!!
                    #print temp,cube.v_axis,freq/1000.0,sigma,shape

    def loadLines(self, band, v_init, v_end, rad_vel):
        # TODO: Read from a database using SQLINE (SS Group)
        v_init_corr=(1 + rad_vel*1000.0/SPEED_OF_LIGHT)*v_init
        v_end_corr=(1 + rad_vel*1000.0/SPEED_OF_LIGHT)*v_end
        location = './votables/band' + band + '.xml'
        tbl = parse_single_table(location)
        #for i in 
        #tbl.where((tbl.frequency >= v_init_corr) & (tbl.frequency <= v_end_corr))
        return tbl


class SynUniverse:
    """Synthetic Cube Generation Class"""

    def __init__(self, log):
        self.conf = defaultUniverse
        self.sources = dict()

    def setConf(self, log, univ):
        self.conf = univ

    def addSource(self, log, name, x_pos, y_pos, rad_vel, temp):
        self.sources[name] = SynSource(log, name, [x_pos, y_pos], rad_vel, temp)

    def addStruct(self, log, name, mol_list, spa_form, spe_form):
        self.sources[name].addStruct(log, mol_list, spa_form, spe_form, self.conf)

    def genCube(self, log, name, x_center, y_center, ang_res, ang_fov, v_center, spe_res, spe_bw, filename):
        cube = SynCube(log, name, self.conf, CubeSpec(x_center, y_center, ang_res, ang_fov, v_center, spe_res, spe_bw))
        for src in self.sources:
            log.write('   * Source: ' + src + '\n')
            self.sources[src].emission(log, cube, self.conf.inten_group, self.conf.inten_values)
        #log.write('   * Adding Noise... \n')
        #cube.addNoise(self.conf.band_noise) # To uncomment by SS group
        cube.saveFits(filename)
        return cube

    def removeSource(self, log, name):
        #TODO: must be implemented by VU group
        return 0

    def genImage(self, log, name, x_center, y_center, ang_res, ang_fov, filename):
        #TODO: must be implemented by the IS group
        return 0


