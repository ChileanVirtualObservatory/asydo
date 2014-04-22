import numpy as np
from database import DataBase
from pylab import *
from collections import namedtuple
import random
import string
from astropy.io import fits
from astropy.io.votable import parse_single_table
from astropy.table import Table
from astropy.io.votable.tree import Table as pTable
from astropy.io.votable.tree import Field as pField
import atpy
from scipy import signal
import matplotlib.animation as animation
import matplotlib.pyplot as plt


SynConf = namedtuple('SynUniverse','profile band_freq band_noise inten_group inten_values iso_abun base_abun base_CO')

CubeSpec = namedtuple('CubeSpec', 'x_center y_center ang_res ang_fov v_center spe_res spe_bw')

defaultUniverse = SynConf('default', \
                          {'3': [88, 116], '4': [125, 163], '6': [211, 275], '7': [275, 373], '8': [385, 500],'9': [602, 720]},\
                          {'3': 0.01, '4': 0.012, '6': 0.02, '7': 0.04, '8': 0.08, '9': 0.16},\
                          [('default'), ('CO'), ('13CO'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')],\
                          [[0.1, 2], [20, 60], [5, 20], [1, 10]],\
                          {'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': '1.0/120', '34S': 1.0 / 30, '33S': 1.0 / 120,'13N': 1.0 / 30, 'D': 1.0 / 30},\
                          [10 ** -5, 10 ** -6],\
                          1.0)

DEG2ARCSEC = 3600.0
S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2)
MAX_CHANNELS = 9000
MAX_BW = 2000000.0 # kHz
SPEED_OF_LIGHT = 299792458.0

class SynStruct:
    """ A synthetic structure """
    def __init__(self, log, code, intens, spa_form, spe_form):
        self.log = log
        self.code = code
        self.intens = intens
        self.spa_form = spa_form
        self.spe_form = spe_form
        # arrays to construct the table
        self.arr_mol = []
        self.arr_chname = []
        self.arr_rest_freq = []
        self.arr_obs_freq = []
        self.arr_fwhm = []
        self.arr_temp = []

    def setTemplate(self, template):
        self.template = template

    def addTransition(self, mol, chname, rest_freq, obs_freq, fwhm, temp):
        self.arr_mol.append(mol)
        self.arr_chname.append(chname)
        self.arr_rest_freq.append(rest_freq)
        self.arr_obs_freq.append(obs_freq)
        self.arr_fwhm.append(fwhm)
        self.arr_temp.append(temp)

    def getImageHDU(self):
        hdu = fits.PrimaryHDU()
        hdu.data = self.template;
        return hdu

    def getTableHDU(self):
        tbhdu = fits.new_table(fits.ColDefs([fits.Column(name='mol',\
                format='20A', array=self.arr_mol), fits.Column(name='chname',\
                format='40A', array=self.arr_chname), fits.Column(name='rest_freq',\
                format='D', array=self.arr_rest_freq), fits.Column(name='obs_freq',\
                format='D', array=self.arr_obs_freq),fits.Column(name='fwhm',\
                format='D', array=self.arr_fwhm), fits.Column(name='temp',\
                format='D', array=self.arr_temp)]))
        return tbhdu
    
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
        log.write('  -> Spatial Center (deg): ra=' + str(spec.x_center)\
                  + ' dec=' + str(spec.y_center) + '\n')
        fact = spec.ang_fov / (DEG2ARCSEC * spec.ang_res)
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
        # self.data = zeros((len(self.x_axis), len(self.y_axis), len(self.v_axis)))
        self.data = zeros((len(self.v_axis), len(self.x_axis), len(self.y_axis)))


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
        # self.im.set_array(self.data[:, :, j])
        self.im.set_array(self.data[j,:,:])
        return self.im,


    def animate(self, inte, rep=True):
        """ Animate the cube.
                inte       : time interval between frames
                 rep[=True] : boolean to repeat the animation
          """
        fig = plt.figure()
        # self.im = plt.imshow(self.data[:, :, 0], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
        #                     extent=(self.x_border[0], self.x_border[1], self.y_border[0], self.y_border[1]))
        self.im = plt.imshow(self.data[0,:,:], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
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
        # self.data += np.random.random((len(self.x_axis), len(self.y_axis), len(self.v_axis))) * noise
        self.data += np.random.random((len(self.v_axis), len(self.x_axis), len(self.y_axis))) * noise


    def getSpectrum(self, x, y):
        """ Return an spectrum in x, y """
        xi = int(round((x - self.x_axis[0]) / self.spec.ang_res))
        yi = int(round((y - self.y_axis[0]) / self.spec.ang_res))
        return self.data[:, xi, yi]


    def getCubeHDU(self):
        """ Write the final FITS file in filename """
        prihdr = fits.Header()
        prihdr['AUTHOR'] = 'Astronomical SYnthetic Data Observatory'
        prihdr['COMMENT'] = "Here's some commentary about this FITS file."
        prihdr['SIMPLE'] = True
        # prihdr['BITPIX'] = 8
        prihdr['NAXIS'] = 3
        hdu = fits.PrimaryHDU(header=prihdr)
        hdu.data = self.data
        return hdu
        # self.hdu.writeto(filename, clobber=True)

    def saveFits(self, sources, filename):
        hdulist = fits.HDUList([self.getCubeHDU()])
        for src in sources:
            for struct in sources[src].structs:
                hdulist.append(struct.getImageHDU())
                hdulist.append(struct.getTableHDU())
        hdulist.writeto(filename)

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
        self.structs.append(SynStruct(log, code, intens, spa_form, spe_form))
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
            struct.setTemplate(self.genSurface(struct.spa_form, cube))
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
                    freq = (1 + self.rad_vel*1000.0/SPEED_OF_LIGHT)*mlin[i]['frequency']
                    log.write('E: ' + mlin[i]['chemicalname'] + ' (' + mlin[i]['molformula'] + ') at ' + str(freq) + ' Mhz, T = exp(-|' \
                              + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(rinte) + ' = ' + str(temp) + ' K\n')
                    window = cube.freqWindow(freq / 1000.0, fwhm)
                    struct.addTransition(mol, mlin[i]['chemicalname'], mlin[i]['frequency'], freq, fwhm, temp)
                    log.write('Window:' + str(window) + '\n')
                    sigma = fwhm / S_FACTOR 
                    for idx in range(window[0], window[1]):
                        # cube.data[:, :, idx] += struct.template * temp * exp(
                        #    (-0.5 * (cube.v_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape)))
                        cube.data[idx,:,:] += struct.template * temp * exp(
                            (-0.5 * (cube.v_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape)))


    def loadLines(self, band, v_init, v_end, rad_vel):
        # TODO: Read from a database using SQLINE (SS Group)
        v_init_corr = (1 + rad_vel*1000.0/SPEED_OF_LIGHT)*v_init
        v_end_corr = (1 + rad_vel*1000.0/SPEED_OF_LIGHT)*v_end
        location = './votables/band' + band + '.xml'
        tbl = parse_single_table(location).to_table()
        print type(tbl)


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
        log.write('   * Adding Noise... \n')
        cube.addNoise() 
        cube.saveFits(self.sources, filename)
        return cube

    def removeSource(self, log, name):
        # TODO: must be implemented by VU group
        return 0

    def genImage(self, log, name, x_center, y_center, ang_res, ang_fov, filename):
        # TODO: must be implemented by the IS group
        return 0

# def loadLines(band):
        # TODO: Read from a database using SQLITE (SS Group)
        # v_init_corr=(1 + rad_vel*1000.0/SPEED_OF_LIGHT)*v_init
        # v_end_corr=(1 + rad_vel*1000.0/SPEED_OF_LIGHT)*v_end
#        location = './votables/band' + band + '.xml'
#        tbl = parse_single_table(location)
#        if isinstance(tbl,pTable):
        # tbl.array contiene los datos
        # tbl.field contiene la metadata
            # db = DataBase()
            # db.loadFields(tbl.fields)
            # db.genTable()
#            c = 0
#            data = tbl.array
#            datas = data._data
#            print type(datas)
#            for i in datas:
#                print i
#            print "hola"
#
#        return tbl



# loadLines("lite")
