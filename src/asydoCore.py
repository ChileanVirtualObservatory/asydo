import numpy as np
import math
from collections import namedtuple
import random
import string
from astropy.io import fits
import sqlite3 as lite
#from astropy.io.votable import parse_single_table
#from astropy.table import Table
#from astropy.io.votable.tree import Table as pTable
#from astropy.io.votable.tree import Field as pField
#import atpy
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
        self.clear()

    def clear(self):
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
        log.write('  -> Angular Coordinates (deg): ra=' + str(spec.alpha)\
                  + ' dec=' + str(spec.delta) + '\n')
        fact = spec.ang_fov / (DEG2ARCSEC * spec.ang_res)
        self.alpha_border = [spec.alpha - fact / 2, spec.alpha + fact / 2]
        self.delta_border = [spec.delta - fact / 2, spec.delta + fact / 2]
        self.alpha_axis = np.linspace(self.alpha_border[0], self.alpha_border[1], spec.ang_fov / spec.ang_res)
        self.delta_axis = np.linspace(self.delta_border[0], self.delta_border[1], spec.ang_fov / spec.ang_res)
        if spec.alpha > 90 or spec.alpha < -90:
            raise Exception('ERROR: invalid coordinate: ra=' + spec.alpha)
        if spec.delta > 90 or spec.delta < -90:
            raise Exception('ERROR: invalid coordinate: dec=' + spec.delta)
        log.write('  -> FOV (arcsec): ra=' + str(self.alpha_border) + ' dec=' + str(self.delta_border) + '\n')
        self.freq_border = [spec.freq - spec.spe_bw / (2.0 * 1000000), spec.freq + spec.spe_bw / (2.0 * 1000000)]
        if spec.spe_bw > MAX_BW:
            log.write('WARNING: max ALMA bandwidth exceeded\n')
        self.channels = round(spec.spe_bw / spec.spe_res)
        if self.channels > MAX_CHANNELS:
            log.write('WARNING: max ALMA channels exceeded\n')
        self.freq_axis = np.linspace(self.freq_border[0], self.freq_border[1], self.channels)
        log.write('  -> Spectral (GHz): center=' + str(spec.freq) + ' bandwidth=' + str(self.freq_border) + '\n')
        log.write('  -> Cube size: ' + str(len(self.alpha_axis)) + ' x ' + str(len(self.delta_axis)) + ' x ' + str(
            len(self.freq_axis)) + ' \n')
        self.band = 'NO_BAND'
        for bnd in conf.band_freq:
            freqs = conf.band_freq[bnd]
            if self.freq_border[0] >= freqs[0] and self.freq_border[1] <= freqs[1]:
                self.band = bnd
                log.write('  -> Band: ' + bnd + '\n')
        if self.band == 'NO_BAND':
            log.write('WARNING: not in a valid ALMA band\n')
        # self.data = zeros((len(self.alpha_axis), len(self.delta_axis), len(self.freq_axis)))
        if self.band == 'NO_BAND':
            noise = 0.0001
        else:
            noise = self.conf.band_noise[self.band]
        self.data = np.random.random((len(self.freq_axis), len(self.alpha_axis), len(self.delta_axis))) * noise


    def freqWindow(self, freq, fwhm):
        """ Frequency window.
                 Given a freq and a fwhm compute returns the range of channels that are affectet
            """
        factor = 2.0;
        ldiff = freq - self.freq_border[0] - factor * fwhm;
        udiff = freq - self.freq_border[0] + factor * fwhm;
        l_w = 0
        if (ldiff > 0):
            l_w = ldiff * 1000000 / self.spec.spe_res
        l_u = udiff * 1000000 / self.spec.spe_res
        if l_u > self.channels:
            l_u = self.channels - 1
        return (int(l_w), int(l_u))


#    def _updatefig(self, j):
#        """ Animate helper function """
#        # self.im.set_array(self.data[:, :, j])
#        self.im.set_array(self.data[j,:,:])
#        return self.im,


#    def animate(self, inte, rep=True):
#        """ Animate the cube.
#                inte       : time interval between frames
#                 rep[=True] : boolean to repeat the animation
#          """
#        fig = plt.figure()
#        # self.im = plt.imshow(self.data[:, :, 0], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
#        #                     extent=(self.alpha_border[0], self.alpha_border[1], self.delta_border[0], self.delta_border[1]))
#        self.im = plt.imshow(self.data[0,:,:], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
#                             extent=(self.alpha_border[0], self.alpha_border[1], self.delta_border[0], self.delta_border[1]))
#        ani = animation.FuncAnimation(fig, self._updatefig, frames=range(len(self.freq_axis)), interval=inte, blit=True,
#                                      repeat=rep)
#        plt.show()




    def getSpectrum(self, x, y):
        """ Return an spectrum in x, y """
        xi = int(round((x - self.alpha_axis[0]) / self.spec.ang_res))
        yi = int(round((y - self.delta_axis[0]) / self.spec.ang_res))
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
        hdulist.writeto(filename, clobber=True)

class SynSource:
    """A source"""

    def __init__(self, log, name, alpha, delta, rad_vel, temp):
        log.write('Source \'' + name + '\' added \n')
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.rad_vel = rad_vel
        self.temp = temp
        self.structs = list()

    def addStruct(self, log, mol_list, spa_form, spe_form, conf):
        code = self.name + '-' + str(len(self.structs) + 1)
        intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(conf.base_abun[0], conf.base_abun[1])
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
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
        r = 3 * math.sqrt(sx ** 2 + sy ** 2)
        alpha_mesh, delta_mesh = np.meshgrid(cube.alpha_axis, cube.delta_axis, sparse=False, indexing='xy')
        Xc = alpha_mesh.flatten() - self.alpha * np.ones(len(cube.alpha_axis) * len(cube.delta_axis))
        Yc = delta_mesh.flatten() - self.delta * np.ones(len(cube.alpha_axis) * len(cube.delta_axis))
        XX = (Xc) * math.cos(theta) - (Yc) * math.sin(theta);
        YY = (Xc) * math.sin(theta) + (Yc) * math.cos(theta);
        u = (XX / sx) ** 2 + (YY / sy) ** 2;
        sol = sx * sy * np.exp(-u / 2) / (2 * math.pi);
        res = np.transpose(np.reshape(sol, (len(cube.delta_axis), len(cube.delta_axis))))
        res=res/res.max()
        return res

    def emission(self, log, cube, inten_group, inten_values):
        log.write('Loading visible lines in band ' + cube.band + ' (rad_vel=' + str(self.rad_vel) + ')\n')
        db=lite.connect('db/lines.db')
        #lines = self.loadLines(cube.band, cube.freq_border[0], cube.freq_border[1], self.rad_vel)
        for struct in self.structs:
            struct.clear()
            log.write(' --> Struct Name: ' + struct.code + '\n')
            struct.setTemplate(self.genSurface(struct.spa_form, cube))
            fwhm = struct.spe_form[1]
            shape = struct.spe_form[2]
            freq_init_corr = cube.freq_border[0]*1000.0/(1 + self.rad_vel*1000.0/SPEED_OF_LIGHT)
            freq_end_corr = cube.freq_border[1]*1000.0/(1 + self.rad_vel*1000.0/SPEED_OF_LIGHT)

            for mol in struct.intens:
                log.write("SQL SENTENCE:\n")
                select="SELECT * FROM Lines WHERE SPECIES like '"+mol+"' AND FREQ > "+str(freq_init_corr)+" AND FREQ < "+str(freq_end_corr)
                log.write(select+'\n')
                resp=db.execute(select)
                linlist=resp.fetchall()
                rinte = inten_values[0]
                for j in range(len(inten_group)):
                    if mol in inten_group[j]:
                        rinte = inten_values[j]
                rinte = random.uniform(rinte[0], rinte[1])
                #for i in range(len(mlin)):
                #    try:
                #        trans_temp = float(mlin[i]['lowerstateenergyK'])
                #    except ValueError:
                #        log.write('WARNING: Strange transition temperature, printing line...\n')
                #        log.write('WARNING: ' + str(mlin[i]) + '\n')
                #        log.write('WARNING: Setting transition temperature to ' + str(self.temp) + '\n')
                #        trans_temp = self.temp
                for lin in linlist:
                    trans_temp=lin[5]
                    temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                    freq = (1 + self.rad_vel*1000.0/SPEED_OF_LIGHT)*lin[3]
                    log.write('E: ' + str(lin[2]) + ' (' + str(lin[1]) + ') at ' + str(freq) + ' Mhz, T = exp(-|' \
                              + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(rinte) + ' = ' + str(temp) + ' K  ')
                    window = cube.freqWindow(freq / 1000.0, fwhm)
                    struct.addTransition(mol,  str(lin[2]) , str(lin[3]) , freq, fwhm, temp)
                    log.write('[W:' + str(window) + ']\n')
                    sigma = fwhm / S_FACTOR
                    distro=list()
                    for idx in range(window[0], window[1]):
                        # cube.data[:, :, idx] += struct.template * temp * exp(
                        #    (-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape)))
                        distro.append(np.exp((-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
                    distro=distro/max(distro)
                    #log.write(str(distro)+'\n')
                    for idx in range(window[0], window[1]):
                        cube.data[idx] = cube.data[idx] + struct.template * temp * distro[idx-window[0]]


#    def loadLines(self, band, freq_init, freq_end, rad_vel):
#        # TODO: Read from a database using SQLINE (SS Group)
#        freq_init_corr = (1 + rad_vel*1000.0/SPEED_OF_LIGHT)*freq_init
#        freq_end_corr = (1 + rad_vel*1000.0/SPEED_OF_LIGHT)*freq_end
#        location = './votables/band' + band + '.xml'
#        tbl = parse_single_table(location).to_table()
#        print type(tbl)


#        return tbl


class SynUniverse:
    """Synthetic Cube Generation Class"""

    def __init__(self, log):
        self.conf = defaultUniverse
        self.sources = dict()

    def setConf(self, log, univ):
        self.conf = univ

    def addSource(self, log, name, alpha , delta, rad_vel, temp):
        self.sources[name] = SynSource(log, name, alpha, delta, rad_vel, temp)

    def addStruct(self, log, name, mol_list, spa_form, spe_form):
        self.sources[name].addStruct(log, mol_list, spa_form, spe_form, self.conf)

    def genCube(self, log, name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw, filename):
        cube = SynCube(log, name, self.conf, CubeSpec(alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw))
        for src in self.sources: 
            log.write('   * Source: ' + src + '\n')
            self.sources[src].emission(log, cube, self.conf.inten_group, self.conf.inten_values)
        log.write('   * Adding Noise... \n')
        cube.saveFits(self.sources, filename)
        return cube

    def removeSource(self, log, name):
        # TODO: must be implemented by VU group
        return 0

    def genImage(self, log, name, alpha, delta, ang_res, ang_fov, filename):
        # TODO: must be implemented by the IS group
        return 0

# def loadLines(band):
        # TODO: Read from a database using SQLITE (SS Group)
        # freq_init_corr=(1 + rad_vel*1000.0/SPEED_OF_LIGHT)*freq_init
        # freq_end_corr=(1 + rad_vel*1000.0/SPEED_OF_LIGHT)*freq_end
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



