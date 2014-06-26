import random
import numpy as np
import math
import sqlite3 as lite
from Structure import *

SPEED_OF_LIGHT = 299792458.0
S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2) #TODO Preguntar Profe que es esto?


class Source:
    """A source of EM Waves"""

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
        self.structs.append(Structure(log, code, intens, spa_form, spe_form))
        log.write('Added to \'' + self.name + '\': molecules ' + str(mol_list) + '\n')
        log.write('  -> Spatial form  = ' + str(spa_form) + '\n')
        log.write('  -> Spectral form = ' + str(spe_form) + '\n')

    def genSurface(self, form, cube):
        "Create a gaussian surface over a mesh created by x and y axes"
        sx = form[1]
        sy = form[2]
        theta = form[3]
        r = 3 * math.sqrt(sx ** 2 + sy ** 2)
        alpha_mesh, delta_mesh = np.meshgrid(cube.alpha_axis, cube.delta_axis, sparse = False, indexing = 'xy')
        Xc = alpha_mesh.flatten() - self.alpha * np.ones(len(cube.alpha_axis) * len(cube.delta_axis))
        Yc = delta_mesh.flatten() - self.delta * np.ones(len(cube.alpha_axis) * len(cube.delta_axis))
        XX = (Xc) * math.cos(theta) - (Yc) * math.sin(theta);
        YY = (Xc) * math.sin(theta) + (Yc) * math.cos(theta);
        u = (XX / sx) ** 2 + (YY / sy) ** 2;
        sol = sx * sy * np.exp(-u / 2) / (2 * math.pi);
        res = np.transpose(np.reshape(sol, (len(cube.delta_axis), len(cube.delta_axis))))
        res = res / res.max()
        return res

    def emission(self, log, cube, inten_group, inten_values):
        log.write('Loading visible lines in band ' + cube.band + ' (rad_vel=' + str(self.rad_vel) + ')\n')
        db = lite.connect('db/lines.db')
        #lines = self.loadLines(cube.band, cube.freq_border[0], cube.freq_border[1], self.rad_vel)
        for struct in self.structs:
            struct.clear()
            log.write(' --> Struct Name: ' + struct.code + '\n')
            struct.setTemplate(self.genSurface(struct.spa_form, cube))
            fwhm = struct.spe_form[1]
            shape = struct.spe_form[2]
            freq_init_corr = cube.freq_border[0] * 1000.0 / (1 + self.rad_vel * 1000.0 / SPEED_OF_LIGHT)
            freq_end_corr = cube.freq_border[1] * 1000.0 / (1 + self.rad_vel * 1000.0 / SPEED_OF_LIGHT)

            for mol in struct.intens:
                log.write("SQL SENTENCE:\n")
                select = "SELECT * FROM Lines WHERE SPECIES like '" + mol + "' AND FREQ > " + str(
                    freq_init_corr) + " AND FREQ < " + str(freq_end_corr)
                log.write(select + '\n')
                resp = db.execute(select)
                linlist = resp.fetchall()
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
                    trans_temp = lin[5]
                    temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                    freq = (1 + self.rad_vel * 1000.0 / SPEED_OF_LIGHT) * lin[3]
                    log.write('E: ' + str(lin[2]) + ' (' + str(lin[1]) + ') at ' + str(freq) + ' Mhz, T = exp(-|' \
                              + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(
                        rinte) + ' = ' + str(temp) + ' K  ')
                    window = cube.freqWindow(freq / 1000.0, fwhm)
                    struct.addTransition(mol, str(lin[2]), str(lin[3]), freq, fwhm, temp)
                    log.write('[W:' + str(window) + ']\n')
                    sigma = fwhm / S_FACTOR
                    distro = list()
                    for idx in range(window[0], window[1]):
                        # cube.data[:, :, idx] += struct.template * temp * exp(
                        #    (-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape)))
                        distro.append(
                            np.exp((-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
                    distro = distro / max(distro)
                    #log.write(str(distro)+'\n')
                    for idx in range(window[0], window[1]):
                        cube.data[idx] = cube.data[idx] + struct.template * temp * distro[idx - window[0]]


#    def loadLines(self, band, freq_init, freq_end, rad_vel):
#        # TODO: Read from a database using SQLINE (SS Group)
#        freq_init_corr = (1 + rad_vel*1000.0/SPEED_OF_LIGHT)*freq_init
#        freq_end_corr = (1 + rad_vel*1000.0/SPEED_OF_LIGHT)*freq_end
#        location = './votables/band' + band + '.xml'
#        tbl = parse_single_table(location).to_table()
#        print type(tbl)


#        return tbl

import sys

log = sys.stdout

casa = Source(log,'P-33SO2', 0.0, 1.0, 150, 300.0)