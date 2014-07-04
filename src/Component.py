from astropy.io import fits
import random
import numpy as np
import math
import sqlite3 as lite

class Component:
    """Abstract Component Model"""
    def __init__(self,log,alpha,delta,z_base):
        self.log=log
        self.alpha=alpha
        self.delta=delta
        self.z_base=z_base

    def project(self,cube):
        pass

    def info(self):
        pass

class IMCM(Component):
    """ Interstellar Molecular Cloud Model """
    def __init__(self, log, alpha, delta,z_base, temp, mol_list , temp, spa_form, spe_form, z_grad):
        super(log,alpha,delta,z_base)
        self.spa_form = spa_form
        self.spe_form = spe_form
        self.z_grad=z_grad
        self.temp=temp
        self.intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(conf.base_abun[0], conf.base_abun[1])
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
                abun += conf.base_CO
            for iso in conf.iso_abun:
                if iso in mol:
                    abun *= conf.iso_abun[iso]
            self.intens[mol] = abun

    def changeIntensities(self,itens):
        '''User defined dictionary in the form {molecule: intensity}'''
        self.intens=intens;

    def info(self):
#TODD...
        pass

    def project(self, cube):
        arr_mol = []
        arr_chname = []
        arr_rest_freq = []
        arr_obs_freq = []
        arr_fwhm = []
        arr_temp = []
        self.log.write('   * Generating spatial form\n') # TODO More info
        T,Tbord=genSurface(spa_form,cube.alpha_axis,cube.delta_axis);
        self.log.write('   * Generating line form\n') #TODO More info
        L,Lbord=genLine(spe_form,cube.freq_axis);
        self.log.write('   * Loading and correcting lines with z=' + str(self.z_base) + ')\n')
        db = lite.connect('db/lines.db')
        freq_init_corr = cube.freq_border[0] * 1000.0 / (1 + self.z_base)
        freq_end_corr = cube.freq_border[1] * 1000.0 / (1 + self.z_base)
        for mol in self.intens:
            log.write("SQL SENTENCE:\n")
            select = "SELECT * FROM Lines WHERE SPECIES like '" + mol + "' AND FREQ > " + str(freq_init_corr) + " AND FREQ < " + str(freq_end_corr)
            log.write(select + '\n')
            resp = db.execute(select)
            linlist = resp.fetchall()
            rinte = inten_values[0]
            for j in range(len(inten_group)): # TODO baaad python...  
                if mol in inten_group[j]:
                    rinte = inten_values[j]
            rinte = random.uniform(rinte[0], rinte[1])
            for lin in linlist:
                trans_temp = lin[5]
                temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                freq = (1 + z) * lin[3]
                log.write('E: ' + str(lin[2]) + ' (' + str(lin[1]) + ') at ' + str(freq) + ' Mhz, T = exp(-|' \
                          + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(
                    rinte) + ' = ' + str(temp) + ' K  ')
                #window = cube.freqWindow(freq / 1000.0, fwhm)
                #log.write('[W:' + str(window) + ']\n')
                #sigma = fwhm / S_FACTOR
                #distro = list()
                #for idx in range(window[0], window[1]):
                #    distro.append(
                #        np.exp((-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
                #distro = distro / max(distro)
                for idx in range(Lbord[0], Lbord[1]):
                    cube.data[idx] = cube.data[idx] + struct.template * temp * S[idx - window[0]]
                arr_mol.append(mol)
                arr_temp.append(temp)
                arr_chname.append(chname)
                arr_rest_freq.append(rest_freq)
                arr_obs_freq.append(obs_freq)
                arr_fwhm.append(fwhm)
        hdu = fits.PrimaryHDU()
        hdu.data = T;
        tbhdu = fits.new_table(fits.ColDefs([fits.Column(name='mol',\
                format='20A', array=self.arr_mol), fits.Column(name='chname',\
                format='40A', array=self.arr_chname), fits.Column(name='rest_freq',\
                format='D', array=self.arr_rest_freq), fits.Column(name='obs_freq',\
                format='D', array=self.arr_obs_freq),fits.Column(name='fwhm',\
                format='D', array=self.arr_fwhm), fits.Column(name='temp',\
                format='D', array=self.arr_temp)]))
        
        return hdu,tbhdu
