from Component import *
from Surface import *
from Line import *
import sqlite3 as lite
import random
from astropy.io import fits
import numpy as np

import math

#TODO: parametrize this!
inten_group=[('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')]
inten_values=[[0.1, 2], [20, 60], [5, 20], [1, 10]]

default_iso_abundance={'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': 1.0/120, '34S': 1.0 / 30, '33S': 1.0 / 120,'13N': 1.0 / 30, 'D': 1.0 / 30}

class IMCM(Component):
    
    """ Interstellar Molecular Cloud Model """
    def __init__(self, log, mol_list, temp, spa_form, spe_form, z_grad,z_base=0.0, abun_max=10 ** -5, abun_min=10 ** -6, abun_CO=1.0, iso_abun=default_iso_abundance):
        Component.__init__(self,log,z_base)
        self.spa_form = spa_form
        self.spe_form = spe_form
        self.z_grad=z_grad
        self.temp=temp
        self.intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(abun_min, abun_max)
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
                abun += abun_CO
            for iso in iso_abun:
                if iso in mol:
                    abun *= iso_abun[iso]
            self.intens[mol] = abun

    def changeIntensities(self,intens):
        '''User defined dictionary in the form {molecule: intensity}'''
        self.intens=intens;

    def info(self):
        #TODO...
        return "z = "+str(self.z)

    def project(self, cube):
        arr_code = []
        arr_mol = []
        arr_chname = []
        arr_rest_freq = []
        arr_obs_freq = []
        arr_fwhm = []
        arr_temp = []
        self.log.write('   * Generating spatial form\n') # TODO More info
        T,Tbord=genSurface(self.spa_form,self.alpha,self.delta,cube.alpha_axis,cube.delta_axis)
        
        G=genGradient(self.z_grad,self.alpha,self.delta,cube.alpha_axis,cube.delta_axis,Tbord)

        self.log.write('   * Generating line form\n') #TODO More info
        self.log.write('   * Loading and correcting lines with z=' + str(self.z) + ')\n')
        db = lite.connect('db/lines.db')
        freq_init_corr = cube.freq_border[0] / (1 + self.z)
        freq_end_corr = cube.freq_border[1] / (1 + self.z)
        counter=0

        for mol in self.intens:
            # For each molecule specified in the dictionary
            # load its spectral lines

            self.log.write("SQL SENTENCE:\n")
            select = "SELECT * FROM Lines WHERE SPECIES like '" + mol + "' AND FREQ > " + str(freq_init_corr) + " AND FREQ < " + str(freq_end_corr)
            self.log.write(select + '\n')
            resp = db.execute(select)

            linlist = resp.fetchall()        # Selected spectral lines for this molecule

            rinte = inten_values[0]
            for j in range(len(inten_group)): # TODO baaad python...  
                if mol in inten_group[j]:
                    rinte = inten_values[j]
            rinte = random.uniform(rinte[0], rinte[1])

            for lin in linlist:
                counter+=1
                trans_temp = lin[5]
                temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                freq = (1 + self.z) * lin[3] # Catalogs must be in Mhz 
                L,Lbord=genLine(self.spe_form,freq,cube.freq_axis)
                self.log.write('E: ' + str(lin[2]) + ' (' + str(lin[1]) + ') at ' + str(freq) + ' Mhz, T = exp(-|' \
                          + str(trans_temp) + '-' + str(self.temp) + '|/' + str(self.temp) + ')*' + str(
                    rinte) + ' = ' + str(temp) + ' K  ' + '\n')
                for idx in range(Lbord[0], Lbord[1]):
                    cube.data[idx] = cube.data[idx] + T * temp * L[idx - Lbord[0]]
                arr_code.append(self.comp_name + '-r' + str(self.alpha) +'-d'+str(self.delta) + "-l" + str(counter))
                arr_mol.append(mol)
                arr_temp.append(temp)
                arr_chname.append(str(lin[2]))
                arr_rest_freq.append(str(lin[3]))
                arr_obs_freq.append(freq)
                arr_fwhm.append(self.spe_form[1])

        hdu = fits.PrimaryHDU()
        hdu.data = T;
        #TODO Add redshift

        #Add Metadada to the FIT
        tbhdu = fits.new_table(fits.ColDefs([
                     fits.Column(name='line_code',format='60A',array=arr_code), 
                     fits.Column(name='mol',format='20A', array=arr_mol),\
                     fits.Column(name='chname',format='40A', array=arr_chname),\
                     fits.Column(name='rest_freq',format='D', array=arr_rest_freq),\
                     fits.Column(name='obs_freq',format='D', array=arr_obs_freq),\
                     fits.Column(name='fwhm',format='D', array=arr_fwhm),\
                     fits.Column(name='temp',format='D', array=arr_temp)]))

        cube.addHDU(hdu)
        cube.addHDU(tbhdu)



