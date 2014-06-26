from astropy.io import fits

class Structure:
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
