import numpy as np
from astropy.io import fits

DEG2ARCSEC = 3600.0
MAX_CHANNELS = 9000
MAX_BW = 2000000.0 # kHz

class Cube:
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
                 Given a freq and a fwhm compute returns the range of channels that are affected
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
