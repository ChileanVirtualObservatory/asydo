import random
import numpy as np
import pylab as plt
import matplotlib.animation as animation
from astropy.io import fits
from scipy import pi, sqrt, exp
import math
from . import db
from scipy.special import erf



# ## Helper constants ###
SPEED_OF_LIGHT = 299792458.0
S_FACTOR = 2.354820045031  # sqrt(8*ln2)
KILO = 1000
DEG2ARCSEC = 3600.0

# ## ALMA Specific Constants ###
MAX_CHANNELS = 9000
MAX_BW = 2000.0  # MHz
ALMA_bands = {'3': [88000, 116000], '4': [125000, 163000], '6': [211000, 275000], '7': [275000, 373000],
              '8': [385000, 500000], '9': [602000, 720000]}
ALMA_noises = {'3': 0.01, '4': 0.012, '6': 0.02, '7': 0.04, '8': 0.08, '9': 0.16}

### IMC Specific Constants ###
inten_group = [('default'), ('COv=0'), ('13COv=0'), ('HCO+, HC3N, CS, C18O, CH3OH, N2H')]
inten_values = [[0.1, 2], [20, 60], [5, 20], [1, 10]]
default_iso_abundance = {'13C': 1.0 / 30, '18O': 1.0 / 60, '17O': 1.0 / 120, '34S': 1.0 / 30, '33S': 1.0 / 120,
                         '13N': 1.0 / 30, 'D': 1.0 / 30}

### Toolbox functions ###



def fwhm2sigma(freq,fwhm):
    """
    Compute the sigma in MHz given a frequency in MHz and a fwhm in km/s
    """
    sigma = (fwhm * 1000 / S_FACTOR) * (freq / SPEED_OF_LIGHT)
    return sigma

def freq_correct(freq,rv):
    freq_new = math.sqrt((1 + rv * KILO / SPEED_OF_LIGHT) / (1 - rv * KILO / SPEED_OF_LIGHT)) * freq
    return freq_new

def freq_window(freq, factor, axis):
    """
    Compute a window centered at freq within an axis.

    freq_window() returns a tuple with the lower and upper indices in the axis for a
    frequency window of freq +- factor. The size is at most 2*factor, but
    is limited by the axis borders. It returns (0,0) or (end,end) if the window is
    out of the axis by the left or right respectively (i.e., end = len(axis)).
    Equispaced axes are assumed.
    @rtype : tuple
    @param freq: center frequency
    @param factor: the window factor (half window size)
    @param axis: a equispaced array of elements
    """
    dlta = axis[1] - axis[0]
    ini = int(round((freq - factor - axis[0]) / dlta))
    end = int(round((freq + factor - axis[0]) / dlta))
    if ini < 0:
        ini = 0
    if end > len(axis):
        end = len(axis)
    return ini, end


def gen_line(spe_form, freq, axis):
    """
    Returns a spectral line distribution and its application window.

    gen_line() generates a normalized distribution (with compact support)
    of a spectral line centered at freq, within an axis. The spe_form
    parameter is a tuple where the first element is the name of the
    distribution, and the next elements are their parameters.

    The available distributions are:
    - ("skew",fwhm,alpha) : a skew-normal distribution, where fwhm is
              the full width at half maximum, and alpha is the curtosis
              parameter. If alpha = 0, it degenerates to a Gaussian
              distribution, if alpha < 0, the distribution is left-biased
              and alpha > 0 means a right bias.
    - TODO: implement more spectral distributions

    The function returns the tuple (distro,window), where distro is the
    spectral distribution (array), and window is a tuple of upper and
    lower indices (see freq_window function). The window factor is 3*sigma
    to include relatively large tails.
    @rtype : tuple
    @param spe_form: a spatial form (see doc)
    @param freq: central frequency
    @param axis: a equispaced array of elements
    """

    def pdf(x):
        return 1 / sqrt(2 * pi) * exp(-x ** 2 / 2)

    def cdf(x):
        return (1 + erf(x / sqrt(2))) / 2

    def skew(x, ep=0, wp=1, ap=0):
        t = (x - ep) / wp
        return 2 / wp * pdf(t) * cdf(ap * t)

    fwhm = spe_form[1]
    a = spe_form[2]

    sigma = fwhm2sigma(freq,fwhm)
    factor = 3 * sigma
    window = freq_window(freq, factor, axis)
    if window[0] > window[1]:
        return False, window

    d = a / sqrt(1.0 + a ** 2)
    w = sigma / sqrt(1.0 - (2.0 / pi) * d ** 2)
    e = freq - w * d * sqrt(2.0 / pi)
    distro = skew(axis[window[0]:window[1] + 1], e, w, a)
    ss = sum(distro)
    if ss != 0:
        distro = distro / ss
    #distro=np.sqrt(2*np.pi)*sigma*distro
    return distro, window


def spatial_window(alpha, delta, spx, spy, alpha_axis, delta_axis):
    """
    Compute a 2D window centered at (alpha,delta) within alpha_axis and delta_axis.

    spatial_window() returns a tuple (alpha_window,delta_window), where each element
    are the lower and upper indices for each axis. The window correspond to a window
    of freq +- spx or freq +- spy depending on the axis. This window is limited by
    the axis borders. Equispaced axes are assumed.
    @rtype : tuple
    @param alpha: RA center
    @param delta: DEC center
    @param spx: Semiaxis X
    @param spy: Semiaxis Y
    @param alpha_axis: a equispaced array of elements
    @param delta_axis: a equispaced array of elements
    """
    dlta_x = alpha_axis[1] - alpha_axis[0]
    dlta_y = delta_axis[1] - delta_axis[0]
    xbord = [int(round((alpha - spx - alpha_axis[0]) / dlta_x)), int(round((alpha + spx - alpha_axis[0]) / dlta_x))]
    ybord = [int(round((delta - spy - delta_axis[0]) / dlta_y)), int(round((delta + spy - delta_axis[0]) / dlta_y))]
    if xbord[0] < 0:
        xbord[0] = 0
    if ybord[0] < 0:
        ybord[0] = 0
    if xbord[1] > len(alpha_axis):
        xbord[1] = len(alpha_axis)
    if ybord[1] > len(delta_axis):
        ybord[1] = len(delta_axis)
    return ybord, xbord


def gen_surface(form, alpha, delta, alpha_axis, delta_axis):
    """
    Returns a spatial distribution and its application window.

    gen_surface() generates a standarized surface distribution
    (with compact support) centered at (alpha,delta), within the axes.
    The spa_form parameter is a tuple where the first element is the
    name of the distribution, and the next elements are their parameters.

    The available distributions are:
    - ("normal",sx,sy,theta) : a 2D Gaussian (normal) distribution,
                with diagonal covariance [[sx,0],[0,sy], rotated by
                theta (in radians!).
    - ("exp",sx,sy,theta) : exponential distribution, with semiaxes
                sx and sy, rotated by theta (in radians!).
    - TODO: implement more spatial distributions

    The function returns the tuple (distro,ybord,xbord), where distro is
    the spatial distribution (matrix), and ybord and xbord are the tuples
    of upper and lower indices for the axes (see spatial_window function).
    If the window is spurious (out of the axes or too small) this function
    returns False.
    """
    stype = form[0]
    sx = form[1] / DEG2ARCSEC
    sy = form[2] / DEG2ARCSEC
    theta = form[3]
    spx = abs(3 * sx * math.cos(theta)) + abs(3 * sy * math.sin(theta))
    spy = abs(3 * sx * math.sin(theta)) + abs(3 * sy * math.cos(theta))
    ybord, xbord = spatial_window(alpha, delta, spx, spy, alpha_axis, delta_axis)
    if xbord[0] > xbord[1] or ybord[0] > ybord[1]:
        return False, [ybord, xbord]
    alpha_axis = alpha_axis[xbord[0]:xbord[1] + 1]
    delta_axis = delta_axis[ybord[0]:ybord[1] + 1]
    alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse=False, indexing='xy')
    Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
    Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
    XX = (Xc) * math.cos(-theta) - (Yc) * math.sin(-theta)
    YY = (Xc) * math.sin(-theta) + (Yc) * math.cos(-theta)
    if stype == 'normal':
        u = (XX / sx) ** 2 + (YY / sy) ** 2
        sol = sx * sy * np.exp(-u / 2) / (2 * math.pi)
    elif stype == 'exp':
        u = sqrt((XX / sx) ** 2 + (YY / sy) ** 2)
        sol = sx * sy * np.exp(-u / sqrt(2))
    else:
        print('!!! ERROR: No such surface type')
        return False, [ybord, xbord]
    mm = max(sol)
    if mm != 0:
        sol = sol / mm
    res = np.reshape(sol, (len(delta_axis), len(alpha_axis)))
    return res, (ybord, xbord)


def gen_gradient(form, alpha, delta, alpha_axis, delta_axis, ybord, xbord):
    """
    Returns a matrix with gradient values within the axes, restricted to y_bord and x_bord.

    gen_gradient() generate a gradient matrix centered at (alpha,delta), for the alpha_axis
    and delta_axis, restricted by the windows ybord and xbord. The form parameter is a tuple
    where the first element is the gradient function, and the next elements are their parameters.
    The available gradient functions are:
    - ("linear", theta, m) : linear gradient, rotated by theta, and intensity m (km/s/arcsec)
    """
    np.set_printoptions(threshold=np.nan)
    gtype = form[0]
    theta = form[1]
    km_sarcs = form[2]
    alpha_axis = alpha_axis[xbord[0]:xbord[1]]
    delta_axis = delta_axis[ybord[0]:ybord[1]]
    alpha_mesh, delta_mesh = np.meshgrid(alpha_axis, delta_axis, sparse=False, indexing='xy')
    Xc = alpha_mesh.flatten() - alpha * np.ones(len(alpha_axis) * len(delta_axis))
    Yc = delta_mesh.flatten() - delta * np.ones(len(alpha_axis) * len(delta_axis))
    XX = Xc * math.cos(-theta) - (Yc) * math.sin(-theta);
    res = np.reshape(km_sarcs * XX, (len(delta_axis), len(alpha_axis)))
    return res


### Core ASYDO Classes ###

class Universe:
    """
    A synthetic universe where to put synthetic objects.
    """

    def __init__(self, log):
        """
        The log parameter is an opened file descriptor for logging purposes
        """
        self.log = log
        self.sources = dict()

    def create_source(self, name, alpha, delta):
        """
        A source needs a name and a spatial position (alpha,delta).
        """
        self.sources[name] = Source(self.log, name, alpha, delta)

    def add_component(self, source_name, model):
        """
        To add a component a Component object must be instantiated (model), and added to
        a source called source_name.
        """
        self.sources[source_name].add_component(model)

    def gen_cube(self, name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw):
        """
        Returns a SpectralCube object where all the sources within the FOV and BW are projected.

        This function needs the following parameters:
        - name    : name of the cube
        - alpha   : right-ascension center
        - delta   : declination center
        - freq    : spectral center (frequency)
        - ang_res : angular resolution
        - ang_fov : angular field of view
        - spe_res : spectral resolution
        - spe_bw  : spectral bandwidth
        """
        cube = SpectralCube(self.log, name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw)
        for src in self.sources:
            self.log.write('*** Source: ' + src + '\n')
            self.sources[src].project(cube)
        return cube

    def save_cube(self, cube, filename):
        """
        Wrapper function that saves a cube into a FITS (filename).
        """
        self.log.write('   -++ Saving FITS: ' + filename + '\n')
        cube.save_fits(self.sources, filename)

    def remove_source(self, name):
        """
        Deletes a source and its components.
        """
        self.log.write('Removing source ' + name)
        return self.sources.remove(name)


class Source:
    """
    A generic source of electromagnetic waves with several components.
    """

    def __init__(self, log, name, alpha, delta):
        """ Parameters:
               * log: logging descriptor
               * name: a name of the source
               * alpha: right ascension
               * delta: declination
        """
        self.log = log
        log.write('+++ Source \'' + name + '\' added\n')
        self.alpha = alpha
        self.delta = delta
        self.name = name
        self.comp = list()

    def add_component(self, model):
        """ Defines a new component from a model.
        """
        code = self.name + '-c' + str(len(self.comp) + 1)  #+ '-r' + str(self.alpha) +'-d'+str(self.delta)
        self.comp.append(model)
        model.register(code, self.alpha, self.delta)
        self.log.write(' |- Component added to \'' + self.name + '\' (' + code + ')\n')
        self.log.write(' ---+ Model: ' + model.info() + '\n')

    def project(self, cube):
        """
        Projects all components in the source to a cube.
        """
        for component in self.comp:
            self.log.write('  |- Projecting ' + component.comp_name + '\n')
            component.project(cube);


class SpectralCube:
    """
    A synthetic spectral cube.
    """

    def __init__(self, log, name, alpha, delta, freq, ang_res, ang_fov, spe_res, spe_bw, band_freq=ALMA_bands,
                 band_noises=ALMA_noises):
        """
        Obligatory Parameters:
        - log	  : descriptor of a log file
        - name    : name of the cube
        - alpha   : right-ascension center
        - delta   : declination center
        - freq    : spectral center (frequency)
        - ang_res : angular resolution
        - ang_fov : angular field of view
        - spe_res : spectral resolution
        - spe_bw  : spectral bandwidth

        Optional Parameters:
        - band_freq   : a diccionary of frequency ranges for the bands (key = band_name, value = (lower,upper))
        - band_noises : a dictionary of noise levels for each band (key = band_name, value = noise)
        """
        self.name = name
        self.alpha = alpha
        self.delta = delta
        self.freq = freq
        self.ang_res = ang_res
        self.ang_fov = ang_fov
        self.spe_res = spe_res
        self.spe_bw = spe_bw
        log.write('[*] Generating cube ' + name + '\n')
        log.write('  |- Angular Coordinates (deg): ra=' + str(alpha) + ' dec=' + str(delta) + '\n')
        fact = ang_fov / DEG2ARCSEC
        self.alpha_border = [alpha - fact / 2, alpha + fact / 2]
        self.delta_border = [delta - fact / 2, delta + fact / 2]
        self.alpha_axis = np.linspace(self.alpha_border[0], self.alpha_border[1], int(ang_fov / ang_res))
        self.delta_axis = np.linspace(self.delta_border[0], self.delta_border[1], int(ang_fov / ang_res))
        if alpha > 90 or alpha < -90:
            raise Exception('!!! ERROR: invalid coordinate: ra=' + alpha)
        if delta > 90 or delta < -90:
            raise Exception('!!! ERROR: invalid coordinate: dec=' + delta)
        log.write('  |- FOV (arcsec): ra=' + str(self.alpha_border) + ' dec=' + str(self.delta_border) + '\n')
        self.freq_border = [freq - spe_bw / 2.0, freq + spe_bw / 2.0]
        if spe_bw > MAX_BW:
            log.write('!!! WARNING: max ALMA bandwidth exceeded\n')
        self.channels = round(spe_bw / spe_res)
        if self.channels > MAX_CHANNELS:
            log.write('!!! WARNING: max ALMA channels exceeded\n')
        self.freq_axis = np.linspace(self.freq_border[0], self.freq_border[1], self.channels)
        log.write('  |- Spectral (MHz): center=' + str(freq) + ' bandwidth=' + str(self.freq_border) + '\n')
        log.write('  |- Cube size: ' + str(len(self.alpha_axis)) + ' x ' + str(len(self.delta_axis)) + ' x ' + str(
            len(self.freq_axis)) + ' \n')
        self.band = 'NO_BAND'
        for bnd in band_freq:
            freqs = band_freq[bnd]
            if self.freq_border[0] >= freqs[0] and self.freq_border[1] <= freqs[1]:
                self.band = bnd
                log.write('  |- Band: ' + bnd + '\n')
        if self.band == 'NO_BAND':
            log.write('!!! WARNING: not in a valid ALMA band\n')
        if self.band == 'NO_BAND':
            self.noise = 0.0001
        else:
            self.noise = band_noises[self.band]
        self.data = (
                        np.random.random(
                            (len(self.freq_axis), len(self.delta_axis), len(self.alpha_axis))) - 0.5 * np.ones(
                            (len(self.freq_axis), len(self.delta_axis), len(self.alpha_axis)))) * 2 * self.noise
        self.hdulist = fits.HDUList([self._get_cube_HDU()])


    def get_spectrum(self, x, y):
        """ Returns the spectrum of a (x,y) position """
        xi = int(round((x - self.alpha_axis[0]) / (self.ang_res / DEG2ARCSEC)))
        yi = int(round((y - self.delta_axis[0]) / (self.ang_res / DEG2ARCSEC)))
        return self.data[:, yi, xi]


    def _get_cube_HDU(self):
        prihdr = fits.Header()
        prihdr['AUTHOR'] = 'Astronomical SYnthetic Data Observatory'
        prihdr['COMMENT'] = "Here's some commentary about this FITS file."
        prihdr['SIMPLE'] = True
        # prihdr['BITPIX'] = 8
        prihdr['NAXIS'] = 3
        prihdr['NAXIS1'] = len(self.alpha_axis)
        prihdr['NAXIS2'] = len(self.delta_axis)
        prihdr['NAXIS3'] = len(self.freq_axis)
        prihdr['CTYPE3'] = "FREQ"
        prihdr['CRVAL3'] = self.freq
        prihdr['CRPIX3'] = int(len(self.freq_axis)/2)
        prihdr['CDELT3'] = self.spe_res
        prihdr['RESTFRQ'] = self.freq
        prihdr['CUNIT3'] = "Hz"
        prihdr['BUNIT'] = "JY/BEAM"
        hdu = fits.PrimaryHDU(header=prihdr)
        hdu.data = self.data
        return hdu

    def _add_HDU(self, hdu):
        pass
        self.hdulist.append(hdu)

    def save_fits(self, sources, filename):
        """ Simple as that... saves the whole cube """
        self.hdulist.writeto(filename, clobber=True)

    def _updatefig(self, j):
        """ Animate helper function """
        self.im.set_array(self.data[j, :, :])
        return self.im,

    def animate(self, inte, rep=True):
        """ Simple animation of the cube.
            - inte       : time interval between frames
            - rep[=True] : boolean to repeat the animation
          """
        fig = plt.figure()
        self.im = plt.imshow(self.data[0, :, :], cmap=plt.get_cmap('jet'), vmin=self.data.min(), vmax=self.data.max(), \
                             extent=(
                                 self.alpha_border[0], self.alpha_border[1], self.delta_border[0],
                                 self.delta_border[1]))
        ani = animation.FuncAnimation(fig, self._updatefig, frames=range(len(self.freq_axis)), interval=inte, blit=True,
                                      repeat=rep)
        plt.show()


class Component:
    """Abstract component model"""

    def __init__(self, log, z_base=0.0):
        """ log: file descriptor for logging
            z_base[=0] : optional parameter to set base redshift (if not, please use set_radial_vel)
        """
        self.log = log
        self.z = z_base
        self.rv = SPEED_OF_LIGHT / KILO * ((self.z ** 2 + 2 * self.z) / (self.z ** 2 + 2 * self.z + 2))

    def set_radial_velocity(self, rvel):
        """Set radial velocity rvel in km/s"""
        self.rv = rvel
        self.z = math.sqrt((1 + self.rv * KILO / SPEED_OF_LIGHT) / (1 - self.rv * KILO / SPEED_OF_LIGHT)) - 1

    def info(self):
        """Print relevant information of the component"""
        return "(none)"

    def register(self, comp_name, alpha, delta):
        """Register the component name and angular position (alpha,delta)"""
        self.comp_name = comp_name
        self.alpha = alpha
        self.delta = delta

    def project(self, cube):
        """Project the component in the cube"""
        pass


class IMCM(Component):
    """ Interstellar Molecular Cloud Model """

    def __init__(self, log, dbpath, mol_list, temp, spa_form, spe_form, z_grad, z_base=0.0, abun_max=10 ** -5,
                 abun_min=10 ** -6, abun_CO=1.0, iso_abun=default_iso_abundance):
        Component.__init__(self, log, z_base)
        self.spa_form = spa_form
        self.spe_form = spe_form
        self.z_grad = z_grad
        self.dbpath = dbpath
        self.temp = temp
        self.intens = dict()
        for mol in mol_list.split(','):
            abun = random.uniform(abun_min, abun_max)
            if mol in ('COv=0', '13COv=0', 'C18O', 'C17O', '13C18O'):
                abun += abun_CO
            for iso in iso_abun:
                if iso in mol:
                    abun *= iso_abun[iso]
            self.intens[mol] = abun

    def change_intensities(self, intens):
        '''User defined dictionary in the form {molecule: intensity}'''
        self.intens = intens;

    def info(self):
        return "mol_list = " + str(self.intens.keys()) + " @ spa_form=" + str(self.spa_form) + ", spe_form=" + str(
            self.spe_form) + ", z=" + str(self.z) + ", grad=" + str(self.z_grad)


    def project(self, cube):
        arr_code = []
        arr_mol = []
        arr_chname = []
        arr_rest_freq = []
        arr_rad_vel = []
        arr_fwhm = []
        arr_temp = []
        self.log.write('   --+ Generating template image\n')  # TODO More info
        T, (ybord, xbord) = gen_surface(self.spa_form, self.alpha, self.delta, cube.alpha_axis, cube.delta_axis)
        if isinstance(T, bool):
            return
        G = gen_gradient(self.z_grad, self.alpha, self.delta, cube.alpha_axis, cube.delta_axis, ybord, xbord)
        self.log.write('   --+ Generating template line distribution\n')  #TODO More info
        self.log.write('   --+ Loading and correcting lines\n')
        dba = db.lineDB(self.dbpath)
        dba.connect()
        freq_init_corr = cube.freq_border[0] / (1 + self.z)
        freq_end_corr = cube.freq_border[1] / (1 + self.z)
        counter = 0
        used = False
        for mol in self.intens:
            # For each molecule specified in the dictionary
            # load its spectral lines

            linlist = dba.getSpeciesLines(mol, freq_init_corr,
                                          freq_end_corr)  # Selected spectral lines for this molecule
            rinte = inten_values[0]
            for j in range(len(inten_group)):  # TODO baaad python...
                if mol in inten_group[j]:
                    rinte = inten_values[j]
            rinte = random.uniform(rinte[0], rinte[1])

            for lin in linlist:
                counter += 1
                trans_temp = lin[5]
                temp = np.exp(-abs(trans_temp - self.temp) / self.temp) * rinte
                if temp < 3 * cube.noise:
                    continue
                freq = (1 + self.z) * lin[3]  # Catalogs must be in Mhz
                self.log.write('      |- Projecting ' + str(lin[2]) + ' (' + str(lin[1]) + ') around ' + str(
                    freq) + ' Mhz, at ' + str(temp) + ' K\n')
                for xp in range(xbord[0], xbord[1]):
                    for yp in range(ybord[0], ybord[1]):
                        freq = freq_correct(lin[3],self.rv + G[yp - ybord[0], xp - xbord[0]])
                        L, Lbord = gen_line(self.spe_form, freq, cube.freq_axis)
                        if isinstance(L, bool):
                            continue
                        cube.data[Lbord[0]:Lbord[1] + 1, yp, xp] = cube.data[Lbord[0]:Lbord[1] + 1, yp, xp] + T[yp - ybord[0], xp - xbord[0]] * temp* L
                        used = True

                arr_code.append(self.comp_name + '-r' + str(self.alpha) + '-d' + str(self.delta) + "-l" + str(counter))
                arr_mol.append(mol)
                arr_temp.append(temp)
                arr_chname.append(str(lin[2]))
                arr_rest_freq.append(str(lin[3]))
                arr_rad_vel.append(self.rv)
                arr_fwhm.append(self.spe_form[1])
        dba.disconnect()
        if not used:
            return
        hduT = fits.PrimaryHDU()
        hduT.data = T;
        hduG = fits.PrimaryHDU()
        hduG.data = G;
        table_columns = fits.ColDefs([
            fits.Column(name='line_code', format='60A', array=arr_code),
            fits.Column(name='mol', format='20A', array=arr_mol), \
            fits.Column(name='chname', format='40A', array=arr_chname), \
            fits.Column(name='rest_freq', format='D', array=arr_rest_freq), \
            fits.Column(name='rad_vel', format='D', array=arr_rad_vel), \
            fits.Column(name='fwhm', format='D', array=arr_fwhm), \
            fits.Column(name='temp', format='D', array=arr_temp)])
        tbhdu = fits.BinTableHDU.from_columns(table_columns)
        cube._add_HDU(hduT)
        cube._add_HDU(hduG)
        cube._add_HDU(tbhdu)
