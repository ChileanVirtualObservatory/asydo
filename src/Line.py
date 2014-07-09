import math
import numpy as np


S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2)

#def freqWindow(freq, fwhm):
#    """ Frequency window.
#             Given a freq and a fwhm compute returns the range of channels that are affected
#        """
#    factor = 2.0;
#    ldiff = freq - self.freq_border[0] - factor * fwhm;
#    udiff = freq - self.freq_border[0] + factor * fwhm;
#    l_w = 0
#    if (ldiff > 0):
#        l_w = ldiff * 1000000 / self.spec.spe_res
#    l_u = udiff * 1000000 / self.spec.spe_res
#    if l_u > self.channels:
#        l_u = self.channels - 1
#    return (int(l_w), int(l_u))

def phi(x,mu,sigma):
    'Cumulative distribution function for the standard normal distribution'
    return (1.0 + erf( (x - mu) /(sigma* sqrt(2.0)))) / 2.0

def genLine(spe_form,freq,freq_axis):
    #TODO: change broadening depending of the frequency
    fwhm=spe_form[1]
    shape=spe_form[2]
    window = [0,len(freq_axis)]
#cube.freqWindow(freq / 1000.0, fwhm)
    sigma = fwhm / S_FACTOR
    distro = list()
    for idx in range(window[0], window[1]):
        #distro.append(np.exp((-0.5 * (freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
        distro.append(np.exp((-0.5 * (freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
    distro = distro / sum(distro)
    return distro,window

