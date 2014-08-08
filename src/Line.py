#from math import sqrt,log,pi,exp
from scipy import linspace
from scipy import pi,sqrt,exp
from scipy.special import erf
import numpy as np

SPEED_OF_LIGHT = 299792458.0
S_FACTOR = 2.3548200450309493820231386529193992754947713787716410 #sqrt(8*ln2)

def freqWindow(ini,end, freq_axis):
    """ Frequency window.
       """
    idx1=0
    idx2=len(freq_axis)
    for i in range(len(freq_axis)):
        if ini > freq_axis[i]:
            idx1=i
        if end < freq_axis[i]:
            idx2=i
    return (idx1,idx2)

def pdf(x):
    return 1/sqrt(2*pi) * exp(-x**2/2)

def cdf(x):
    return (1 + erf(x/sqrt(2))) / 2

def skew(x,e=0,w=1,a=0):
    t = (x-e) / w
    return 2 / w * pdf(t) * cdf(a*t)

#def phi(x,mu,sigma):
#    'Cumulative distribution function for the standard normal distribution'
#    return (1.0 + erf( (x - mu) /(sigma* sqrt(2.0)))) / 2.0

def genLine(spe_form,freq,freq_axis):
    #TODO: change broadening depending of the frequency
    fwhm=spe_form[1]
    a=spe_form[2]
    
    #window = [0,len(freq_axis)]
    sigma = (fwhm*1000 / S_FACTOR) * (freq/SPEED_OF_LIGHT)
    factor=3*sigma
    window=freqWindow(freq - factor,freq + factor,freq_axis)

    #distro = zeros(window[0],window[1])
    
    d = a/sqrt(1.0 + a**2)
    #for idx in range(window[0], window[1]):
    #    distro.append(np.exp((-0.5 * (freq_axis[idx] - freq) ** 2) / (sigma ** 2)))
    w = sigma/sqrt(1.0 - (2.0/pi)*d**2)
    e = freq - w*d*sqrt(2.0/pi)
    distro = skew(freq_axis[window[0]:window[1]],e,w,a)
    distro = distro / sum(distro)
    return distro,window

