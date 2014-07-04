
window = cube.freqWindow(freq / 1000.0, fwhm)
log.write('[W:' + str(window) + ']\n')
sigma = fwhm / S_FACTOR
distro = list()
for idx in range(window[0], window[1]):
distro.append(
np.exp((-0.5 * (cube.freq_axis[idx] - freq / 1000.0) ** 2) / (sigma ** (2 * shape))))
distro = distro / max(distro)

